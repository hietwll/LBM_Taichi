import taichi as ti
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

ti.init(arch=ti.gpu)

nx = 256
ny = 256
dx = 1.0
dy = 1.0
rho0 = 1.0
u_c = 0.1
re = 1000.0
lx = dx*float(nx-1)
ly = dy*float(ny-1)
niu = u_c*lx/re
tau = 3.0*niu+0.5


f_old = ti.Vector(9, dt=ti.f32, shape=(nx, ny))
f_new = ti.Vector(9, dt=ti.f32, shape=(nx, ny))

rho = ti.var(dt=ti.f32, shape=(nx, ny))
vel = ti.Vector(2, dt=ti.f32, shape=(nx, ny))

show_var = ti.var(dt=ti.f32, shape=(nx, ny))

w = ti.Vector([4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,
               1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0])

e = ti.Vector([[0,0], [1,0], [0,1], 
    [-1,0], [0, -1], [1,1], [-1,1], [-1,-1], [1,-1]])


# compute the distribution function at equilibrium state
@ti.func
def f_eq(ek0, ek1, wk, rho, u, v):
    eu = ti.cast(ek0,ti.f32)*u + ti.cast(ek1,ti.f32)*v
    uv = u*u + v*v
    return wk*rho*(1.0+3.0*eu+4.5*eu**2-1.5*uv)


@ti.kernel
def init():
    for i,j in rho:
        vel[i,j][0] = 0.0
        vel[i,j][1] = 0.0
        rho[i,j] = rho0
        for k in ti.static(range(9)):
             f_new[i,j][k] = f_eq(e[k,0],e[k,1], w[k], rho[i,j], vel[i,j][0], vel[i,j][1])
             f_old[i,j][k] = f_new[i,j][k]

@ti.kernel             
def collide_and_stream():
    for i,j in ti.ndrange((1,nx-1),(1,ny-1)):
        for k in ti.static(range(9)):
            ip = i - e[k,0]
            jp = j - e[k,1]
            f_new[i,j][k] = (1.0-1.0/tau)*f_old[ip,jp][k] + \
            f_eq(e[k,0],e[k,1], w[k], rho[ip,jp], vel[ip,jp][0], vel[ip,jp][1])/tau

@ti.kernel
def update_macro_var():
    for i,j in ti.ndrange((1,nx-1),(1,ny-1)):
        rho[i,j]= 0.0
        vel[i,j][0] = 0.0
        vel[i,j][1] = 0.0
        for k in ti.static(range(9)):
            f_old[i,j][k] = f_new[i,j][k]
            rho[i,j] += f_new[i,j][k]
            vel[i,j][0] += (ti.cast(e[k,0],ti.f32)*f_new[i,j][k])
            vel[i,j][1] += (ti.cast(e[k,1],ti.f32)*f_new[i,j][k])

        vel[i,j][0] /= rho[i,j]
        vel[i,j][1] /= rho[i,j]


@ti.kernel
def get_show_var():
    for i,j in ti.ndrange(nx,ny):
        show_var[i,j] = vel[i,j][0]/u_c

@ti.kernel
def apply_bc():
    # left and right
    for j in ti.ndrange(1,ny-1):
        for k in ti.static(range(9)):

            #left
            ibc = 0
            jbc = j
            inb = 1
            jnb = j

            rho[ibc,jbc] = rho[inb,jnb]
            f_old[ibc,jbc][k] = f_eq(e[k,0],e[k,1], w[k], rho[ibc,jbc], vel[ibc,jbc][0], vel[ibc,jbc][1]) - \
                            f_eq(e[k,0],e[k,1], w[k], rho[inb,jnb], vel[inb,jnb][0], vel[inb,jnb][1]) + \
                            f_old[inb,jnb][k]

            #right
            ibc = nx-1
            jbc = j
            inb = nx-2
            jnb = j

            rho[ibc,jbc] = rho[inb,jnb]
            f_old[ibc,jbc][k] = f_eq(e[k,0],e[k,1], w[k], rho[ibc,jbc], vel[ibc,jbc][0], vel[ibc,jbc][1]) - \
                            f_eq(e[k,0],e[k,1], w[k], rho[inb,jnb], vel[inb,jnb][0], vel[inb,jnb][1]) + \
                            f_old[inb,jnb][k]

    # top and bottom
    for i in ti.ndrange(nx):
        for k in ti.static(range(9)):

            # bottom 
            ibc = i
            jbc = 0
            inb = i
            jnb = 1

            rho[ibc,jbc] = rho[inb,jnb]
            f_old[ibc,jbc][k] = f_eq(e[k,0],e[k,1], w[k], rho[ibc,jbc], vel[ibc,jbc][0], vel[ibc,jbc][1]) - \
                            f_eq(e[k,0],e[k,1], w[k], rho[inb,jnb], vel[inb,jnb][0], vel[inb,jnb][1]) + \
                            f_old[inb,jnb][k]

            # top
            ibc = i
            jbc = ny-1
            inb = i
            jnb = ny-2

            rho[ibc,jbc] = rho[inb,jnb]
            vel[ibc,jbc][0] = u_c
            f_old[ibc,jbc][k] = f_eq(e[k,0],e[k,1], w[k], rho[ibc,jbc], vel[ibc,jbc][0], vel[ibc,jbc][1]) - \
                            f_eq(e[k,0],e[k,1], w[k], rho[inb,jnb], vel[inb,jnb][0], vel[inb,jnb][1]) + \
                            f_old[inb,jnb][k]





if __name__ == '__main__':

    gui = ti.GUI('LBM-Fluid', (nx, ny))

    init()

    for i in range(20000):
        collide_and_stream()
        update_macro_var()
        get_show_var()
        apply_bc()

        if(i%1000==0):
            print('Step: {:}'.format(i))
            img = show_var.to_numpy()
            gui.set_image(cm.rainbow((img+0.5)/1.5))
            gui.show()

    img = show_var.to_numpy()

    y_ref, u_ref = np.loadtxt('data/ghia1982.dat',unpack=True,skiprows=2,usecols=(0,2))

    fig,axes = plt.subplots(nrows=1,ncols=1,figsize=(4,3),dpi=200)

    axes.plot(np.linspace(0,1.0,ny),img[nx//2,:],'b-',label='LBM')
    axes.plot(y_ref,u_ref,'rs',label='Ghia et al. 1982')
    axes.legend()
    axes.set_xlabel(r'Y')
    axes.set_ylabel(r'U')

    plt.tight_layout()

    plt.show()


