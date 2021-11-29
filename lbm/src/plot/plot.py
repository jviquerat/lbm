# Generic imports
import os
import math
import numpy              as np
import matplotlib.pyplot  as plt

### ************************************************
### Output 2D flow amplitude
def plot_norm(lattice, val_min, val_max, output_it, dpi):

    # Compute norm
    v = np.sqrt(lattice.u[0,:,:]**2+lattice.u[1,:,:]**2)

    # Mask obstacles
    v[np.where(lattice.lattice > 0.0)] = -1.0
    vm = np.ma.masked_where((v < 0.0), v)
    vm = np.rot90(vm)

    # Plot
    plt.clf()
    fig, ax = plt.subplots(figsize=plt.figaspect(vm))
    fig.subplots_adjust(0,0,1,1)
    plt.imshow(vm,
               cmap = 'RdBu_r',
               vmin = val_min*lattice.u_lbm,
               vmax = val_max*lattice.u_lbm,
               interpolation = 'spline16')

    filename = lattice.png_dir+'u_norm_'+str(output_it)+'.png'
    plt.axis('off')
    plt.savefig(filename, dpi=dpi)
    plt.close()

### ************************************************
### Output 2D flow contour
def plot_contour(lattice, output_it, dpi):

    # Plot
    plt.clf()
    fig, ax = plt.subplots(figsize=plt.figaspect(vm))
    fig.subplots_adjust(0,0,1,1)
    x  = np.linspace(0, 1, lattice.nx)
    y  = np.linspace(0, 1, lattice.ny)
    ux = lattice.u[0,:,:].copy()
    uy = lattice.u[1,:,:].copy()
    uy = np.rot90(uy)
    ux = np.rot90(ux)
    uy = np.flipud(uy)
    ux = np.flipud(ux)
    vm = np.sqrt(ux**2+uy**2)
    plt.contour(x, y, vm, cmap='RdBu_r',
                vmin=0.0, vmax=1.5*lattice.u_lbm)
    filename = lattice.png_dir+'u_ctr_'+str(output_it)+'.png'
    plt.axis('off')
    plt.savefig(filename, dpi=dpi)
    plt.close()

### ************************************************
### Output 2D streamlines
def plot_streamlines(lattice, dpi):

    # The outputted streamplot is rotated and flipped...
    plt.clf()
    fig, ax = plt.subplots(figsize=plt.figaspect(vm))
    fig.subplots_adjust(0,0,1,1)
    ux = lattice.u[0,:,:].copy()
    uy = lattice.u[1,:,:].copy()
    uy = np.rot90(uy)
    ux = np.rot90(ux)
    uy = np.flipud(uy)
    ux = np.flipud(ux)
    vm = np.sqrt(ux**2+uy**2)
    vm = np.rot90(vm)
    x  = np.linspace(0, 1, lattice.nx)
    y  = np.linspace(0, 1, lattice.ny)
    nn = 20
    u  = np.linspace(0, 1, nn)
    str_pts = []
    for a in range(nn):
        for b in range(nn):
            str_pts.append([u[a],u[b]])
            plt.streamplot(x, y, ux, uy,
                           linewidth    = 0.2,
                           cmap         = 'RdBu_r',
                           arrowstyle   = '-',
                           start_points = str_pts,
                           density      = 20)

    filename = lattice.output_dir+'u_stream.png'
    plt.axis('off')
    plt.savefig(filename, dpi=dpi)
    plt.close()
