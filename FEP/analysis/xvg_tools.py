import os, sys
import numpy as np


######## Functions #############


def get_distances(xvgfile, skip=None, stride=1,  gmx_v512=False):
    """Extract time points (ps) and distances (nm) from a pullx.xvg file

    OPTIONS
    skip	time in ps.  Ignore data before this time

    RETURNS
    times	time points in ps
    distances	distances in nm
    """

    if not os.path.exists('tmp'):
        os.mkdir('tmp')

    os.system("cat %s | grep -v '@' | grep -v '#' > tmp/pullx.dat"%xvgfile)
    pullx = np.loadtxt('tmp/pullx.dat')
    #print('pullx', pullx)
    """
    # time (ps) dX (nm)         dY (nm)         dZ (nm)
    0.0000      0.236968        0.955994        -0.792
    1.0000      0.274845        0.997569        -0.834498
    2.0000      0.293425        1.07143 -0.82839
    """

    print('pullx', pullx)
    # compute the distances
    if skip != None:
        Ind = (pullx[:,0] >= skip)  # skip the first few snapshots
        times = pullx[Ind,0]
        if gmx_v512:
            distances = pullx[Ind,1]
        else:
            distances = (pullx[Ind,1]**2 + pullx[Ind,2]**2 + pullx[Ind,3]**2)**0.5
    else:
        times = pullx[:,0]
        if gmx_v512:
            distances = pullx[:,1]
        else:
            distances = (pullx[:,1]**2 + pullx[:,2]**2 + pullx[:,3]**2)**0.5

    return times[::stride], distances[::stride]

def get_dhdl(xvgfile, skip=None, ignore_last_column=False):
    """Extract time points and all the \Delta energies from a dhdl.xvg file.

    OPTIONS
    skip        time in ps.  Ignore data before this time

    RETURNS
    times       	time points in ps
    thermo_indices	index of the thermodynamic ensemble being sampled from
    energies    	in kJ/mol. columns in thermo index (lambda)
    """

    if not os.path.exists('tmp'):
        os.mkdir('tmp')
    os.system("cat %s | grep -v '@' | grep -v '#' > tmp/dhdl.dat"%xvgfile)
    dhdl = np.loadtxt('tmp/dhdl.dat')

    if skip != None:
        Ind = (dhdl[:,0] >= skip)  # skip the first few snapshots
        times = dhdl[Ind,0]
        thermo_indices = dhdl[Ind,1].astype(int)
        if ignore_last_column:
            energies = dhdl[Ind,3:-1]  # ignore pV term in last column
        else:
            energies = dhdl[Ind,3:] 

    else:
        times = dhdl[:,0]
        thermo_indices = dhdl[:,1].astype(int)
        if ignore_last_column:
            energies = dhdl[:,3:-1]
        else:
            energies = dhdl[:,3:]

    return times, thermo_indices, energies 


