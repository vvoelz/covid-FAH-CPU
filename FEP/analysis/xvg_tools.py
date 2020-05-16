import os, sys
import numpy as np

import codecs


######## Functions #############


def get_distances(xvgfile, skip=None, stride=1,  gmx_v512=False, verbose=False):
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
    """
    # time (ps) dX (nm)         dY (nm)         dZ (nm)
    0.0000      0.236968        0.955994        -0.792
    1.0000      0.274845        0.997569        -0.834498
    2.0000      0.293425        1.07143 -0.82839
    """

    if verbose:
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

def get_dhdl(xvgfile, skip=None, ignore_last_column=False, verbose=False):
    """Extract time points and all the \Delta energies from a dhdl.xvg file.

    OPTIONS
    skip        time in ps.  Ignore data before this time

    RETURNS
    times       	time points in ps
    thermo_indices	index of the thermodynamic ensemble being sampled from
    energies    	in kJ/mol. columns in thermo index (lambda)
    """


    # Parse the header to figure out which columns need to be selected.

    ## HERE is an example header:
    """
    $ head /home/server/server2/data/SVR2616698070/PROJ14399/RUN0/CLONE0/results1/dhdl.xvg
# This file was created Mon Apr 13 20:19:35 2020
# Created by:
# GROMACS:      GROMACS, VERSION 5.0.4-20191026-456f0d636-unknown
# GROMACS is part of G R O M A C S:
#
# Groningen Machine for Chemical Simulation
#
@    title "dH/dxl\\f{} and xD\\f{}H"
@    xaxis  label "Time (ps)"
@    yaxis  label "dH/dxl\f{} and xD\f{}H (kJ/mol [xl\f{}]\S-1N)"
base) server@vav16:~/git/covid-FAH-CPU/FEP/analysis$ head -30 /home/server/server2/data/SVR2616698070/PROJ14399/RUN0/CLONE0/results1/dhdl.xvg
# This file was created Mon Apr 13 20:19:35 2020
# Created by:
# GROMACS:      GROMACS, VERSION 5.0.4-20191026-456f0d636-unknown
# GROMACS is part of G R O M A C S:
#
# Groningen Machine for Chemical Simulation
#
@    title "dH/dxl\f{} and xD\f{}H"
@    xaxis  label "Time (ps)"
@    yaxis  label "dH/dxl\f{} and xD\f{}H (kJ/mol [xl\f{}]\S-1N)"
@TYPE xy
@ subtitle "T = 298.15 (K) "
@ view 0.15, 0.15, 0.75, 0.85
@ legend on
@ legend box on
@ legend loctype view
@ legend 0.78, 0.8
@ legend length 2
@ s0 legend "Thermodynamic state"
@ s1 legend "Total Energy (kJ/mol)"
@ s2 legend "dH/dxl\f{} fep-lambda = 0.0000"
@ s3 legend "dH/dxl\f{} coul-lambda = 1.0000"
@ s4 legend "dH/dxl\f{} vdw-lambda = 0.3000"
@ s5 legend "dH/dxl\f{} restraint-lambda = 1.0000"
@ s6 legend "xD\f{}H xl\f{} to (0.0000, 0.0000, 0.0000, 1.0000)"
@ s7 legend "xD\f{}H xl\f{} to (0.0000, 0.0500, 0.0000, 1.0000)"
@ s8 legend "xD\f{}H xl\f{} to (0.0000, 0.1000, 0.0000, 1.0000)"
@ s9 legend "xD\f{}H xl\f{} to (0.0000, 0.1500, 0.0000, 1.0000)"
@ s10 legend "xD\f{}H xl\f{} to (0.0000, 0.2000, 0.0000, 1.0000)"
...
...
@ s44 legend "xD\f{}H xl\f{} to (0.0000, 1.0000, 0.9200, 1.0000)"
@ s45 legend "xD\f{}H l\f{} to (0.0000, 1.0000, 1.0000, 1.0000)"
0.0000   23 -580557.38 0.0000000 -81.241417 -34.527348 0.0000000 518.31707 514.25540 510.19295 506.13076 502.06874 498.00662 493.94458 489.88258 485.82046 481.75832 477.69633 473.63427 469.57224 465.51012 461.44803 457.38596 453.32390 449.26181 445.19975 441.13766 437.07559 75.244266 13.898085 0.0000000 2.5590450 7.0771922 13.019033 20.049904 27.930075 32.989875 38.260201 43.715070 49.331923 55.091383 60.976425 66.972302 73.065810 81.322855 89.709427 106.79246
1.0000   23 -580577.00 0.0000000 -112.39775 4.6928134 0.0000000 266.55224 260.93236 255.31252 249.69236 244.07251 238.45287 232.83283 227.21304 221.59312 215.97329 210.35339 204.73353 199.11352 193.49368 187.87383 182.25397 176.63403 171.01419 165.39432 159.77443 154.15454 36.679521 6.3133123 0.0000000 4.9197292 9.9308764 16.144522 23.309848 31.236908 36.296356 41.550738 46.977366 52.556435 58.270683 64.104687 70.045036 76.079544 84.253844 92.554440 109.46005
2.0000   23 -583739.06 0.0000000 -82.951981 -20.952242 0.0000000 711.73853 707.59054 703.44297 699.29541 695.14779 691.00046 686.85277 682.70484 678.55752 674.40971 670.26239 666.11458 661.96702 657.81939 653.67189 649.52432 645.37664 641.22907 637.08145 632.93388 628.78626 61.850649 11.049285 0.0000000 3.3994724 8.1834587 14.396025 21.733522 29.969914 35.269312 40.798614 46.531270 52.444097 58.516851 64.731493 71.072275 77.524991 86.281342 95.188484 113.36784
"""
    # The first column is the "x"-axis (time in ps)
    # The remaining columns are "y"-axes, each labeled in the header with '@ s0', etc.

    fin = open(xvgfile, 'r')
    xvg_lines = fin.readlines()
    fin.close()

    header_lines = [line for line in xvg_lines if ((line[0:3] == '@ s') and (line.count('subtitle')==0))]

    # pop off header lines until we find the first "@ s3 legend "xD\f{}H xl\f{} to (0.0000, 0.0000, 0.0000, 1.0000)"
    while header_lines[0].count("to (") == 0:
        header_lines.pop(0)
    if verbose:
        print('header_lines[0]', header_lines[0])

    # The number after the s (i.e. "@ s5") means that index 6 is where the energy data starts!!! 
    energy_start_index = int(header_lines[0][3]) + 1
    if verbose:
        print('energy_start_index', energy_start_index)

    if not os.path.exists('tmp'):
        os.mkdir('tmp')
    os.system("cat %s | grep -v '@' | grep -v '#' > tmp/dhdl.dat"%xvgfile)
    dhdl = np.loadtxt('tmp/dhdl.dat')

    if skip != None:
        Ind = (dhdl[:,0] >= skip)  # skip the first few snapshots
        times = dhdl[Ind,0]
        thermo_indices = dhdl[Ind,1].astype(int)
        if ignore_last_column:
            energies = dhdl[Ind,energy_start_index:-1]  # ignore pV term in last column
        else:
            energies = dhdl[Ind,energy_start_index:] 

    else:
        times = dhdl[:,0]
        thermo_indices = dhdl[:,1].astype(int)
        if ignore_last_column:
            energies = dhdl[:,energy_start_index:-1]
        else:
            energies = dhdl[:,energy_start_index:]

    return times, thermo_indices, energies 


