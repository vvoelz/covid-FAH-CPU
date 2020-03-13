import os, sys
import numpy as np

#from math import *
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis



usage = """Usage:    python mbar_pmf_slower.py datadir outfile

    INPUT
    datadir  - the name of directoy containing the pull_slower protocol LapG binder umbrella simulations
               NOTE: must contain pullx*.xvg and energy*.xvg files

    OUTPUT
    outfile  -  name of the output file to save pmf results

    Try:
    $ python mbar_pmf_slower.py pull_series_slower_35 pull_series_slower_35.pmf.dat

    """


if len(sys.argv) < 3:
    print usage
    sys.exit(1)


datadir = sys.argv[1]
outfile = sys.argv[2]


######## Functions #############


def get_distances(xvgfile, skip=None, stride=1,  gmx_v512=True):
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
    #print 'pullx', pullx
    """
    # time (ps) dX (nm)         dY (nm)         dZ (nm)
    0.0000      0.236968        0.955994        -0.792
    1.0000      0.274845        0.997569        -0.834498
    2.0000      0.293425        1.07143 -0.82839
    """

    print 'pullx', pullx
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

def get_dhdl(xvgfile, skip=None):
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
        energies = dhdl[Ind,4:-1]  # ignore pV term in last column
    else:
        times = dhdl[:,0]
        thermo_indices = dhdl[:,1].astype(int)
        energies = dhdl[:,4:-1]

    return times, thermo_indices, energies 



def distance2bin(x, pmf_distances):
    """Given an array of pmf distance edges, return the bin index that the distance is closest to."""

    return  np.argsort(np.abs(pmf_distances - x))[0]
    

###########################
# Main

# Physical constants (in kJ)
kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K


rvalues = [(1.25 + i*0.1) for i in range(0,30)]  # for slower pulls
#rvalues = [(1.25 + i*0.1) for i in range(0,28)]  # for slower pulls
print 'rvalues', rvalues



#### For testing, we only have the first  ###  restraint sims
#rvalues = rvalues[0:60]
#rvalues = rvalues[0:32]

unique_rvalues = np.unique( rvalues )
nthermo_indices = len(unique_rvalues)

r_indices = np.arange(len(rvalues))+2 # indexing starts at 2
print 'r_indices', r_indices

kvalues = [200.0 for r in rvalues]  # kJ/mol/nm^2



# MBAR parameters

# Parameters
temperature = 300.
K = len(rvalues) # number of umbrellas VAV: it's okay to have some of the umbrellas the same
# maximum number of snapshots/simulation:
N_max = 2001 # for pull_series_slower
# N_max = 1001 # for pull_series
T_k = np.ones(K,float)*temperature # inital temperatures are all equal 
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
print 'beta',  beta

# GOAL: build a pmf over these bin centers
pmf_distances = np.arange(1.20, 3.90, 0.1)
print 'pmf_distances', pmf_distances
nbins = len(pmf_distances)


# Allocate storage for simulation data
N_k = np.zeros([K], np.int32) # N_k[k] is the number of snapshots from umbrella simulation k
K_k = np.zeros([K], np.float64) # K_k[k] is the spring constant (in kJ/mol/nm**2) for umbrella simulation k
beta_k = np.zeros([K], np.float64) # beta_k[k] is the beta = 1/kT of simulation k 
# rvalues[k] is the spring center location (in nm) for umbrella simulation k
x_kn = np.zeros([K,N_max], np.float64) # x_kn[k,n] is the Val122_CA-TRP_CA distance (in nm) for snapshot n from umbrella simulation k
u_kn = np.zeros([K,N_max], np.float64) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
g_k = np.zeros([K],np.float32);


# read in the umbrella sampling data
skip_time = 0.0
stride = 1   
# VAV: NOTE the energies are sampled every 2 ps for pull_series; use stride=2
# VAV: NOTE the energies are sampled every 1 ps for pull_series_slower; use stride=1


if (1):
    for k in range(len(rvalues)):

          r  = rvalues[k]
          r_index = r_indices[k]
          K_k[k] = kvalues[k]
          beta_k[k] = beta

          # read in the pullx distance data
          pullfile = os.path.join(datadir,'pullx%d.xvg'%r_index)
          pull_times, distances = get_distances(pullfile, skip=skip_time, stride=stride, gmx_v512=False)

          print "pull_times.shape, distances.shape", pull_times.shape, distances.shape
          nmax = distances.shape[0]
          for n in range(nmax):
              x_kn[k,n] = distances[n]

          # Read energies
          filename = os.path.join(datadir,'energy%d.xvg'%r_index)

          print "Reading %s..." % filename
          infile = open(filename, 'r')
          lines = infile.readlines()
          infile.close()
          # Parse data.
          n = 0
          for line in lines:
              if line[0] != '#' and line[0] != '@':
                  tokens = line.split()
                  u_kn[k,n] = beta_k[k] * (float(tokens[2]) - float(tokens[1])) # reduced potential energy without umbrella restraint
                  n += 1

          # print 'nmax', nmax, 'n', n
          N_k[k] = n


          # Compute correlation times for potential energy 
          g_k[k] = timeseries.statisticalInefficiency(u_kn[k,0:N_k[k]])
          print "Correlation time for set %5d is %10.3f" % (k,g_k[k])
          indices = timeseries.subsampleCorrelatedData(u_kn[k,0:N_k[k]], g=g_k[k])

          # Subsample data.
          N_k[k] = len(indices)
          u_kn[k,0:N_k[k]] = u_kn[k,indices]
          x_kn[k,0:N_k[k]] = x_kn[k,indices]
          print 'Subsampling so that N_k[k=%d] ='%k,N_k[k]


    N_max = np.max(N_k) # shorten the array size
    u_kln = np.zeros([K,K,N_max], np.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l
  
    # Set zero of u_kn -- this is arbitrary.
    u_kn -= u_kn.min()

    print "Binning data..."
    # Bin data
    bin_kn = np.zeros([K,N_max], np.int32)
    for k in range(K):
        for n in range(N_k[k]):
            # Compute bin assignment.
            bin_kn[k,n] =  distance2bin(x_kn[k,n], pmf_distances)

    # Evaluate reduced energies in all umbrellas
    print "Evaluating reduced potential energies..."
    for k in range(K):
        for n in range(N_k[k]):
          for l in range(K):
            # Compute energy of snapshot n from simulation k in umbrella potential l
            u_kln[k,l,n] = u_kn[k,n] + beta_k[k] * (K_k[l]/2.0) * (x_kn[k,n] - rvalues[l])**2

        print 'k', k, 'u_kln', u_kln[k,:,:N_k[k]]

    # Initialize MBAR.
    print "Running MBAR..."
    mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive')

    # Are there any empty bins? 
    bincounts, bin_edges = np.histogram(bin_kn, bins=range(len(pmf_distances)))
    print len(pmf_distances), pmf_distances
    print len(bincounts), bincounts
     

    # Compute PMF in unbiased potential (in units of kT).
    (f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins)

    # reset zero of PMF to be the min in the first few bins 
    f_i -= min(f_i[0:4])

    # print out PMF
    print "PMF (in units of kT)"
    print "%8s %8s %8s" % ('distance (nm)', 'f (kT)', 'df (kT)')
    for i in range(nbins):
        print "%8.3f %8.3f %8.3f" % (rvalues[i], f_i[i], df_i[i])

    # write PMF to file
    fout = open(outfile, 'w')
    fout.write('# distance (nm)\tf (kT)\tdf (kT)\n')
    for i in range(nbins):
        fout.write("%8.4f\t%8.4f\t%8.4f\n" % (rvalues[i], f_i[i], df_i[i]) )
    fout.close()
    print 'Wrote', outfile, '.'
  
    # plt.errorbar(pmf_distances, f_i, yerr=df_i, fmt='--o', label=datadir)


"""
from matplotlib import pyplot as plt
plt.figure()

# read in  pmf_distances, f_i, df_i
#plt.errorbar(pmf_distances, f_i, yerr=df_i, fmt='--o', label=datadir)

plt.xlabel('distance (nm)')
#plt.xlabel('bin index')
plt.ylabel('free energy (units kT)')
#plt.plot(pmf_distances,tram_obj.stationary_distribution, '-')
plt.legend(loc='best')
#plt.show()
plt.savefig('junk.pdf')
"""



