import os, sys, glob
import numpy as np

import xvg_tools

sys.path.append('../scripts')
from expanded_v3 import *

import argparse, textwrap

import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis


def calc_dG_rest(run_clone_datadir, mdp_file, outdir=None, maxlastgens=5, temperature=298.15, verbose=False):
    """
    Returns the free energy of harmonic restraint in kT.
    
    INPUT
    run_clone_datadir  - the pathname of the PROJ*/RUN*/CLONE* data in ~/server2/data.
                         Example: /home/server/server2/data/SVR2616698070/PROJ14727/RUN0/CLONE0
    mdp_file           - the pathname of the *.mdp file that contains the force constant equilibrium restraint distance

    PARAMETERS  
    maxlastgens        - the max last results* gens to use in the calculation.  For instance, 
                         if the trajectory goes out to results35, and maxlastgens = 3,  then results33-35 will be used.               
    outdir             - If specified, the pathname of a directory to write the energies, states, and distances.
                         If outdir=None (default), then no output files will be written.
    temperature        - the temperature in Kelvin

    RETURNS
    dG_rest            - the free energy of harmonic restraint in kT.
    dG_rest_sigma     - the uncertainty (std dev) of dG_rest in kT

    
    REQUIREMENTS

    pymbar    - see installation instructions at  https://github.com/choderalab/pymbar
          NOTE: as of 2020-505-14, `$ conda install -c omnia pymbar` fails, but `$ pip install pymbar` works

    """

    if outdir != None:
        # Make an output directory if it doesn't exist
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    # Find all the results* directories in the PROJ*/RUN*/CLONE* 
    
    #runclonefor each clone, grab the data in the dhdl.xvg and pullx.xvg in each gen
    resultdirs = []
    resultdirs1 = glob.glob( os.path.join(run_clone_datadir, 'results?') )
    resultdirs1.sort()
    resultdirs += resultdirs1

    resultdirs2 = glob.glob( os.path.join(run_clone_datadir, 'results??') )
    resultdirs2.sort()
    resultdirs += resultdirs2

    resultdirs3 = glob.glob( os.path.join(run_clone_datadir, 'results???') )
    resultdirs3.sort()
    resultdirs += resultdirs3

    # test to see if the last results directory has any *.xvg files in it...
    # if not, remove it -- it's empty!
    if len(resultdirs) > 0:
        if len(glob.glob(os.path.join(resultdirs[-1], '*.xvg'))) == 0:
            resultdirs.pop()

    # We will only analyze the last maxlastgens results
    while len(resultdirs) > maxlastgens:
        resultdirs.pop(0)

    if verbose:
        for resultdir in resultdirs:
            print('\t',resultdir)

    ###### ENERGIES ########

    result_dhdlfiles = [os.path.join(resultdir,'dhdl.xvg') for resultdir in resultdirs if os.path.exists(os.path.join(resultdir,'dhdl.xvg'))]
    if len(result_dhdlfiles) > 0:

        ### Scrape and collect all the dhdl energies!
        for resultdir in resultdirs[0:1]:
            dhdl_xvgfile =  os.path.join(resultdir, 'dhdl.xvg')
            time_in_ps, states, energies = xvg_tools.get_dhdl(dhdl_xvgfile)
            if verbose:
                print(resultdir)
                print('\ttime_in_ps', time_in_ps)
                print('\tstates', states)
                print('\tenergies.shape', energies.shape)
                print('\tenergies', energies)

        for resultdir in resultdirs[1:]:
            dhdl_xvgfile =  os.path.join(resultdir, 'dhdl.xvg')
            more_time_in_ps, more_states, more_energies = xvg_tools.get_dhdl(dhdl_xvgfile)
            time_in_ps = np.concatenate( (time_in_ps, more_time_in_ps[1:]), axis=0)
            states     = np.concatenate( (states, more_states[1:]), axis=0 )
            energies   = np.concatenate( (energies, more_energies[1:,:]), axis=0 )

        if outdir != None:
            states_outfile = os.path.join(outdir, 'states.npy')
            np.save(states_outfile, states)
            print('Wrote:', states_outfile)

            energies_outfile = os.path.join(outdir, 'energies.npy')
            np.save(energies_outfile, energies)
            print('Wrote:', energies_outfile)

        ### Scrape and collect all the pullx distances!
        for resultdir in resultdirs[0:1]:
            pullx_xvgfile =  os.path.join(resultdir, 'pullx.xvg')
            times, distances = xvg_tools.get_distances(pullx_xvgfile)
            if verbose:
                print('distances', distances)

        for resultdir in resultdirs[1:]:
            pullx_xvgfile =  os.path.join(resultdir, 'pullx.xvg')
            more_times, more_distances = xvg_tools.get_distances(pullx_xvgfile)
            times = np.concatenate( (times, more_times[1:]), axis=0)
            distances = np.concatenate( (distances, more_distances[1:]), axis=0)

        if outdir != None:
            distances_outfile = os.path.join(outdir, 'distances.npy')
            np.save(distances_outfile, distances)
            print('Wrote:', distances_outfile)

    if verbose:
        print('Are all 40 intermediates sampled adequately?')
        print('np.histogram(states, bins=np.arange(40))', np.histogram(states, bins=np.arange(40)))

    ##########################################
    ### Do the MBAR calculation here!
    ###########################################

    # Get the temperature and equilibrium restraint distance r0 from the mdp_file
    e = expanded_ensemble_mdpfile()
    e.read_parms_from_mdpfile(mdp_file, VERBOSE=verbose)

    kvalue = float(e.pull_coord1_k)
    r0 = float(e.pull_coord1_init)
    if verbose:
        print('kvalue', kvalue, 'r0', r0) 

    dG_rest, sigma_dG_rest = estimate_free_energy(states, energies, distances, kvalue=kvalue, r0=r0, temperature=temperature, verbose=verbose)

    return dG_rest, sigma_dG_rest


def estimate_free_energy(states, energies, distances,
                         kvalue=800.0, r0=0.4, temperature=300., verbose=False):
    """Use MBAR to estimate the free energy vs. lambda.
    N is the number of samples
    K is the number of thermodynamic states

    INPUTS
    states    - state indices in a numpy array of shape (N,)
    energies  - a numpy array of shape (N, K) with the dhdl info
    distances - restraint distances numpy array of shape (N,)

    PARAMETERS
    kvalue      - the harmonic force constant in kJ/nm^2 (Default: 400.0)
    r0          - the restraint distance umbrella center, in nm.
    temperature - in K (Default: 300.0)
    verbose     - print verbose output

    OUTPUT
    """

    ###########################
    # Main

    # Physical constants (in kJ)
    kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K

    # In addition to the given K thermo ensembles with indices 0,..K-1,
    # there is one more -- the *unbiased* ensemble with no harmonic restraint.
    #  -- let's make it state index K
    K = energies.shape[1] + 1
    unbiased_state_index = K - 1


    # maximum number of snapshots/simulation:
    N_max = energies.shape[0]
    ### TO DO in the future collect *all* run and clone energies and flatten

    T_k = np.ones(K,float)*temperature # inital temperatures are all equal
    beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
    if verbose:
        print('beta',  beta)

    # Allocate storage for simulation data
    N_k = np.zeros([K], np.int32) # N_k[k] is the number of snapshots from umbrella simulation k
    # rvalues[k] is the spring center location (in nm) for umbrella simulation k
    x_kn = np.zeros([K,N_max], np.float64) # x_kn[k,n] is the Val122_CA-TRP_CA distance (in nm) for snapshot n from umbrella simulation k
    u_kn = np.zeros([K,N_max], np.float64) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
    g_k = np.zeros([K],np.float32);


    ### To do MBAR, we need to convert to data to u_kn format
    if verbose:
        print('np.argsort(states)', np.argsort(states))

    Ind = np.argsort(states)
    energies_sorted = energies[Ind,:]
    states_sorted = states[Ind]
    distances_sorted = distances[Ind]
    for k in range(K-1):

        # Count how many snapshots belong to each k
        N_k[k] = np.where(states_sorted == k, 1, 0).sum()

        # fill the energies
        u_kn[k, :] = energies_sorted[:, k]

    # for the last (unbiased) ensemble (with no samples), subtract the harmonic potential from the
    # state index 0 (totally coupled) eneries
    u_kn[K-1, :] = u_kn[0, :] - beta * (kvalue/2.0) * (distances_sorted - r0)**2
    if verbose:
        print('u_kn', u_kn)
        print('N_k', N_k)

    # Initialize MBAR.
    if verbose:
        print('Running MBAR...')
    mbar = pymbar.MBAR(u_kn, N_k, verbose=verbose)  #, maximum_iterations=100000, initialize='BAR')  
    
    # MBAR(u_kn, N_k, maximum_iterations=10000, relative_tolerance=1e-07, verbose=False, initial_f_k=None, solver_protocol=None, initialize='zeros', x_kindices=None, **kwargs))

    # Set zero of u_kn -- this is arbitrary.
    u_kn -= u_kn.min()

    results = mbar.getFreeEnergyDifferences()
    Deltaf_ij, dDeltaf_ij = results[0], results[1]
    if verbose:
        print('Deltaf_ij, dDeltaf_ij', Deltaf_ij, dDeltaf_ij)

    df, sigma_df = np.zeros(K), np.zeros(K)
    for i in range(K-1):
        #  print('Deltaf_%d,%d = '%(i,i+1), Deltaf_ij[i,i+1], '+/-', dDeltaf_ij[i,i+1])
        df[i+1] = df[i] + Deltaf_ij[i,i+1]
        sigma_df[i+1] = dDeltaf_ij[0,i+1]

    dG_rest       = df[0] - df[-1]      # THIS should be the same as +f[-1] of the MBAR object!
    sigma_dG_rest = sigma_df[-1] 

    if verbose:
        print('Delta f (norest, lam=1 -> rest, lam=1) =', df[0] - df[-1])
        print('dDelta f (norest, lam=1 -> rest, lam=1) =', sigma_df[-1])

    return dG_rest, sigma_dG_rest


#############################################

# Main

if __name__ == '__main__':

    import argparse, textwrap

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\

    Computes an MBAR estimate of the free energy of *adding* the harmonic restraint, $\Delta G_rest

    EXAMPLE
    $ python harmonic_rest_correction.py ~/server2/projects/p14727 ~/server2/data/SVR2616698070 14727 0 0 

      ''' ))

    parser.add_argument('setupdir', type=str, help='The setup dir in which RUN0etup')
    parser.add_argument('datadir', type=str, help='The direcory where projdata is saved')
    parser.add_argument('projnum', type=int, help='The project number')
    parser.add_argument('run', type=int, help='The run number')
    parser.add_argument('clone', type=int, help='The clone number')
    parser.add_argument('-o', '--outdir', required=False)
    parser.add_argument('-n', '--maxlastgens', type=int, required=False, help='Limit data analyzed to the last n gens. (Default: 5)')

    parser.add_argument('--verbose', dest='verbose', action='store_true',
                    help='If specified, print output verbosely')
    args = parser.parse_args()

    # process the input arguments 
    index_file = os.path.join(args.setupdir, 'RUN%d/index.ndx'%args.projnum)
    mdp_file = os.path.join(args.setupdir, 'RUN%d/prod.mdp'%args.run)
    run_clone_datadir = os.path.join(args.datadir, 'PROJ%d/RUN%d/CLONE%d'%(args.projnum, args.run, args.clone))

    if args.verbose:
        print('### INPUTS from argparse ###')
        print('args.setupdir', args.setupdir)
        print('args.datadir', args.datadir)
        print('args.projnum', args.projnum)
        print('args.run', args.run)
        print('args.clone', args.clone)
        print('args.outdir', args.outdir)
        print('args.maxlastgens', args.maxlastgens)
        print('args.verbose', args.verbose)
        print('--------')
        print('\tindex_file =', index_file)
        print('\tmdp_file =', mdp_file)
        print('\trun_clone_datadir =', run_clone_datadir)


    if len(sys.argv) < 4:
        print(usage)
        sys.exit(1)

    # Calculate the free energy of restraint
    dG_rest, sigma_dG_rest = calc_dG_rest(run_clone_datadir, mdp_file, outdir=args.outdir, verbose=args.verbose, maxlastgens=args.maxlastgens) 

    print('# PROJ', args.projnum)
    print('# RUN', args.run)
    print('# CLONE', args.clone)
    print('dG_rest = %6.6f +/- %6.6f kT.'%(dG_rest, sigma_dG_rest))


