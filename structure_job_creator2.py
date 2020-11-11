#!/usr/bin/env python

###Imports###
import numpy as np
import subprocess
import os

###Definitions###
reps = 5 # per k
min_k = 1
max_k = 1
burnin = 100000
mcmc_reps = 500000
num_inds = 27
num_loci = 30087
str_file = '/fs/project/PAA0202/Duckett/Chapter2/RAD/structure/microtus_27/microtus_29_m50_sm75_mmd107_rm2.recode.1snp.str'
walltime = '30:00:00'
nodes = 1
ppn = 28
account = 'PAA0202'

k_values = np.arange(min_k, max_k + 1).tolist()


###Functions###

def make_mainparams(pops, rep, burnin, mcmc_reps, num_inds, num_loci):
    out_name = 'mainparams_k%s_rep%s' % (str(pops), str(rep)) # name of param file
    out_file = 'out_k%s_rep%s' % (str(pops), str(rep)) # name of structure output
    with open(out_name, 'at') as out:
        out.writelines('''\
KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
FILE extraparams.


"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean 
        (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!) 


Basic Program Parameters

#define MAXPOPS    {0}      // (int) number of populations assumed
#define BURNIN    {1}   // (int) length of burnin period
#define NUMREPS   {2}   // (int) number of MCMC reps after burnin

Input/Output files

#define INFILE   {3}   // (str) name of input data file
#define OUTFILE  {4}  //(str) name of output data file

Data file format

#define NUMINDS    {5}    // (int) number of diploid individuals in data file
#define NUMLOCI    {6}    // (int) number of loci in data file
#define PLOIDY       2    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line


#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   1     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says
                            whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier

#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data 
                            before the genotype data start.

#define MARKERNAMES      0  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances 
                            // between loci


Advanced data file options

#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data



Command line options:

-m mainparams
-e extraparams
-s stratparams
-K MAXPOPS 
-L NUMLOCI
-N NUMINDS
-i input file
-o output file
-D SEED



'''.format(str(pops), str(burnin), str(mcmc_reps), str_file, out_file, str(num_inds), str(num_loci)))



def make_extraparams(pops, rep, seed):
    out_name = 'extraparams_k%s_rep%s' % (str(pops), str(rep)) # name of param file
    with open(out_name, 'at') as out:
        out.writelines('''\
EXTRA PARAMS FOR THE PROGRAM structure.  THESE PARAMETERS CONTROL HOW THE
PROGRAM RUNS.  ATTRIBUTES OF THE DATAFILE AS WELL AS K AND RUNLENGTH ARE 
SPECIFIED IN mainparams.

"(int)" means that this takes an integer value.
"(d)"   means that this is a double (ie, a Real number such as 3.14).
"(B)"   means that this variable is Boolean 
        (ie insert 1 for True, and 0 for False).

PROGRAM OPTIONS

#define NOADMIX     0 // (B) Use no admixture model (0=admixture model, 1=no-admix)
#define LINKAGE     0 // (B) Use the linkage model model 
#define USEPOPINFO  0 // (B) Use prior population information to pre-assign individuals
                            to clusters
#define LOCPRIOR    1 //(B)  Use location information to improve weak data

#define FREQSCORR   1 // (B) allele frequencies are correlated among pops
#define ONEFST      0 // (B) assume same value of Fst for all subpopulations.

#define INFERALPHA  1 // (B) Infer ALPHA (the admixture parameter)
#define POPALPHAS   1 // (B) Individual alpha for each population
#define ALPHA     1.0 // (d) Dirichlet parameter for degree of admixture 
                            (this is the initial value if INFERALPHA==1).

#define INFERLAMBDA 0 // (B) Infer LAMBDA (the allele frequencies parameter)
#define POPSPECIFICLAMBDA 0 //(B) infer a separate lambda for each pop 
                    (only if INFERLAMBDA=1).
#define LAMBDA    1.0 // (d) Dirichlet parameter for allele frequencies 




PRIORS

#define FPRIORMEAN 0.01 // (d) Prior mean and SD of Fst for pops.
#define FPRIORSD   0.05  // (d) The prior is a Gamma distribution with these parameters

#define UNIFPRIORALPHA 1 // (B) use a uniform prior for alpha;
                                otherwise gamma prior
#define ALPHAMAX     10.0 // (d) max value of alpha if uniform prior
#define ALPHAPRIORA   1.0 // (only if UNIFPRIORALPHA==0): alpha has a gamma 
                            prior with mean A*B, and 
#define ALPHAPRIORB   2.0 // variance A*B^2.  


#define LOG10RMIN     -4.0   //(d) Log10 of minimum allowed value of r under linkage model
#define LOG10RMAX      1.0   //(d) Log10 of maximum allowed value of r
#define LOG10RPROPSD   0.1   //(d) standard deviation of log r in update
#define LOG10RSTART   -2.0   //(d) initial value of log10 r


USING PRIOR POPULATION INFO \\(USEPOPINFO\\)

#define GENSBACK    2  //(int) For use when inferring whether an indiv-
                        idual is an immigrant, or has an immigrant an-
                        cestor in the past GENSBACK generations.  eg, if 
                        GENSBACK==2, it tests for immigrant ancestry 
                        back to grandparents. 
#define MIGRPRIOR 0.01 //(d) prior prob that an individual is a migrant 
                            (used only when USEPOPINFO==1).  This should 
                            be small, eg 0.01 or 0.1.
#define PFROMPOPFLAGONLY 0 // (B) only use individuals with POPFLAG=1 to update	P.
                                This is to enable use of a reference set of 
                                individuals for clustering additional "test" 
                                individuals.

LOCPRIOR MODEL FOR USING LOCATION INFORMATION

#define LOCISPOP      1    //(B) use POPDATA for location information 
#define LOCPRIORINIT  1.0  //(d) initial value for r, the location prior
#define MAXLOCPRIOR  20.0  //(d) max allowed value for r




OUTPUT OPTIONS

#define PRINTNET     1 // (B) Print the "net nucleotide distance" to screen during the run
#define PRINTLAMBDA  1 // (B) Print current value(s) of lambda to screen
#define PRINTQSUM    1 // (B) Print summary of current population membership to screen

#define SITEBYSITE   0  // (B) whether or not to print site by site results. 
                    (Linkage model only) This is a large file!
#define PRINTQHAT    1  // (B) Q-hat printed to a separate file.  Turn this 
                        on before using STRAT.
#define UPDATEFREQ   100  // (int) frequency of printing update on the screen.
                                Set automatically if this is 0.
#define PRINTLIKES   0  // (B) print current likelihood to screen every rep
#define INTERMEDSAVE 0  // (int) number of saves to file during run

#define ECHODATA     1  // (B) Print some of data file to screen to check
                            that the data entry is correct.
(NEXT 3 ARE FOR COLLECTING DISTRIBUTION OF Q:)
#define ANCESTDIST   1  // (B) collect data about the distribution of an-
                            cestry coefficients (Q) for each individual
#define NUMBOXES   1000 // (int) the distribution of Q values is stored as 
                            a histogram with this number of boxes. 
#define ANCESTPINT 0.90 // (d) the size of the displayed probability  
                            interval on Q (values between 0.0--1.0)



MISCELLANEOUS

#define COMPUTEPROB 1     // (B) Estimate the probability of the Data under 
                            the model.  This is used when choosing the 
                            best number of subpopulations.
#define ADMBURNIN  500    // (int) [only relevant for linkage model]: 
                            Initial period of burnin with admixture model (see Readme)
#define ALPHAPROPSD 0.025 // (d) SD of proposal for updating alpha
#define STARTATPOPINFO 0  // Use given populations as the initial condition 
                            for population origins.  (Need POPDATA==1).  It 
                            is assumed that the PopData in the input file 
                            are between 1 and k where k<=MAXPOPS.
#define RANDOMIZE      0  // (B) use new random seed for each run 
#define SEED        {0}  // (int) seed value for random number generator 
                        (must set RANDOMIZE=0) 
#define METROFREQ    10   // (int) Frequency of using Metropolis step to update
                            Q under admixture model (ie use the metr. move every
                            i steps).  If this is set to 0, it is never used.
                            (Proposal for each q^(i) sampled from prior.  The 
                            goal is to improve mixing for small alpha.)
#define REPORTHITRATE 0 //   (B) report hit rate if using METROFREQ
'''.format(str(seed)))


def make_pbs(walltime, nodes, ppn, pops, rep, account):
    out_name = 'structure_k%s_rep%s.pbs' % (str(pops), str(rep)) # name of pbs file
    with open(out_name, 'at') as out:
        out.writelines('''\
#!/usr/bin/bash
#
#PBS -l walltime={0}
#PBS -l nodes={1}:ppn={2}
#PBS -N k{3}_r{4}
#PBS -j oe
#PBS -A {5}

cd $PBS_O_WORKDIR

/fs/project/PAA0202/Programs/STRUCTURE/structure -m mainparams_k{3}_rep{4} -e extraparams_k{3}_rep{4} -o out_k{3}_r{4} -k {3}
'''.format(walltime, str(nodes), str(ppn), str(pops), str(rep), account))


###Main###
for k in range(0, len(k_values)):
    k_val = k_values[k] # get k value
    cmd_mkdir = 'mkdir k%s' % k_val
    #print(cmd_mkdir)
    subprocess.call(cmd_mkdir, shell = True) # create directory
    os.chdir('k%s' % k_val)
    for i in range(0,reps):
        seed = int(np.random.randint(0, 1000000, 1)) # generate random seed
        make_mainparams(k_val, i, burnin, mcmc_reps, num_inds, num_loci) # make mainparams file
        make_extraparams(k_val, i, seed) # make extraparams file
        make_pbs(walltime, nodes, ppn, k_val, i, account) # make pbs file
        cmd_submit = 'qsub structure_k%s_rep%s.pbs' % (k_val, i)
        print(cmd_submit)
        #subprocess.call(cmd_submit, shell = True)
    os.chdir('..')


