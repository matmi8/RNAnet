### libraries...
### python libraries

import os, sys

os.environ["MKL_NUM_THREADS"] = "5" 
os.environ["NUMEXPR_NUM_THREADS"] = "5" 
os.environ["OMP_NUM_THREADS"] = "5" 

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as mpl
import random as rm
import glob

### home-made libraries
print(os.getcwd())
sys.path.append("Functions/")
import Network_functions as NF
import Susceptibility_functions as SF
import Observables_functions as OF
import Structure_functions as STF


PLOT = 0
N_ARG = 7

ER_ON = 0
REAL_ON = 1

def  requested_argv():
    print("\nError! Wrong number of argument.\n")
    print("python Cv_analyser.py net_file  param_file label N_rip savedir")
    print("- net_file, file where original network is stored;")
    print("- param_file, file where simulation parameters are stored;")
    print("- label, an identification for the output;")
    print("- savedir, directory where to save networks.")
    print("- Tr, Cv of mirnas and cernas m0, and mu0.")
    print("- Bd, 0, 1, 2 to regulate binding heterogs.")
    print("\n")


    
    
def main():
    
    if len(sys.argv) != N_ARG:
        requested_argv()
        sys.exit(1)
        
    net_file = sys.argv[1]
    param_file = sys.argv[2]
    Lab = sys.argv[3]   
    savedir = sys.argv[4]

    Tr = float(sys.argv[5])
    Bd = int(sys.argv[6])
    
    os.system("mkdir %s"%savedir)

    ### reading options (parameters) from file
    Opt = NF.Options(param_file)        

    Opt.Bd = Bd
    Opt.Cv_beta = Tr
    Opt.Cv_b = Tr

    N_step = 25 #int(Opt.N_step)
    N_rip = 10 #int(Opt.N_rip)
        
    ####################################
    ### Reading df_links...
    #
    # #,gene_m,mir_n,weigth
    #
    # 1,0,219,3
    # 2,0,1299,3
    # 3,0,2177,3
    # -,-,----,-
    #
    ########################

    
    links =  NF.LinksFunc(net_file, Opt.m0, Opt.Bd)
   
    
    #### checking number of of rows and columns...
    sys.stderr.write("links--> %d, %d\n"%(np.shape(links)))
    if(np.shape(links)[1] !=3):
        sys.stderr.write("Error! Links file has more that three columns!\n")
        exit()


    
    ### computing number of mirna and cernas...
    M, N = NF.SpeciesNumbers(links) 

    ### initializing beta vector..
    beta_vet = np.logspace(np.log(Opt.betaMin), np.log(Opt.betaMax), N_step, base= np.e)
    

    #################
    #################
    ###
    ### Network observables
    ###
    #################
    #################
   

    ### intializing res matrix... 
    Obs_REAL = np.zeros((N_step, 3))
    Obs_REAL[:,0] = beta_vet

    Obs_ER = np.zeros((N_step, 3))
    Obs_ER[:,0] = beta_vet

    Obs_REAL_dev = np.zeros((N_step,2))
    Obs_ER_dev = np.zeros((N_step, 2))
    
    for j in range(N_rip):

	sys.stdout.write("Rip -> %i\n"%j)
	sys.stdout.flush()        
	
	if(j > 0):    
            Obs_ER[:,1:] *= (1.*j)
            Obs_REAL[:,1:] *= (1.*j) 
	    Obs_REAL_dev[:,:] *= (1.*j)
	    Obs_ER_dev[:,:] *= (1.*j)

        links_er =  STF.ER_rand(links)
	### start iterating over betas...
        for i in range(N_step):

            ### initializing beta..
            beta = beta_vet[i]
                
            #links_er = BF.ER_rand(links)

            if(ER_ON):
                Xii_mean, Xii_max = OF.CompObs_Xii(links_er, beta, Opt)
                Obs_ER[i,1:] += Xii_mean, Xii_max
                Obs_ER_dev[i,:] +=  Xii_mean**2, Xii_max**2

            if(REAL_ON):
                Xii_mean, Xii_max = OF.CompObs_Xii(links, beta, Opt)
                Obs_REAL[i,1:] += Xii_mean, Xii_max
                Obs_REAL_dev[i,:] +=  Xii_mean**2, Xii_max**2
        
        if(j > 0):    
            Obs_ER[:,1:] /= float(j+1) 
            Obs_REAL[:,1:] /= float(j+1) 
	    Obs_REAL_dev[:,:] /= float(j+1)
	    Obs_ER_dev[:,:] /= float(j+1)

        if(ER_ON):
            savefile = "ObsX_%s_ER_Tr-%.2e_Bd-%i.dat"%(Lab, Tr, Bd)
	    
            Res = np.column_stack([Obs_ER, np.sqrt(Obs_ER_dev - Obs_ER[:,1:]**2)])    
            np.savetxt(savefile, Res, header = "beta, Xii_mean, Xii_max, dev_Xii_mean, dev_Xii_max")

        if(REAL_ON):
            savefile = "ObsX_%s_REAL_Tr-%.2e_Bd-%i.dat"%(Lab, Tr, Bd)
    
	    Res = np.column_stack([Obs_REAL, np.sqrt(Obs_REAL_dev - Obs_REAL[:,1:]**2)])
            np.savetxt(savefile,Res, header = "beta, Xii_mean, Xii_max, dev_Xii_mean, dev_Xii_max")


    sys.stdout.write("Done!\n")

main()
