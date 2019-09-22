import os, sys

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as mpl
import random as rm

sys.path.append("Functions/")

import Network_functions as NF
import Observables_functions as OF
import Structure_functions as STF

PLOT = 0
N_ARG = 7

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

    N_step = int(Opt.N_step)
    N_rip = int(Opt.N_rip)
    
    
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

    Cv_real = OF.Cv_Func(links, Opt)
    sys.stdout.write("Real Cv computed!\n")

    
    Cv_ER = np.zeros((N_step, 3))
    
    for i in range(N_rip):

        print(N_rip - i)

        links_er = STF.ER_rand(links)
           
        Cv_ER += OF.Cv_Func(links_er, Opt)

    Cv_ER /= float(N_rip)
    
    sys.stdout.write("ER Cv computed!\n")

    Res = np.column_stack((Cv_real, Cv_ER[:,1],Cv_ER[:,2]))
    savefile = "Cv_%s_Tr-%.2e_Bd-%i.dat"%(Lab, Tr, Bd)
    
    np.savetxt(savefile, Res,header = "beta Cv_real SNR_real Cv_ER SNR_rand")

    if(PLOT):
        mpl.figure()
        mpl.xscale("log")
        mpl.plot(Cv_real[:,0], Cv_real[:,1], label = "Real")
        mpl.plot(Cv_ER[:,0], Cv_ER[:,1], label = "ER")
        mpl.xlabel("$\\beta$")
        mpl.ylabel("$C_v$")
        mpl.legend()
        mpl.title("%s"%savefile)
        mpl.show()
        

main()
