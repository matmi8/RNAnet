import os, sys

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as mpl
import random as rm



import Network_functions as NF
import Building_functions as BF

DEFAULT = 2
N_ARG = 5

def  requested_argv():
    print("\nError! Wrong number of argument.\n")
    print("python2 Random_Graph.py template_file  label N_rip savedir")
    print("- template_file, file where original network is stored;")
    print("- Label, name for the generated networks: Label_links_ER_n.dat")
    print("- N_rip, number of random network to be generated;")
    print("- savedir, directory where to save networks.")
    print("\n")


def main():

    if len(sys.argv) != N_ARG:
        requested_argv()
        sys.exit(1)
    
    net_file = sys.argv[1]

    Lab = sys.argv[2]
    N_rip = int(sys.argv[3]) 
   
    savedir = sys.argv[4]
   
    os.system("mkdir %s"%savedir)
    
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
   
    
    print("ER Randomization")
    for i in range(N_rip):
        sys.stdout.write("ER -> %i\n"%i)
        new_link = BF.ER_rand(links)
        name = "%s/%s_links_ER_%04d.dat"%(savedir,Lab,i)
        np.savetxt(name, new_link, fmt='%d')

    sys.stdout.write("Done!\n")

main()
