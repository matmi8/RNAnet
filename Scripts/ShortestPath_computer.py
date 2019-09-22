import os, sys

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import numpy as np
from numpy import ndarray
import scipy.sparse.csgraph as sp

sys.path.append("Functions/")

import Network_functions as NF
import Structure_functions as STF


N_ARG = 3
PATH_PRINT = 0



def  requested_argv():
    print("\nError! Wrong number of argument.")
    print("\npython2 shortestPaths.py network_file label")
    print("- network_file, the file where to find the network links;")
    print("- label, a label reffering to the network.")
    print("\n")
    
 
def main():

    if(len(sys.argv) != N_ARG):
        requested_argv()
        sys.exit(1)
        
    filename = sys.argv[1]
    label = sys.argv[2]
    
    ### loading network...
    links = np.loadtxt(filename)
    #print(links[:10,:])
    
    
    ### checking network file shape... 
    q,w = np.shape(links)
    
    #### checking number of of rows and columns...
    sys.stderr.write("links--> %d, %d\n"%(np.shape(links)))
    if(np.shape(links)[1] !=3):
        sys.stderr.write("Error! Links file has more that three columns!\n")
        exit()

    ### computing network cerna and mirna numbers...
    m,n = NF.SpeciesNumbers(links)
    
    d1 = STF.A_ai_matrix(links)
    
    d_ = np.zeros((m+n,m+n))
    
    d_[:m,m:] = np.transpose(d1)
    d_[m:,:m] = d1

    ### STF scritp is too slow!
    ### creating mXm contact matrix of 1 miRNA distance...
    #d2 = STF.A_ij_matrix(links)
    #l = STF.ShortestPath_func(links)

    print("Computing shortest paths!")
    l_tmp = sp.shortest_path(d_)

    
    l = np.triu(l_tmp[:m,:m],1)
    
    ### saving shortest path matrix...
    np.savetxt("SP_%s.dat"%label, l.astype(np.int8), fmt = "%i")
    print("The end")
    


main()


