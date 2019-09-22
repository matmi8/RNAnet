
#coding: utf-8

import numpy as np
import matplotlib.pyplot as mpl
import os, sys

sys.path.append("../Functions/")
import Network_functions as NF
import Susceptibility_functions as SF
import Structure_functions as STF
import Observables_functions as OF


if(len(sys.argv) !=3):
	print("Error! Invalid number of arguments!")
	print("")
	print("Insert:")
	print("python EdgeSwapGen.py linkfile savefile")
	print("- linkfile, the initial network file.")
	print("- savefile, the file where to save the edgeswapped network.")
	print("")
	exit()

linkname = sys.argv[1]
savename = sys.argv[2]

links = np.loadtxt(linkname)

M, N = NF.SpeciesNumbers(links=links)

print(M,N)
print("")
print("Intial network is connected? %r"%np.bool(STF.ConnectivityCheck(links)))
print("")



def EdgeSwap(links):
   
    new_links = links.copy()
    
    N = np.shape(new_links)[0]
    
    a = np.random.randint(0,N)
    b = np.random.randint(0,N)
    
    while(b == a):
         b = np.random.randint(0,N)
    
    tmp = new_links[a, 1]  
    new_links[a,1] = new_links[b,1]
    new_links[b,1] = tmp
    
    check_c = STF.ConnectivityCheck(new_links)
    check_l = LinkNumberCheck(links,new_links)
    
    if(check_c and check_l):
        return(new_links,1)
    else:
        return(links, 0)

    
def NewtorkSimilarity(links,new_links):
    a = links[:,0] + 1j*links[:,1]
    b = new_links[:,0] + 1j*new_links[:,1]
    
    Na = len(a)
    Nb = len(b)
    
    if(Na != Nb):
        print("Attention! Networks have different numbers of links!")
    
    s = len(np.intersect1d(a,b))
    
    return(s/float(Na))
    

def LinkNumberCheck(links,new_links):
    
    a = links[:,0] + 1j*links[:,1]
    b = new_links[:,0] + 1j*new_links[:,1]
    
    N = len(np.unique(a))
    M = len(np.unique(b))
    
    if(N == M):
        return 1
    else:
        return 0



## coping original network link 
ll = links.copy()

ss = []
cc = []

Nstep = 100000


for i in range(Nstep):
    
    sys.stderr.write("\r Progress %.1f %%"%(i*100./Nstep))
    sys.stderr.flush()
    
    l,c= EdgeSwap(ll)
    
    cc.append(c)
    
    ll = l.copy()



print( "Final network similarity",NewtorkSimilarity(ll, links))


np.savetxt("%s_ES.dat"%savename, ll,fmt="%d")

