### libraries...
### python libraries
import numpy as np
import matplotlib.pyplot as mpl
import os,sys
import glob

#####
##### This file contains functions that process a given miRNA-ceRNA network
#####

#--> LinksFunc(filename, m0, Bd)
#--> SpeciesNumbers(links)
#--> SampleDistFunc(m, s , choice, N)
#--> ConcFinder(Links, vet_star_n, vet_star_m, lam)      
#--> f_x(x, A, a)                                        



### defining class Options with all parameters read from a file...
class Options:
    def __init__(self, filename):
        self.m0, self.d, self.delta, self.lambda_, self.meanb, self.betaMin, self.betaMax, self.Cv_beta, self.Cv_b, self.Bd, self.N_rip, self.N_step =  np.loadtxt(filename, usecols=[2])
        ### N.B. lam is  redefined as the product of the participation ratio
        ### lambda and the ratio d/delta
        ### So that we avoid to modify all other fucntions in the scripts
        self.lam = (1.*self.d/self.delta)*self.lambda_
        
        #self.sigmab = self.Cv_b*self.meanb        
    
    




### The function reads and properly modified df_links...
### returning a network of the form:
###
### cerna mirna m0*f with f = 0.3,0.4, 0.6, 1.  corresponding to
### none, none9nt, nc or kmer binding modes. 

def LinksFunc(filename, m0, Bd):
    
    links = np.loadtxt(filename) #, delimiter = ",")
    
    #### checking number of of rows and columns...
    print("links-->",np.shape(links))
    if(np.shape(links)[1] !=3):
        print("Error! Links file has more that three columns!")
        exit()
        
    links = abs(links)

    ### Assigning weights to links...
 
    if(Bd == 0):
        print("delta in mu0")
        links[:,2] = m0*2.
    elif(Bd == 1):
        print("Bimodal heterogeneity")
        for i in range(0,np.shape(links)[0]):
            #prob = np.random.random()
            seed_type = links[i,2]
            if(seed_type == 3 or seed_type == 4 or seed_type == 5 or seed_type == 6):
	        links[i,2] =  m0                 #m0
            else:
                links[i,2] =  m0*2.              #m0/0.3 
    elif(Bd == 2):
        print("Real heterogeneity")
        for i in range(0,np.shape(links)[0]):
            seed_type = links[i,2]
            if(seed_type == 5 or seed_type == 6):
	        links[i,2] =  m0/2.          #m0
            elif(seed_type == 4):
                links[i,2] =  m0             #m0/0.6
            elif(seed_type == 2 or seed_type == 3):
                links[i,2] =  m0*2.          #m0/0.4
            else:
                links[i,2] =  m0*4.          #m0/0.3
                
    return links


#### this function returns the number of cernas (M) and of mirnas (N) taking in input the links network
def SpeciesNumbers(links):
    
    M = np.shape(np.unique(links[:,0]))[0]  ### m --> genes
    N = np.shape(np.unique(links[:,1]))[0]  ### n --> mirna

    return(M,N)


#### this fucntion returns a vector of length N of number sampled from the "choice" distribution:
####  "lognorm" for a lognormal distribution with mean m and sigma s
####  "norm" for a normal distribution with mean m and sigma s
####  "uniform" for a uniform sampling in  m+-s
####  "constant" for a vector of m values

def SampleDistFunc(m, s , choice, N):
    
    
    if(choice == "lognorm"):
    
        m_eff = np.log(m**2/np.sqrt(s**2 + m**2))
        s_eff = np.sqrt(np.log(s**2/m**2+1))

        n_ = np.random.lognormal(m_eff, s_eff,N)
    elif(choice == "norm"):
        n_ = np.random.normal(m, s,N)
    elif(choice == "constant"):
        n_ = np.ones(N)*m
    elif(choice == "uniform"):
	a = m*(1- np.sqrt(3)*s/m)
	b = m*(1+ np.sqrt(3)*s/m)
	if(a < 0):
		a = 0
		b = 2*m
		cv = 1/np.sqrt(3)
		print("CV is too large!\n")
		print("Max CV is %.3e\n"%cv)
	
        n_ = a +  np.random.random(N)*(b-a)       
    else:
        print("Error, choice of distribution not allowed!")
    
    return n_


### This function computes the concentrations of mirna and mrna, given their network (links) and their maximal concentrations vet star m and n and lam;
def ConcFinder(Links, vet_star_n, vet_star_m, lam):

    PRECISION_THRES = 1e-9
    
    links = np.copy(Links)

    linkNum = np.shape(links)[0]
    
    m,n = SpeciesNumbers(links)

    vet_m = np.ones(m)
    vet_n = np.ones(n)

    ### The mirna matri...
    M = np.zeros((m,n))

    ### The mrna matrix...
    N = np.zeros((n,m))
    
    for i in range(0,linkNum):

        c = links[i,2] ### weight
        
        M[int(links[i,0]),int(links[i,1])] = 1./(c*lam)  # mu
        N[int(links[i,1]), int(links[i,0])] = 1./c       # m

    chi = [1]
    
    while(abs(chi[-1]) > PRECISION_THRES):

        tmp_n = vet_n
        tmp_m = vet_m
        
        vet_n = f_x(vet_m, N, vet_star_n)
        vet_m = f_x(vet_n, M, vet_star_m)

        chi__ = abs(np.sum((abs(vet_n)-abs(tmp_n)))) + abs(np.sum((abs(vet_m)-abs(tmp_m))))
        chi.append(chi__)

    return (vet_n, vet_m)

#### this fucntion performs a operation on matrix required for ConcFinder..
def f_x(x, A, a):
    tmp = np.dot(A,x)
    tmp = tmp + 1
    return a/tmp
