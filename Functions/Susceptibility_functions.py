### libraries...
### python libraries
import numpy as np
import matplotlib.pyplot as mpl
import os,sys
import glob

### hand-made libraries
import Network_functions as NF


TIME = 0 ### showing progress in computing susceptibilities...


#####
##### This file contains functions that compute susceptibility-related stuff...
##### 


#--> X_K_ij_func(links, n,m, n_,m_, lam, choiceK = 1)            TBC
#--> X_K_ab_func(links, n,m, n_,m_, lam, choiceK = 1)            TBC
#--> X_ia_func(links, X_ab, m, m_)                               TBC
#--> W_K_func(links, n,m, n_,m_, lam, M, N, choiceK = 1)         TBC



##### this function compute the X_ij susceptibility matrix and the K_ij matrix quantifing the near-neigbor interaction strengh
##### K_ij computation can be omitted by setting choiceK = 0

def X_K_ij_func(links, n,m, n_,m_, lam, choiceK = 1):

    M, N = NF.SpeciesNumbers(links) ### cerna (M) and mirna (N)

    if(choiceK):
        W_ij, K_ij = W_K_func(links, n,m, n_,m_, lam, M, N, 1)
    else:
        W_ij = W_K_func(links, n,m, n_,m_, lam, M, N, 0)

    I = np.zeros((M,M))
    np.fill_diagonal(I,1)

    ### computing X = (1-W)^-1 diag(m/m^*)
    X_ij = (np.dot(np.linalg.inv(I-W_ij), np.diag(m/m_)))  
    I= []

    if(choiceK):
        return(X_ij, K_ij)
    else:
        return(X_ij)

def X_K_ab_func(links, n,m, n_,m_, lam, choiceK = 1):
    
    tmp_links = np.zeros(np.shape(links))

    tmp_links[:,0] = links[:,1]
    tmp_links[:,1] = links[:,0]
    tmp_links[:,2] = links[:,2]

    ###### in order to compute X_ab we can use X_ij func switching both the links file columns and the concetrations (m,n, m_,n_)
    if(choiceK):
        X_ab,K_ab = X_K_ij_func(tmp_links, m,n, m_,n_, lam, 1)
        return(X_ab,K_ab)
    else:
        X_ab = X_K_ij_func(tmp_links, m,n, m_,n_, lam, 0)
        return(X_ab)
    

#### this function computes the X_ia susceptibility matrix...
def X_ia_func(links, X_ab, m, m_):
    
    # computing number of mirnas and cernas...
    M, N = NF.SpeciesNumbers(links) 


    # creating distance matrix... 
    d = np.zeros((M,N))

    x = (links[:,0]).astype(int)
    y = (links[:,1]).astype(int)
    z = (links[:,2]).astype(float)
    
    d[x,y] = 1./z
    
    for i in range(M):
        d[i,:] = d[i,:]*-1*m[i]**2/m_[i]
    
    
    # initialing X_ia matrix...
    X_ia = np.zeros((M,N))
    
    for i in range(M):
        X_ia[i,:] = np.dot(np.transpose(X_ab), d[i,:])
    
    return(X_ia)



def W_K_func(links, n,m, n_,m_, lam, M, N, choiceK = 1):
    
    #W_ij = np.zeros((M,M), dtype=np.float32)
    W_ij = np.zeros((M,M), dtype=np.float64)

    if(choiceK):
        K_ij = np.zeros((M,M),dtype=np.float32)
        

    #### We had to compute W_ij = m_i^2/m_i^* sum_alpha mu_a^2/mu_a^* * 1/(m_ja0*mu_ia0)

    ### the W can be computed faster in the case of big X matrix, computing a row at times,as W_ix = m_i^2/m_i^* dot(J, vet(mu_a^2/mu_a^*/mu_ia))
    #### where J is a matrix of dim (M,N_ai) and N_ai is the number of mirnas that contact cerna i. J elements are 0 or 1/(m_jam) if cerna j contacts cerna i through mirna alpha.
    
    #Computing W_ij... 
    for i in range(M):
        if(TIME):
            sys.stderr.write("\rSimulation in progress : %d" % int(i*100./M))
	    sys.stderr.flush()
            
        m_i = m[i] # m_i
        m_i_ = m_[i] # m_i*

        ### remebering that mu_ia = lambda*m_ia ...
        ### N.B. links[:,2] = m_0!
        
        mu_ia = links[links[:,0] == i,2]*lam # mu_ia

        if(choiceK == 1):

            ### mu_ia = d/k^+(1+ k^-/(sigma-kappa)), if k^- = 0
            ### mu_ia = d/k^+
            
            
            k_ia = 1./mu_ia # k_ia/d
            ### the d term symplifies in the rho calculation so it can be ignored!
            #k_ia = 1./links[links[:,0] == i,2] # k_ia
 
        ### Selecting mirnas that bind cerna i...
        a = links[links[:,0] == i,1]

        
        dim_a = np.shape(a)[0]

        ### creating matrix of all possible cerna j (row)  and mirnas a (col), so that we compute in one shot all W_ix with x in (1,M)..  
        J = np.zeros((M, dim_a))

        for q in range(dim_a):

            ### selecting cernas j that contact cerna i througth mirna alpha
            j = links[links[:,1] == a[q],0]
            m_ja = links[links[:,1]== a[q],2]
	    J[j.astype(int),q] = 1/m_ja

                
        ### Creating the vector of mirnas concentrations mu and mu^*
        m_a = n[a.astype(int)]
        m_a_= n_[a.astype(int)]
       
        
        J_a = np.dot(J, (m_a**2/m_a_)/mu_ia)

        W_ij[i,:] = J_a*m_i**2/m_i_
        
        if(choiceK == 1):
            K_ij[i,:] = np.dot(J/lam,k_ia)/float(N)

    if(choiceK):
        return(W_ij, K_ij)
    else:
        return(W_ij)

################ end W_K_func #############################
