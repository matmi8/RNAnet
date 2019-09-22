### libraries...
### python libraries
import numpy as np
import matplotlib.pyplot as mpl
import os,sys
import glob

### hand-made libraries
import Network_functions as NF
import Structure_functions as STF
import Susceptibility_functions as SF

DEB_RHO = 0

EPSILON = 1e-8 ### threshold for considering X_ij null..
#### NB: it is chosen as the precision threshold of double floating machine approximation...



#####
##### This file contains functions that compute userful observables in miRNA-ceRNA networks...
##### 



#--> Obs_conc_func(m, n, m_, n_)                                              TBC
#--> Obs_X_MeanMax_func(X_ij)                                                 TBC
#--> Obs_X_Delta_func(X_ij)                                                   TBC
#--> Obs_X_rho_func(X_ij, K_ij)                                               TBC
#--> Obs_Sel_func(X_ij)                                                       TBC
#--> Cv_Func(links, Opt,  save_profile = False, savefile_name = "nowhere")    TBC
#--> CompObs(links, beta, Opt)        TBC
#--> CompObs_Xii(links, beta, Opt)        TBC




def Obs_conc_func(m, n, m_, n_):

    ### computing number of cernas (M) and mirnas (N)...
    M = float(len(m))
    N = float(len(n))

    ### computing m, n tot... 
    m_tot = np.sum(m)
    n_tot = np.sum(n)

    ### comouting m/m_ and n/n_ means...
    mean_m = np.sum(m/m_)/float(M)
    mean_n= np.sum(n/n_)/float(N)

    return(m_tot, n_tot, mean_m, mean_n)


def Obs_X_MeanMax_func(X_ij):

       
    ### computing number of cernas (M) and mirnas (N)...
    #M, N = NF.SpeciesNumbers(links) 
    M = np.shape(X_ij)[0]
    
    X, Y = np.where(X_ij < EPSILON)   

    X_ij[X,Y] = 0

    ### removing diagonal values of X...
    np.fill_diagonal(X_ij, 0)

    ### computing means quantities related to X...
    meanX_ij = np.sum(X_ij)/float(M*(M-1))
    

    ### computing max of off diagonal X...
    maxX_ij = np.max(X_ij)

    return(meanX_ij, maxX_ij)



def Obs_X_Delta_func(X_ij):


    M = np.shape(X_ij)[0]
    
    ### computing Delta (directionality-asimmetry proxy)...

    ### computing X_ij transpose...
    X_ji = np.transpose(X_ij)

    ### computing X_ij +- X_ji...
    dX_plus = (X_ij + X_ji)
    dX_minus = (X_ij - X_ji)

    ### old countin of non nan elements...
    #X_phi, Y_phi = np.where(np.abs(dX_minus) > EPSILON)
    #Phi = len(X_phi)/2.
    #X_phi = []
    #Y_phi = []
    

    Delta_ij = (dX_minus/dX_plus)**2

    ### finding pos where X_ij is null...
    X, Y = np.where(X_ij < EPSILON)   

    Delta_ij[X,Y] = 0
    Delta_ij[Y,X] = 0

    Delta_ij = np.sum(Delta_ij)/float(M*(M-1))

    return(Delta_ij)


def Obs_X_rho_func(X_ij, K_ij):

    M = np.shape(X_ij)[0]
    
    ### to check whether now is still needed the
    ### introduction of the STDX cut-off
    ### (possible error in the previuos code...)
    STDX = 1e-6

    ### removing diagonal values of K in order to erase self terms...
    np.fill_diagonal(K_ij, 0)

    X, Y = np.where(X_ij < EPSILON)   
    X_ij[X,Y] = 0.
    
    ### removing diagonal values of X...
    np.fill_diagonal(X_ij, 0)



    XK = X_ij*K_ij
    
    ### computing means quantities related to K...
    meanK_ij = np.sum(K_ij)/float(M*(M-1))
    meanK2_ij = np.sum(K_ij**2)/float(M*(M-1))
   
    ## mean XK has diagonal values put to zeros thanks to K...
    meanXK_ij = np.sum(XK)/float(M*(M-1)) 

    ### computing means quantities related to X...
    meanX_ij = np.sum(X_ij)/float(M*(M-1))
    meanX2_ij = np.sum(X_ij**2)/float(M*(M-1))

   
    ### computing Pearson coefficient...
    sigmaX = np.sqrt( meanX2_ij- meanX_ij**2)
    sigmaK = np.sqrt( meanK2_ij- meanK_ij**2)

    if(sigmaX < EPSILON):
        sigmaX = EPSILON
        
    Cor = (meanXK_ij -  meanX_ij*meanK_ij)
   
    rho_ij =   Cor/(sigmaX*sigmaK)

    if(DEB_RHO):
        return(rho_ij,meanXK_ij, meanX2_ij, meanK2_ij, meanX_ij, meanK_ij )
    else:
        return(rho_ij)


def Obs_Sel_func(X_ij):

    M = np.shape(X_ij)[0]

    ### computing selectivities...
    g_i = np.sum(X_ij**2, axis = 1)/np.sum(X_ij, axis = 1)**2
    h_j = np.sum(X_ij**2, axis= 0)/np.sum(X_ij, axis = 0)**2

    ### think about the nan = 1 passage...
    g_i[np.isnan(g_i)] = 1
    h_j[np.isnan(h_j)] = 1
    
    
    S_in = np.sum(g_i)/float(M)
    g_i = [] ## freeing memory..
    S_out = np.sum(h_j)/float(M)
    h_j = [] ## freeing memory..

    return(S_in, S_out)




#### this function computes the Cv and SNR trends as a function of beta for a given miRNA-ceRNA network, with pameters given in Opt... 
def Cv_Func(links, Opt,  save_profile = False, savefile_name = "nowhere"):


    N_step = int(Opt.N_step)
    N_rip = int(Opt.N_rip)
    
    
    ### computing number of mirna and cernas...
    M, N = NF.SpeciesNumbers(links) 

    ### initializing beta vector..
    beta_vet = np.logspace(np.log(Opt.betaMin), np.log(Opt.betaMax), Opt.N_step, base= np.e)
    
    
    ### intializing res matrix... 
    Res = np.zeros((N_step, 3))
    
    ### start iterating over betas...
    for i in range(N_step):

        ### initializing beta..
        beta = beta_vet[i]
        
        m_1 = np.zeros(M)
        m_2 = np.zeros(M)

        if(save_profile):
            M_profile = np.zeros((M, N_rip))

        for j in range(N_rip):
    	
	    #### computing m0 vet...
    	    if(Opt.Cv_b == 0):
            	m_ = NF.SampleDistFunc(Opt.meanb, 0, "constant", M)/float(Opt.d) 
            else:     
        	m_ = NF.SampleDistFunc(Opt.meanb, Opt.Cv_b*Opt.meanb,  "lognorm", M)/float(Opt.d)
    
            #### computing mu0 vet...
            if(Opt.Cv_beta == 0):
                n_ = NF.SampleDistFunc(beta, 0, "constant", N)/float(Opt.delta) 
            else:     
                n_ = NF.SampleDistFunc( beta, Opt.Cv_beta*beta, "lognorm", N)/float(Opt.delta)
            
            ### computing mirna and cerna concetrations...
            n,m = NF.ConcFinder(links, n_, m_, Opt.lam) 

            if(save_profile):
                M_profile[:,j] = m
                
            m_1 += m
            m_2 += m**2
        
        if(save_profile):
            np.savetxt(savefile_name, M_profile, fmt= "%.4e", header= "M x N_rip")
        
        m_1 /= N_rip
        m_2 /= N_rip
    
        s_m = np.sqrt(m_2 - m_1**2)
    
        Cv = s_m/m_1
        SNR = m_1/s_m
    
    
        Res[i, 0] = beta
        Res[i,1] = np.mean(Cv)
        Res[i,2] = np.mean(SNR)
    
        

    return(Res)







def CompObs(links, beta, Opt):

    
    ### computing number of mirna and cernas...
    M, N = NF.SpeciesNumbers(links) 

    
    #### computing m0 vet...
    if(Opt.Cv_b == 0):
        m_ = NF.SampleDistFunc(Opt.meanb, 0, "constant", M)/float(Opt.d) 
    else:     
        m_ = NF.SampleDistFunc(Opt.meanb, Opt.Cv_b*Opt.meanb,  "lognorm", M)/float(Opt.d)
        
    #### computing mu0 vet...
    if(Opt.Cv_beta == 0):
        n_ = NF.SampleDistFunc(beta, 0, "constant", N)/float(Opt.delta) 
    else:     
        n_ = NF.SampleDistFunc( beta, Opt.Cv_beta*beta, "lognorm", N)/float(Opt.delta)
        
    ### computing mirna and cerna concetrations...
    n, m = NF.ConcFinder(Links=links, vet_star_m= m_, vet_star_n= n_, lam= Opt.lam)
    
    
    ### computing susceptibility and inverse binding weigth proxy...
    X_ij, K_ij = SF.X_K_ij_func(links, n,m, n_,m_, Opt.lam, 1)
    
    mm, nn, m_norm, n_norm = Obs_conc_func(m, n, m_, n_)
    Xmean, Xmax = Obs_X_MeanMax_func(X_ij)
    Delta = Obs_X_Delta_func(X_ij)
    rho = Obs_X_rho_func(X_ij, K_ij)
    Sin, Sout = Obs_Sel_func(X_ij)
    
    return(m_norm, n_norm, Xmean, Xmax, Delta, rho, Sin, Sout)



def CompObs_Xii(links, beta, Opt):

    
    ### computing number of mirna and cernas...
    M, N = NF.SpeciesNumbers(links) 

    
    #### computing m0 vet...
    if(Opt.Cv_b == 0):
        m_ = NF.SampleDistFunc(Opt.meanb, 0, "constant", M)/float(Opt.d) 
    else:     
        m_ = NF.SampleDistFunc(Opt.meanb, Opt.Cv_b*Opt.meanb,  "lognorm", M)/float(Opt.d)
        
    #### computing mu0 vet...
    if(Opt.Cv_beta == 0):
        n_ = NF.SampleDistFunc(beta, 0, "constant", N)/float(Opt.delta) 
    else:     
        n_ = NF.SampleDistFunc( beta, Opt.Cv_beta*beta, "lognorm", N)/float(Opt.delta)
        
    ### computing mirna and cerna concetrations...
    n, m = NF.ConcFinder(Links=links, vet_star_m= m_, vet_star_n= n_, lam= Opt.lam)
    
    
    ### computing susceptibility and inverse binding weigth proxy...
    X_ij, K_ij = SF.X_K_ij_func(links, n,m, n_,m_, Opt.lam, 1)
    
    Xii_mean = np.mean(np.diag(X_ij))
    Xii_max = np.max(np.diag(X_ij))
     
    return(Xii_mean, Xii_max)

