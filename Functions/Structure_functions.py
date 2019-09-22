import os, sys
import numpy as np

from scipy.optimize import fsolve
import matplotlib.pyplot as mpl
import random as rm



import Network_functions as NF

PATH_PRINT = 0
### a big value used as indicator of uncomputed path...
unset_value = 100


#####
##### This file contains functions that build and randomize miRNA-ceRNA networks...
##### 

#-->  ShortestPath_func(links)      TBC
#-->  A_ij_matrix(links)
#-->  A_ai_matrix(links)
#-->  ConnectivityCheck(Links)      TBC
#-->  ER_rand(links)                TBC


### This function computes the shortest paths between first group of nodes in a bipartite network...
### N.B. the function works but the scipy one is better up to CLASH networks (TargetScan bho!).
def ShortestPath_func(links):

    ### computing network cerna and mirna numbers...
    m,n = NF.SpeciesNumbers(links)

    ### creating mXm contact matrix of 1 miRNA distance...
    d2 = A_ij_matrix(links)
    
    
    ### freeing memory...
    #links = [] 
    
    ### filling uncomputed paths with big value (bigger that paths lenghts)...
    print(np.sum(d2))
    
    d2[np.where(d2==0)] = unset_value

    print(np.sum(d2))
    
    ### excluding lower matrix triangle from calculations...
    d2 = np.triu(d2, 0).astype(np.int8)
    
    
    ### finding d=2 paths positions...
    x2, y2 = np.where(np.triu(d2,0) == 2)

    print(len(x2))
    ### turning positions in complex number to better use them...
    c2 = x2 + 1j*y2
    
    v2, w2 = np.where(np.triu(d2,0) == unset_value)
    print(v2)
    void = v2 + 1j*w2
    
    n_void = len(void)
    print("len void -->%d"%n_void)
    
    d = 2
    
    l = len(c2)
    frac = int(l/4.)

    ### saving d=2 paths...
    if(PATH_PRINT): 
     	np.savetxt("%s_path_%d.dat"%(label,d),c2)
        
    ### a  way in order not to have big vectors...
    ##  we split in four parts...
    if(l>4):
        cd1 = np.copy(c2[:frac])
        cd2 = np.copy(c2[frac:2*frac])
        cd3 = np.copy(c2[2*frac:3*frac])
        cd4 = np.copy(c2[3*frac:])
    else:
        cd = np.copy(c2)
     
    while(n_void > 0):
        
        d +=2
        cd_tmp = []
        sys.stderr.write("Computing %d-paths...\n"%d)
        sys.stderr.flush()
        if(l > 4):
            ### cycling over cernas...
            for i in range(m):
                sys.stdout.write("\rProgress %d %%"%(i*100./m))
                sys.stderr.flush()
        
                ### finding d long paths from d-2 and 2 paths...
                ## the idea: for each i we isolate 2-path (a,i) or (i,a) and (d-2)-path (i,b)  or (b,i)
                ## then we merge them to obtain (a,b) or (b,a) that are the path of length d.
                    
                x = np.append(np.real(c2[np.where(np.imag(c2)==i)]), np.imag(c2[np.where(np.real(c2)==i)]))

          
                y1 = np.append(np.real(cd1[np.where(np.imag(cd1)==i)]), np.imag(cd1[np.where(np.real(cd1)==i)]))
                y2 = np.append(np.real(cd2[np.where(np.imag(cd2)==i)]), np.imag(cd2[np.where(np.real(cd2)==i)]))
                y3 = np.append(np.real(cd3[np.where(np.imag(cd3)==i)]), np.imag(cd3[np.where(np.real(cd3)==i)]))
                y4 = np.append(np.real(cd4[np.where(np.imag(cd4)==i)]), np.imag(cd4[np.where(np.real(cd4)==i)]))

                y = np.append(y1,y2)
                y = np.append(y,y3)
                y = np.append(y,y4)

                ### creating a submatrix of all found positions...
                tmp = (d2[np.ix_(x.astype(int),y.astype(int))] == unset_value).astype(np.int8)

                ### inserting d-paths in d2 matrix...
                d2[np.ix_(x.astype(int),y.astype(int))] =d2[np.ix_(x.astype(int),y.astype(int))] - (tmp)*(unset_value-d)
                    
        else:
            for i in range(m):
                
                sys.stdout.write("\rProgress %d %%"%(i*100./m))
                sys.stderr.flush()
        
                    
                x = np.append(np.real(c2[np.where(np.imag(c2)==i)]), np.imag(c2[np.where(np.real(c2)==i)]))
                y = np.append(np.real(cd[np.where(np.imag(cd)==i)]), np.imag(cd[np.where(np.real(cd)==i)]))

                tmp = (d2[np.ix_(x.astype(int),y.astype(int))] == unset_value).astype(np.int8)
                d2[np.ix_(x.astype(int),y.astype(int))] =d2[np.ix_(x.astype(int),y.astype(int))] - (tmp)*(unset_value-d)


        ### creating vector of all d path positions..
        x,y = np.where(d2 == d)
        cd = x+1j*y
        l = len(cd)
        frac = int(l/4.)

        ### saving it...
	if(PATH_PRINT):
            np.savetxt("%s_path_%d.dat"%(label,d), cd)

     
        if(l > 4):
            
            cd1 = np.copy(cd[:frac])
            cd2 = np.copy(cd[frac:2*frac])
            cd3 = np.copy(cd[2*frac:3*frac])
            cd4 = np.copy(cd[3*frac:])
          
        v2, w2 = np.where((d2) == unset_value) 
        void = v2 + 1j*w2
        n_void_tmp = len(void)

        ### checking whether there are still paths to be found...
        if(n_void == n_void_tmp):
	    break
        n_void = n_void_tmp

    sys.stdout.write("Done!\n")

    return(d2)
########### end of shortest-path function


### This function computes the MxM adiacence matrix A_ij,
### that has 2 (0) at position i,j if ceRNA i and ceRNA j are (are not) linked by a path of 1 miRNAs
def A_ij_matrix(links):

    ### computing network cerna and mirna numbers...
    m,n = NF.SpeciesNumbers(links)

    ### creating nXm contact matrix...
    d0 = A_ai_matrix(links)


    ### initializing cerna-cerna contact matrix through one mirna..
    ## (distance between cerna i and j in terms of mirnas = 0)
    d2 = np.zeros((m,m))
    
    
    ### computing d2 matrix elements...
    for i in range(m):
        ### selecting submatrix of all mirna (n') that contact cerna i, tmp --> n'Xm
        tmp = d0[d0[:,i] == 1]
        ### checking which cernas contact cerna i through one mirna...
        vet = np.sum(tmp,axis=0)
        d2[i,:] = vet	
        d2[:,i] = vet
        
    ### freeing memory...
    d0 = []
    
    
    ### all cells with a value higher than zeros are in contact at distace 2 (1 mirna)
    d2 = (d2 > 0).astype(np.int8)*2
    
    ### filling diagonal with ones in order to exclude it from further calculations...
    np.fill_diagonal(d2,1)

    return(d2)

     
    
### This function computes the rectangular adiacence matrix A_ai,
### that has 1 (0) at position a,i if miRNA a and ceRNA i are (are not) linked 
def A_ai_matrix(links):

     ### computing network cerna and mirna numbers...
     m,n = NF.SpeciesNumbers(links)

     ##  creating contacts position lists... 
     y = list((links[:,0]).astype(int))
     x = list((links[:,1]).astype(int))

     ### initializing mirna-cerna contacts matrix...
     d0 = np.zeros((n,m))

     ### filling contact matrix...
     d0[(x,y)] = 1

     ### freeing memory...
     x = []
     y = []

     return(d0)




### This function checks the connectivity of the graph
def ConnectivityCheck(Links):

    ### creating a copy of the links matrix
    links = np.copy(Links)

    ### initializing inital "contagion" node...
    start_node = links[0,0]

    ### initializing col to the column opposite to the one of start_node... 
    col = 1

    ### finding nodes in contact with start_node...
    nodes = links[links[:,np.abs(col-1)] == start_node, col]

    ### erasing visited links...
    new_links = links[links[:,np.abs(col-1)] != start_node]
    links = new_links

    ### computing number of contacted nodes..
    d = len(nodes)
   
    
    ### starting iteration, jumping on one colunm to the other and eliminating the reached links...
    while(d != 0):

        ###finding all nodes contacted by "nodes"...
        tmp_nodes = []
        for i in range(d):
            tmp_nodes= np.concatenate((links[links[:,col] == nodes[i], np.abs(col-1)], tmp_nodes))
            new_links =   links[links[:,col] != nodes[i]]
            links = new_links

        ### switcing column!
        if(col == 1):
            col = 0
        elif(col == 0):
            col = 1

        ### uploading nodes...
        nodes = tmp_nodes

        ### uploading number of contacted nodes...
        d = len(nodes)

    ### once d is zero, we check if the dimension of the remaining links matrix is zero (network is connected) or not (it is present a disjoined component!
    if(np.shape(links)[0] == 0):
        ## network is fully connected!
        return 1
    else:
        ## network is disjoined!
        return 0


### the function takes a network links matrix and returns a randomized version of it where only number of nodes and numbers of links are preserved (Erdos-Renyi randomization)
def ER_rand(links):

    ### finding number of cernas and mirnas...
    m,n = NF.SpeciesNumbers(links)

    ### finding number of links..
    d = np.shape(links)[0]

    ### creating a copy of network link table...
    ### it will be needed for preserving the weight column..
    new_link = np.copy(links)


    ### initializing the number of unique rows
    l = 1

    ### initializing new 
    s = []

    ### extracting m (number of cernas) mirna nodes (allowing repetitions)
    N = np.random.randint(0,n,m)
    ### creating a vector of all cerna nodes...
    M = np.arange(0,m,1)

    ### creating an initial link file where each cerna node
    ### is linked to one (random) mirna node
    c = M +1j*N
    s = np.append(c,s)

    ### repeating the above procedure with switched species..
    M = np.random.randint(0,m,n)
    N = np.arange(0,n,1)

    c = M +1j*N
    s = np.append(s,c)

    #### The above procedure assures that network is fully connected, since all mirnas and cernas are linked to one rna of the other group at least one time.


    #### initializing the number of rows of the link table needed to reach the desidered d dimension..
    tmp_d = d - m - n


    ### starting iteration...
    while(l != d):

        ### generating random links...
        x = np.random.randint(0,m,tmp_d)
        y = np.random.randint(0,n,tmp_d)

        ### creating tmp table...
        c = x+1j*y

        ### adding tmp table to complete table...
        s = np.append(c,s)

        ### parsing table eliminating identical rows..
        s = np.unique(s)

        ### recomputing number of different rows...
        l = len(s)

        ### uploading number of row left to add...
        tmp_d = (d-l)

    
    ### creating final random link matrix..
    new_link[:,0] = np.real(s)
    new_link[:,1] = np.imag(s)


    ### finding number of cernas and mirnas...
    M, N = NF.SpeciesNumbers(new_link)

    ### checking whether the ER network is indeed connected! 
    check = ConnectivityCheck(new_link) 
    
    
    if( (M-m) == 0 and (N-n) == 0 and check == 1):
        return(new_link)
    else:
        print("Error! Network has problems!")
        print("Mirna number: %d -- cerna number: %d -- connectivity: %d"%((M-m), (N-n), check))
        return(ER_rand(links))
        
    
