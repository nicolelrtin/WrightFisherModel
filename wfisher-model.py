#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
sns.set(style="whitegrid")
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


#given an array, find ID of mutants
def find_mutants(array, cutoff):
    def find(x):
          if x < cutoff:
            return True
          else:
            return False
    mutants = filter(find, array)
    
    thing = np.array([])
    for x in mutants:
        thing = np.append(thing, x)
    return thing

#are all the cells greater than the cutoff value?
# in other words, are all the cells wt?
# if array still contains mutants, add another day
def is_extinct(array, cutoff):
    hh = all(x >= cutoff for x in array)
    if hh: val = 0
    else: val = 1
    return val


# In[11]:


# pmutants, dilfactor are perctentages

def wfisher(N, nflasks, pmutants, gens, dilfactor):
    
    cutoff = int(N*pmutants)
    pick = int(N*dilfactor)
    
    #initialize ancestral population, give all non-unique ID
    IDs = np.random.choice(np.arange(N*1.2), N, replace = False)

    #total count of mutants in each gen (for plotting)
    mutants = np.array([])
    
    #population census (for records & manipulation)
    dils = np.random.choice(IDs, pick)

             
   #count of mutants at t =0 (to be appended to mutants)
    count = len(find_mutants(dils, cutoff))
    #split ancestral population into 9+1 flasks
    for i in np.arange(nflasks-1):
        select = np.random.choice(IDs, pick)
        dils = np.vstack([dils, select])
        count = np.append(count, len(find_mutants(select, cutoff)))
    
    #Initialize flasks for all generations. This is now a 3-tensor (all_ijk denotes the ith generation, the jth flask and the kth cell)
    all = np.array([dils])    
    #print(all)
    
    mutants = np.append(mutants, count)

    #propagate cultures gen number of times
    for i in np.arange(gens): 
        #print('Generation %i ' %(i+1))
        #initialize first flask; scale up most recent dilution
        lastdil = np.tile(all[-1], int(N/pick))
        #choose selection
        gendil = np.random.choice(lastdil[0], pick)
        #find num of mutants
        genmutants = len(find_mutants(gendil, cutoff))
        for j in np.arange(nflasks-1):
            btlneck = np.random.choice(lastdil[j+1], pick)
            gendil = np.vstack([gendil, btlneck])
            genmutants = np.append(genmutants, len(find_mutants(btlneck, cutoff)))
        all = np.vstack([all,[gendil]])
        mutants = np.vstack([mutants, genmutants])
        
    #binomial distribution
    binom = np.array([count])
    for i in range(gens):
        p = binom[-1]/pick
        z = np.random.binomial(pick, p)
        binom = np.vstack([binom, z])
    
    #poisson dist
    poi = np.array([count])
    for i in range(gens):
        lam = np.random.poisson(poi[-1])
        for (i, item) in enumerate(lam):
            if item > N:
                lam[i] = N
        poi = np.vstack([poi, lam])
    
    return N, nflasks, pmutants, gens, dilfactor, mutants, binom, poi


#wfisher(N, nflasks, pmutants, gens, dilfactor)
np.random.seed()
# wfisher(100000, 10, 0.01, 300, 0.01)

N, nflasks, pmutants, gens, dilfactor, mutants, binom, poi = wfisher(100000, 50, 0.01, 250, 0.01)


# In[12]:


def STATS(N, nflasks, pmutants, gens, dilfactor, mutants, binom, poi):
    #extinction
    mext= np.array([])
    bext= np.array([])
    pext= np.array([])
    N, nflasks, pmutants, gens, dilfactor, mutants, binom, poi = wfisher(N, nflasks, pmutants, gens, dilfactor)
    for j in range(nflasks):
        mext = np.append(mext, np.count_nonzero(mutants[:,j]))
        bext = np.append(bext, np.count_nonzero(binom[:,j]))
        pext = np.append(pext, np.count_nonzero(poi[:,j]))
    print('mean extinction time for simulated values:',np.mean(mext),'binomial dist:',np.mean(bext),'poisson dist:', np.mean(pext))

    #VISUALIZE DATA
    pick = N*dilfactor
    fig = plt.figure(figsize=(11,7)) #initialize a figure
    fig.subplots_adjust(wspace=.7,hspace=0.4)
    plt.subplot(111)
    xrange = range(gens+1)
    for j in range(nflasks):
        plt.plot(xrange, mutants[:,j]/pick, ':r', linewidth=1)
        plt.plot(xrange, binom[:,j]/pick, ':', c='skyblue', linewidth=1)
        plt.plot(xrange, poi[:,j]/pick, ':', c='orchid', linewidth=1)
    #plot MEAN f
    plt.plot(xrange, np.mean(mutants/pick,axis=1), '-r', linewidth=1)
    plt.plot(xrange, np.mean(binom/pick,axis=1), c='skyblue', linewidth=1)
    plt.plot(xrange, np.mean(poi/pick,axis=1), 'orchid', linewidth=1)
    plt.title('Mutant frequency')
    plt.xlabel('transfer number')
    plt.ylabel('mutant frequency')
    plt.show()
    
    #plot MEAN and VARIANCE n
    # mean should be constant and variance changing
    fig = plt.figure(figsize=(11,7)) #initialize a figure
    fig.subplots_adjust(wspace=.7,hspace=0.4)
    plt.subplot(111)
    xrange = range(gens+1)
    plt.plot(xrange, np.mean(mutants,axis=1), '-r', linewidth=1)
    plt.plot(xrange, np.mean(binom,axis=1), 'skyblue', linewidth=1)
    plt.plot(xrange, np.mean(poi,axis=1), 'orchid', linewidth=1)
    
    plt.plot(xrange, np.sqrt(np.var(mutants,axis=1)), ':r', linewidth=1)
    plt.plot(xrange, np.sqrt(np.var(binom,axis=1)), ':', c='skyblue', linewidth=1)
    plt.plot(xrange, np.sqrt(np.var(poi,axis=1)),':', c= 'orchid', linewidth=1)
    plt.title('Mean and Variance (Count)')
    plt.xlabel('transfer number')
    plt.ylabel('mutant variance')
    plt.show() 
    
    # DRAW WITH NEGATIVE BINOMIAL
    t=10 #time points
    r=nflasks #replicates
    phi=50 #dispersion factor
    sim=np.zeros((t,r))
    sim[0,:]= mutants[0]

    for j in range(1,t):
        p=phi/(phi+sim[j-1,:])
        sim[j,:]=np.random.negative_binomial(phi,p)
        
    fig = plt.figure(figsize=(11,7))
    fig.subplots_adjust(wspace=.7,hspace=0.4)
    plt.subplot(111)
    for i in range(r):
        plt.plot(range(t), sim[:,i])
    plt.plot(range(t), np.mean(sim, axis=1))
    plt.title('Trajectory of Mutants using a Negative Binomial Model')
    plt.xlabel('transfer number')
    plt.ylabel('mutant variance')
    plt.show()

STATS(N, nflasks, pmutants, gens, dilfactor, mutants, binom, poi)


#error: 'division by zero'
def increaseN(nflasks, pmutants, gens, dilfactor):
    tens = [30, 75, 150, 300, 600, 1200, 2400, 5000]
    pick = []
    for ten in tens:
        picked = ten*dilfactor
        pick.append(picked)
        
    #plot each loop
    fig = plt.figure(figsize=(11,7)) #initialize a figure
    fig.subplots_adjust(wspace=.7,hspace=0.4)
    plt.subplot(111)
    xrange = range(gens+1)
    for n in tens:
        Nn, nf, pm, gen, dilf, m, b, p = wfisher(n, nflasks, pmutants, gens, dilfactor)
        for j in range(nflasks):
            plt.plot(xrange, m[:,j], ':r', linewidth=1)
        plt.plot(xrange, np.mean(m,axis=1), '-r', linewidth=1)
    plt.show()

# increaseN(10, 0.01, 150, 0.01)


# In[ ]:




