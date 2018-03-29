import numpy as np
import multiprocessing
from primerProbability import *

def de(fobj, seq_len, mut=0.8, crossp=0.7, popsize=20, its=1000, cores=10):
    '''
    Testing an implementation of Differential Evolution optimization algorithm
    '''
    workers = multiprocessing.Pool(cores)
    pop = np.random.rand(popsize, seq_len*4)  ## Makes the population
    # min_b, max_b = np.asarray(bounds).T
    # diff = np.fabs(min_b - max_b)   ## Functions that operate element by element on whole arrays.
    # pop_denorm = min_b + pop * diff
    pop_denorm = np.vstack([el for el in map(reset_costraints, pop)])
    fitness = np.asarray([fobj(ind) for ind in pop_denorm])
    best_idx = np.argmin(fitness)   ## Returns the indices of the minimum values along an axis.
    best = pop_denorm[best_idx]
    for i in range(its):
        for f,j in workers.imap_inordered(de_mutation, [ (fobj,j,popsize,mut,crossp) for j in range(popsize)])
            if f < fitness[j]:
                fitness[j] = f
                pop[j] = trial
            if f < fitness[best_idx]:
                best_idx = j
                best = trial_denorm
    yield best, fitness[best_idx]

def de_mutation(params):
    '''
    Function for mutation of one vector in the population of vectors
    '''
    fobj,j,popsize,mut,crossp = params
    idxs = [idx for idx in range(popsize) if idx != j]
    a, b, c = pop[np.random.choice(idxs, 3, replace = False)]
    mutant = np.clip(a + mut*(b - c), 0, 1)
    cross_points = np.random.rand(seq_len*4) < crossp
    if not np.any(cross_points):    ## If there are no cross_points
        cross_points[np.random.randint(0, dimensions)] = True
    trial = np.where(cross_points, mutant, pop[j])
    trial_denorm = reset_costraints(trial)
    f = fobj(trial_denorm)
    return(f,j)

    if f < fitness[j]:
        fitness[j] = f
        pop[j] = trial
    if f < fitness[best_idx]:
        best_idx = j
        best = trial_denorm

arr=np.asarray(random_ppm(0.01))
vec = arr.reshape([1,24])

arr2 = np.random.rand(6,4).reshape(1,24)
arr2 /= arr2.sum(axis=1)[:,np.newaxis]
arr2.reshape(1,24)

def reset_costraints(prob_vector):
    '''
    Input: vector of 4xn elements in [0,1], where n is the number of positions in the sequence
    Output: vector with costraints respected (the sum of every group of 4 is < 1)
    '''
    # if pop.shape[0]!=1:
    #     return(print('Input needs to be a vector (one row)!!'))
    a = prob_vector.reshape(int(prob_vector.size/4),4)
    a /= a.sum(axis=1)[:,np.newaxis]
    return(a.reshape(1,prob_vector.size))

def test_function(prob):
    '''
    Minimizes the probability of having CGCGCG hexamer
    '''
    hexProb = prob_from_ppm(from_vec_to_ppm(prob), all_hexamers())
    return(hexProb.loc['CGCGCG'][0])
