import numpy as np
import multiprocessing
import json
import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/optimization')
from primerProbability import *
from predictCovBs import *

def de(fobj, fun_params, seq_len, mut=0.8, crossp=0.7, popsize=20, its=1000, cores=10):
    '''
    Testing an implementation of Differential Evolution optimization algorithm
    '''
    workers = multiprocessing.Pool(cores)
    seqs,dgMat,genomeAb,target = fun_params
    pop = np.random.rand(popsize, seq_len*4)  ## Makes the population
    pop_denorm = np.vstack([el for el in map(reset_costraints, pop)])
    fitness = np.asarray(workers.map(fobj, [(ind, seqs,dgMat,genomeAb, target) for ind in pop_denorm])) # <-- takes too long
    best_idx = np.argmin(fitness)   ## Returns the indices of the minimum values along an axis.
    best = pop_denorm[best_idx]
    performanceVal = []
    performanceMat = []
    for i in range(its):
        print("--- Iteration no. "+ str(i)+" ---", flush=True)
        for f,trial,j in workers.imap_unordered(de_mutation, [ (fobj,j,pop_denorm,fun_params,seq_len,popsize,mut,crossp) for j in range(popsize)]):
            if f < fitness[j]:
                fitness[j] = f
                pop_denorm[j] = trial
            if f < fitness[best_idx]:
                best_idx = j
                best = trial
                print("Best score: ", flush=True)
                print(round(fitness[best_idx],6), flush=True)
                print(from_vec_to_ppm(best), flush=True)
        performanceVal.append(fitness[best_idx])
        performanceMat.append(best)
    yield performanceVal, np.asarray(performanceMat)

def de_mutation(params):
    '''
    Function for mutation of one vector in the population of vectors:
    selects 3 vectors from the population and constructs the mutated vector from them.
    '''
    fobj,j,pop,fun_params,seq_len,popsize,mut,crossp = params
    seqs,dgMat,genomeAb,target = fun_params
    idxs = [idx for idx in range(popsize) if idx != j]
    a, b, c = pop[np.random.choice(idxs, 3, replace = False)]
    mutant = np.clip(a + mut*(b - c), 0, 1)
    cross_points = np.random.rand(seq_len*4) < crossp
    if not np.any(cross_points):    ## If there are no cross_points
        cross_points[np.random.randint(0, dimensions)] = True
    trial = np.where(cross_points, mutant, pop[j])
    trial_denorm = reset_costraints(trial)
    f = fobj((trial_denorm, seqs,dgMat,genomeAb, target))
    return(f,trial_denorm,j)

def reset_costraints(prob_vector):
    '''
    Input: vector of 4xn elements in [0,1], where n is the number of positions in the sequence
    Output: vector with costraints respected (the sum of every group of 4 is < 1)
    '''
    # if pop.shape[0]!=1:
    #     return(print('Input needs to be a vector (one row)!!'))
    p_vec=prob_vector.copy()
    a = p_vec.reshape(4, int(p_vec.size/4))
    a /= a.sum(axis=0)[np.newaxis, :]
    return(a.reshape(1,p_vec.size))

def new_test_function(params):
    '''
    Maximizes the probability of having CGCGCG hexamer
    '''
    prob_vec,seqs,dgMat,genomeAb,genomeAb = params
    ppm = from_vec_to_ppm(prob_vec)
    hexProb = prob_from_ppm(ppm, all_hexamers())
    return(1-hexProb.loc['CGCGCG'][0])

def coverage_function(params):
    '''
    The function to be minimized. May the force be with us.
    Target needs to be a series!
    '''
    prob_vec, seqs, dgMat, genomeAb, target = params
    ppm = from_vec_to_ppm(prob_vec)
    primer_prob = prob_from_ppm(ppm, seqs)
    coverage = predictCoverage_setProbs(dgMat, genomeAb, primer_prob)
    coverage.index=coverage.template
    rho = pd.concat([coverage.exp, target], axis=1).corr(method='spearman')
    return(1-rho['exp'][1])

def format_performance_matrix(performanceMat):
    performanceDf = pd.DataFrame()
    for vec in performanceMat:
        itDf = pd.DataFrame(vec)
        performanceDf = performanceDf.append(itDf, ignore_index=True)
    performanceDf.columns = [nuc+'.'+str(pos) for nuc,pos in list(it.product("ATCG",range(1,7)))]
    return(performanceDf)

def save_output_json(outdic, filename):
    '''
    Save output of optimization
    '''
    with open(filename, 'w') as f:
        print(json.dumps(outdic), file=f)

def save_de_output(outlist, fileprefix):
    '''
    Save output of optimization
    '''
    with open(fileprefix+'.DE.rho.txt', 'w') as f:
        for score in outlist[0][0]:
            print(score, file=f)
    performanceMat = outlist[0][1]
    performanceDf = format_performance_matrix(performanceMat)
    performanceDf.to_csv(fileprefix+'.DE.matrix.csv')
    return('Output saved!')

def run_DE(deltaGfile, abundanceFile, outfileprefix, popsize=20, its=1000, cores=5):
    seqs = all_hexamers()
    dgMat = pd.read_csv(deltaGfile, index_col=0)
    abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
    genomeAb = abundance[1]
    res = de(coverage_function,(seqs,dgMat,genomeAb,genomeAb), 6, popsize=popsize, its=its, cores=cores)
    outlist = list(res)
    save_de_output(outlist, outfileprefix)
    # outdic = {'score':p[0][0], 'mat':p[0][1].tolist()}
    # save_output_json(outdic,outfile)

deltaGfile = '/hpc/hub_oudenaarden/edann/crypts_bs/VAN2408/CM1_tr2_R1_bismark_bt2_ptDg_qual.csv'
abundanceFile = "/hpc/hub_oudenaarden/edann/hexamers/genomes_kmers/mm10.kmerAbundance.csv"
# run_DE(deltaGfile, abundanceFile, 'match_genomeAb.pop20.it1000.json', cores=10)
run_DE(deltaGfile, abundanceFile, '/hpc/hub_oudenaarden/edann/hexamers/DEotimization/even_cov/test_VAN2408_CM1', popsize=20, its=1000, cores=10)
