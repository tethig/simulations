#!/usr/bin/env python3

__author__ = "Ben Dickins"
__email__ = "benjamin.dickins@ntu.ac.uk"
__status__ = "Draft"
__version__ = "2.3.0"

import sys, copy
import numpy as np
import scipy.stats as stats
from copy import deepcopy
# from numba import jit
# not good for methods - and really speeds up native code so need to make generic outer function for mutate and jit decorate
# http://stackoverflow.com/questions/25683592/whenever-i-try-to-use-jit-on-my-class-method-i-get-indentationerror-unexpecte

# define distribution of mutational effect sizes as truncated normal (Vale et al)
lower, upper, mu, sigma = 0, 2, 0.71, 0.31
dfem = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)

# define the viral class (no methods allowed for parallel code)
class Virus:
    def __init__(self, size, genome=None, phenotype=None, mutations=None, fitness=None):
        self.size = size
        if genome is None:
            self.genome = np.zeros(shape=size, dtype=int)
        if phenotype is None:
            self.phenotype = np.ones(shape=size, dtype=float)
        if mutations is None:
            self.mutations = np.sum(self.genome)
        if fitness is None:
            self.fitness = np.prod(self.phenotype)

    def __str__(self): # for debugging
        return "Mutations: {}; Fitness: {}".format(self.mutations, self.fitness)

    def reproduce(self, fecundity): # this method is still too slow
        familysize = int(round((fecundity * self.fitness)))
        sprogs = []
        for v in range(familysize):
            kid = Virus(L)
            kid.genome = np.copy(self.genome)
            kid.phenotype = np.copy(self.phenotype)
            kid.mutations = np.copy(self.mutations)
            kid.fitness = np.copy(self.fitness)
            sprogs.append(kid)
        return sprogs

#    def reproduce(self, fecundity): # this method is very slow
#        familysize = int(round((fecundity * self.fitness)))
#        return [deepcopy(self) for v in range(familysize)]

# for debugging
def print_attributes(obj):
    for attr in obj.__dict__:
        print(attr, getattr(obj, attr))

# mutation function (relatively fast)
# not convinced of a speedup from @jit(nogil=True)
def mutate(virion, rate):
    newmuts = min(virion.size, int(np.random.poisson(rate, 1))) # truncate number of mutations to genome length
    positions = np.random.choice(range(virion.size), size=newmuts, replace=False)
    virion.genome[positions] = abs(virion.genome[positions] + 1 - 2) # change 0 to 1 and 1 to 0
    virion.phenotype[positions] = ((dfem.rvs(newmuts) - 1) * virion.genome[positions]) + 1 # if reversion (to 0) then fitness=1, else from dfem
    virion.mutations = np.sum(virion.genome)
    virion.fitness = np.prod(virion.phenotype)

# recombination function (executes quickly)
def recombine(dam, sire):
    breakpoint = int(np.random.choice(range(L), 1))
    offspring = Virus(L)
    offspring.genome = np.concatenate([ dam.genome[0:breakpoint], sire.genome[breakpoint:L] ])
    offspring.phenotype = np.concatenate([ dam.phenotype[0:breakpoint], sire.phenotype[breakpoint:L] ])
    offspring.mutations = np.sum(offspring.genome)
    offspring.fitness = np.prod(offspring.phenotype)
    return offspring

# main program loop
if __name__ == '__main__':

    # fetch user inputs
    B = int(input("What is the starting population size? [B=100] ") or 100)
    L = int(input("What is the genome length? [L=20] ") or 20)
    U = float(input("What is the genomic mutation rate (rem. capped to L)? [U=1] ") or 1)
    X = float(input("What is the recombination rate? [X=0.01] ") or 0.01)
    R = int(input("What is the average burst size? [R=10] ") or 10)
    S = int(input("Maximum # of generations for simulation? [S=30] ") or 30)

    # open connections and print headers
    maffile = open("mafs.txt", "w")
    mutfile = open("muts.txt", "w")
    print("StartPop{} GenomLen{} GenMutRate{} RecombRate{} BurstSize{} SimLimit{}".format(B,L,U,X,R,S), file=maffile)
    print("Gen\tPop", "\t".join([str(x) for x in range(1,L+1)]), sep="\t", file=maffile)
    print("Gen\tMut\tFit", file=mutfile)

    # main simulation loop
    population = [Virus(L) for x in range(B)]
    for gen in range(S):
        nextgen, mutnum, fitnum = [], [], []
        for organism in population:
            children = organism.reproduce(R)
            for child in children: # convert to list comprehension
                mutate(child, U)
            nextgen.extend(children)
            mutnum.extend(child.mutations for child in children)
            fitnum.extend(child.fitness for child in children)

        # print data to file
        for mut, fit in zip(mutnum,fitnum):
            print(*(gen+1,mut,fit), sep="\t", file=mutfile)

        # report lethal mutagenesis
        popsize = len(nextgen)
        if popsize == 0:
            print("Lethal mutagenesis achieved...")
            break

        # recombination
        if X > 0:
            recnum = round(X * popsize)
            recnum = recnum + ( recnum % 2 ) # rounding up to even number
            sexindex = np.random.choice(range(popsize), size=recnum, replace=False)
            parents = [nextgen[i] for i in sexindex]
            # pick parents and make children
            for i in range(0, recnum, 2): # be careful here
                mum, dad = parents[i], parents[i+1]
                lovechild = recombine(mum, dad)
                nextgen.append(lovechild)
            # remove parents
            sexrevsd = np.sort(sexindex)[::-1] # reverse ensures accurate traversal
            for r in sexrevsd: # because pop removes elements
                nextgen.pop(r)
                mutnum.pop(r)
                fitnum.pop(r)

        # recycle population
        population = nextgen # next is now current also

        # collect statistics for standard out
        print("\nPopulation Size:",popsize)
        print("# of Mutation-Free Viruses:",sum([x==0 for x in mutnum]))
        print("Median # of Mutations:",np.median(mutnum))
        print("Mean # of Mutations:",np.mean(mutnum))
        print("Max # of Mutations:",max(mutnum))
        print("Median Fitness:",np.median([x.fitness for x in nextgen]))
        print("Mean Fitness:",np.mean([x.fitness for x in nextgen]))
        print("Maximum Fitness:",max([x.fitness for x in nextgen]))
        print(gen+1, popsize, sep='\t', end='', file=maffile)
        for base in range(L):
            hits = sum([x.genome[base] for x in nextgen])
            print("\t{}".format(str(hits/popsize)), end='', file=maffile)
        print("", file=maffile)

    maffile.close()
    mutfile.close()
    print("End of simulation")
