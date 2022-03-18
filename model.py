#!/usr/bin/env python3

import random
from deap import base, creator, tools, algorithms
from src.amino_acid import AminoAcid
from scipy.spatial.distance import hamming
from tqdm import tqdm

#TODO All these parameters -- do these valules make sense?
#Probability of mating (binding) FIXME what's a good value?
CXPB = 0.5
#How similar an epitope/antibody have to be to bind
HAM_CUTOFF = 0.5
#Probability of a cell replicating
REPLICATE_PB = 0.4
#If a cell replicates, probability that the offspring will have a mutation
MUTATE_PB = 0.6
#Number of timesteps
NUM_GEN=40

creator.create("FitnessMin", base.Fitness, weights=(-1.0, -1.0)) #TODO weights correct? what are these for?
creator.create("Epitope", list, fitness=creator.FitnessMin)
creator.create("Antibody", list, fitness=creator.FitnessMin)

def gen_epitope():
    num_amino_acids = random.randint(6,8)
    amino_acids = [AminoAcid().gen_random() for _ in range(num_amino_acids)]
    epitope_bits = []
    for amino_acid in amino_acids:
        epitope_bits.extend(amino_acid.bits)
    return epitope_bits

toolbox = base.Toolbox()
toolbox.register("attr_epitope", gen_epitope)
toolbox.register("epitope", tools.initIterate, creator.Epitope, gen_epitope)
toolbox.register("epitope_pop", tools.initRepeat, list, toolbox.epitope)

ep_pop = toolbox.epitope_pop(n=100)

#For now, antibodies are represented the same way as epitopes
toolbox.register("antibody", tools.initIterate, creator.Antibody, gen_epitope)
toolbox.register("antibody_pop", tools.initRepeat, list, toolbox.antibody)
ant_pop = toolbox.antibody_pop(n=100)


def mate(ep, ant):
    '''Determines if the epitope and antibody should bind'''
    if ep is None:
        print("EP was NONE")
    if ant is None:
        print("ANT was NONE")
    while len(ep) > len(ant):
        ep = ep[:-1]
    while len(ant) > len(ep):
        ant = ant[:-1]
    ham_dist = hamming(ep, ant)
    if ham_dist < HAM_CUTOFF:
        return True

def mutate(cell):
    #Flip a random bit
    ind = random.randint(0, len(cell) - 1)
    if cell[ind] == 0:
        cell[ind] = 1
    else:
        cell[ind] = 0

def should_break(ep_pop, ant_pop):
    if len(ep_pop) == 0:
        print("We are completely cured! All epitopes are gone. Breaking")
        return True
    if len(ant_pop) == 0:
        print("We have died. All of our antibodies are gone and we still have epitopes. Breaking")
        return True
    return False

#Run the algorithm
for _ in tqdm(range(NUM_GEN)):
    #With some probability, have antibodies bind to epitopes (i.e. "mate")
    for ant_ind, ant in enumerate(ant_pop):
        if random.random() < CXPB:
            ep_ind = random.randint(0, len(ep_pop) - 1)
            did_mate = mate(ep_pop[ep_ind], ant)
            if did_mate:
                ep_pop[ep_ind] = None
                ep_pop.remove(None)
                if should_break(ep_pop, ant_pop): break
                ant_pop[ant_ind] = None
                ant_pop.remove(None)

    if should_break(ep_pop, ant_pop): break

    #With some probability, each cell should mutate/replicate
    new_antibodies = []
    for ant in ant_pop:
        if random.random() < REPLICATE_PB:
            clone = ant.copy()
            if random.random() < MUTATE_PB:
                mutate(clone) #modify in place
            new_antibodies.append(clone)
    ant_pop.extend(new_antibodies)

    #Do the same for epitopes
    new_epitopes = []
    for ep in ep_pop:
        if random.random() < REPLICATE_PB:
            clone = ep.copy()
            if random.random() < MUTATE_PB:
                mutate(clone) #modify in place
            new_epitopes.append(clone)
    ep_pop.extend(new_epitopes)

    print(f"# of antibodies={len(ant_pop)}\n# of epitopes={len(ep_pop)}")
print(f"# of antibodies={len(ant_pop)}\n# of epitopes={len(ep_pop)}")
