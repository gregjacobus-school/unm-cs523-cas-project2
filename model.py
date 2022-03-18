#!/usr/bin/env python3

import random
from deap import base, creator, tools, algorithms
from src.amino_acid import AminoAcid
from scipy.spatial.distance import hamming

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
    if len(ep) != len(ant):
        return False
    ham_dist = hamming(ep, ant)
    if ham_dist < 0.5: #TODO is this a good param?
        return True

def remove_none_from_list(x):
    pruned = []
    for val in x:
        if val is not None:
            pruned.append(val)
    return pruned

#Probability of mating (binding) FIXME what's a good value?
cxpb = 0.1

#Run the algorithm
NUM_GEN=40
for _ in range(NUM_GEN):
    #With some probability, have antibodies bind to epitopes
    for ant_ind, ant in enumerate(ant_pop):
        if random.random() < cxpb:
            ep_ind = random.randint(0, len(ep_pop) - 1)
            did_mate = mate(ep_pop[ep_ind], ant)
            if did_mate:
                ep_pop[ep_ind] = None
                ep_pop.remove(None)
print(f"# of antibodies={len(ant_pop)}\n# of epitopes={len(ep_pop)}")
