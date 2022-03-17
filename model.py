#!/usr/bin/env python3

import random
from deap import base, creator, tools, algorithms
from src.amino_acid import AminoAcid

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

toolbox.register("antibody", tools.initRepeat, creator.Antibody)
toolbox.register("antibody_pop", tools.initRepeat, list, toolbox.antibody)
#ant_pop = toolbox.antibody_pop(n=100)

#Run the algorithm
t = 0
while True:
    #Do things
    t += 1
    if t > 10:
        break
