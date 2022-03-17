#!/usr/bin/env python3

import random
from deap import base, creator, tools, algorithms
from src.amino_acid import AminoAcid

EPITOPE_SIZE = 8

creator.create("FitnessMin", base.Fitness, weights=(-1.0, -1.0)) #TODO weights correct?
creator.create("Epitope", list, fitness=creator.FitnessMin)
creator.create("Antibody", list, fitness=creator.FitnessMin)

def gen_epitope():
    num_amino_acids = random.randint(6,8)
    amino_acids = [AminoAcid for _ in range(num_amino_acids)]
    return amino_acids

toolbox = base.Toolbox()
toolbox.register("attr_epitope", gen_epitope)
toolbox.register("epitope", tools.initIterate, creator.Epitope, gen_epitope)
#toolbox.register("antibody", tools.initRepeat, creator.Antibody)

