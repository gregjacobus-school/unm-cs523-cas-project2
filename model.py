#!/usr/bin/env python3

from deap import base, creator, tools, algorithms

EPITOPE_SIZE = 8

creator.create("FitnessMin", base.Fitness, weights=(-1.0, -1.0)) #TODO weights correct?
creator.create("Epitope", list, fitness=creator.FitnessMin)
creator.create("Antibody", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
toolbox.register("epitope", tools.initRepeat, creator.Epitope, n=EPITOPE_SIZE)
toolbox.register("antibody", tools.initRepeat, creator.Antibody, n=EPITOPE_SIZE)
