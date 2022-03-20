#!/usr/bin/env python3

import random
from deap import base, creator, tools, algorithms
from scipy.spatial.distance import hamming
from tqdm import tqdm
from src.generators import gen_antigen, gen_antibody
from src.germinal_center import GerminalCenter
from src.constants import *

'''
virus HAS antigens
antigens HAVE epitopes
antibodies HAVE proteins

initialilze the virus 
seed the antibodies based on epitopes

once it converges, do a mutation of the virus

each GC gets a single (random) antigen

virus is random long bit string
antigens are portions of virus, sent to GC

start with a virus:
    split it into antigens, send each to GC
        evolve antibodies to fight the antigen
        do this until convergence (with threshold)

evolve the virus:
    seed with antibodies from previous virus
    see how long this takes
'''

#TODO weights correct? what are these for?
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Virus", list)
creator.create("Antigen", list)
creator.create("Antibody", list, fitness=creator.FitnessMin)
creator.create("Epitope", list)

toolbox = base.Toolbox()
toolbox.register("attr_virus", random.randint, 0, 1)
toolbox.register("virus", tools.initRepeat, creator.Virus, toolbox.attr_virus, n=VIRUS_LEN)
virus = toolbox.virus()
toolbox.register("antigen", tools.initIterate, creator.Antigen, gen_antigen(virus))
toolbox.register("antigens", tools.initRepeat, list, toolbox.antigen, n=ANTIGENS_PER_VIRUS)
toolbox.register("antibody", tools.initIterate, creator.Antibody, gen_antibody)
toolbox.register("antibodies", tools.initRepeat, list, toolbox.antibody)
toolbox.register("init_antibody_fitness", lambda x: (0,))

antigens = toolbox.antigens()

gcs = [GerminalCenter(toolbox, antigen, toolbox.antibodies(n=NUM_ANTIBODIES)) for antigen in antigens]

for gc in gcs:
    while True:
        done = gc.evolve()
        if done:
            break
