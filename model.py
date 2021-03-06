#!/usr/bin/env python3

import random
from deap import base, creator, tools, algorithms
from scipy.spatial.distance import hamming
from tqdm import tqdm
from src.generators import gen_antigen, gen_antibody
from src.germinal_center import GerminalCenter
from src.constants import *

from antibody_map.antibody import Antibody
from antibody_map.create_map import create_map

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
random.seed(12345)

#Trying to minimize the hamming distance of antibodies from the virus
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

#Modify in place
def mutate(virus):
    num_bits_to_mutate = int(VIRUS_MUTATION_PCT * len(virus))
    bits_to_mutate = random.sample(list(range(len(virus))), num_bits_to_mutate)
    for ind in bits_to_mutate:
        if virus[ind] == 0:
            virus[ind] = 1
        else:
            virus[ind] = 0

gcs = [GerminalCenter(toolbox, antigen, toolbox.antibodies(n=NUM_ANTIBODIES)) for antigen in antigens]
antigen = antigens[0]
gcs = [GerminalCenter(toolbox, antigen, toolbox.antibodies(n=NUM_ANTIBODIES)) for _ in range(ANTIGENS_PER_VIRUS)]

def fight_virus(gc):
    t = 0
    while True:
        t += 1
        done = gc.evolve()
        if done:
            break
    print(f"\t\ttook {t} timesteps to converge")

def create_antibody_map(gcs):
    for gc_id, gc in enumerate(gcs):
        antibodies = list()
        for variant, bit_arrays in gc.best_evolved.items():
            for generation, bit_array in enumerate(bit_arrays):
                antibody = Antibody(
                        bit_string="".join([str(i) for i in bit_array]), 
                        germinal_center_id=gc_id,
                        generation=generation,
                        evolved_against=variant
                    )
                antibodies.append(antibody)
        create_map(antibodies, gc_id)

def main():
    for gc_num, gc in enumerate(gcs):
        print(f"GC {gc_num}:")
        for variant in range(NUM_VIRUS_VARIANTS):
            print(f"\tVirus {variant}:")
            gc.variant = variant
            fight_virus(gc)
            mutate(gc.antigen)
            gc.refresh()
    
    create_antibody_map(gcs)

    # map best from each germinal center
    merged_antibodies = []
    for gc in gcs:
        merged_antibodies.extend(gc.antibodies)

    super_gc = GerminalCenter(toolbox, antigen, merged_antibodies)
    print("SUPER GC:")
    for variant in range(NUM_VIRUS_VARIANTS):
        print(f"\tVirus {variant}:")
        fight_virus(gc)
        mutate(gc.antigen)
        gc.refresh()

if __name__ == "__main__":
    main()
