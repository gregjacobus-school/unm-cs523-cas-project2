#! /usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import random
import os
from src.constants import *

from pymds import DistanceMatrix

try:
    from antibody_map.antibody import Antibody
except:
    from antibody import Antibody

try:
    from antibody_map.util import hamming_distance, levenshtein_distance
except:
    from util import levenshtein_distance

# Assume antibodies can be 6-8 amino acids
def get_antibodies():
    """Returns a list of antibody objects"""
    antibodies = list()
    for i in range(10):
        antibody = Antibody()
        antibody.generation = str(i)
        antibodies.append(antibody)
    return antibodies

def get_distance_matrix(antibody_list):
    print(f"Generating map for best {len(antibody_list)} evolved antibodies")
    size = len(antibody_list)
    distance_matrix = np.zeros((size, size))
    for idx1, e1 in enumerate(antibody_list, start=0):
        for idx2, e2 in enumerate(antibody_list, start=0):
            if idx1 <= idx2:
                continue # just solve the bottom half then copy since we're symmetric
            distance_matrix[idx1][idx2] = levenshtein_distance(e1, e2)
    return np.triu(distance_matrix.T,1) + distance_matrix

def _get_color_from_variant(variant):
    if type(variant) == int:
        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
                "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]
        return colors[variant]
    else:
        if variant == "Alpha":
            return "red"
        elif variant == "Delta":
            return "blue"
        else:
            return "green"

def create_map(antibody_list, gc_num):
    os.makedirs('./output', exist_ok=True)
    dist_matrix = get_distance_matrix(antibody_list)
    df = pd.DataFrame(dist_matrix)
    dm = DistanceMatrix(df)
    print("Solving dimensionality reduction...")
    projection = dm.optimize(n=2)
    coords = projection.coords.iterrows()
    # pymds creates a plot for us, so we're going to get the points, clear it, and
    # then apply our own coloring
    plt.cla()
    for points, antibody in zip(coords, antibody_list):
        x, y = points[1][0], points[1][1]
        variant = antibody.evolved_against
        plt.scatter(x, y, color=_get_color_from_variant(variant), s=10)
        plt.grid(True)
        '''
        if random.uniform(0, 1) > 0.8:
            plt.annotate(f"{antibody.germinal_center_id}:{antibody.generation}", 
                    xy=(x, y), xytext=(1,1), textcoords="offset points")
        '''
    #plt.show()
    handles = []
    for variant in range(NUM_VIRUS_VARIANTS):
        color = _get_color_from_variant(variant)
        patch = mpatches.Patch(color=color, label=f'Variant {variant}')
        handles.append(patch)
    plt.legend(handles=handles)
    plt.savefig(f'output/antibody_map-gc-{gc_num}.eps')

def main():
    ab = antibodies = get_antibodies()
    create_map(ab)

def _test_distance_metric():
    for i in range(0, len(antibodies), 2):
        a1 = antibodies[i]
        a2 = antibodies[i+1]
        bit_dist = levenshtein_distance(a1.bit_string, a2.bit_string)
        base_dist = levenshtein_distance(a1.base_string, a2.base_string)
        if base_dist > bit_dist:
            print(bit_dist, base_dist)
            print(a1)
            print(a2)

if __name__ == "__main__":
    main()
