#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pymds import DistanceMatrix

from antibody import Antibody
from util import hamming_distance, levenshtein_distance

# Assume antibodies can be 6-8 amino acids
def get_antibodies():
    """Returns a list of antibody objects"""
    antibodies = list()
    for i in range(10):
        antibody = Antibody()
        antibody.label = str(i)
        antibodies.append(antibody)
    return antibodies

def get_distance_matrix(antibody_list):
    distance_matrix = list()
    for i in range(len(antibody_list)):
        distance_matrix.append([])
        for j in range(len(antibody_list)):
            distance_matrix[i].append(0)
    for idx1, e1 in enumerate(antibody_list, start=0):
        for idx2, e2 in enumerate(antibody_list, start=0):
            if idx1 == idx2:
                continue
            distance_matrix[idx1][idx2] = levenshtein_distance(e1, e2)
    return np.array(distance_matrix)

def _get_color_from_variant(variant):
    if variant == "Alpha":
        return "red"
    elif variant == "Delta":
        return "blue"
    else:
        return "green"

def create_map(antibody_list):
    dist_matrix = get_distance_matrix(antibody_list)
    print(dist_matrix)
    df = pd.DataFrame(dist_matrix)
    dm = DistanceMatrix(df)
    projection = dm.optimize(n=2)
    coords = projection.coords.iterrows()
    # pymds creates a plot for us, so we're going to get the points, clear it, and
    # then apply our own coloring
    plt.cla()
    for points, antibody in zip(coords, antibody_list):
        x, y = points[1][0], points[1][1]
        plt.scatter(x, y, color=_get_color_from_variant(antibody.evolved_against))
        plt.annotate(antibody.label, xy=(x, y), xytext=(1,1), textcoords="offset points")
    plt.show()

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
