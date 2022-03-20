import matplotlib.pyplot as plt

import pandas as pd
from pymds import DistanceMatrix

# Distances between the vertices of a right-angled triangle
dist = pd.DataFrame({
    'a': [0, 1, 1, 1],
    'b': [1, 0, 1, 1],
    'c': [1, 1, 0, 1],
    'd': [1, 1, 1, 0]},
    index=['a', 'b', 'c', 'd'])

# Make an instance of DistanceMatrix
dm = DistanceMatrix(dist)

# Embed vertices in two dimensions
projection = dm.optimize(n=2)

# Plot projection
projection.plot(c='black', edgecolor='white', s=50)

plt.show()
