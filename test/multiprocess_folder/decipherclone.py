import numpy as np
from multiclone import the_multiplicator


def very_complex_calculation():
    list_of_samples = []
    list_of_factors = []
    for x in np.random.randint(10, 100, size=40):
        list_of_samples.append(np.random.randint(0, x, (40, 40)))
        list_of_factors.append(np.random.rand(40, 40))

    y = the_multiplicator(list_of_samples, list_of_factors)
    for sub in y:
        print(sub.mean())
