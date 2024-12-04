from multiprocessing import Pool
import numpy as np
from tqdm import tqdm


def some_func(lst1, lst2):
    return lst1, lst2


if __name__ == '__main__':
    with Pool(processes=4) as pool:
        def the_multiplicator(l1, l2):
            none_list_results = list(tqdm(map(some_func, l1, l2), total=len(l1)))
