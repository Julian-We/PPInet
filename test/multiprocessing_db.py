from multiprocessing import Pool
import numpy as np
from tqdm import tqdm
from aflib.decipher import helper_generate


if __name__ == '__main__':
    with Pool(processes=4) as pool:
        none_list_results = list(tqdm(map(helper_generate, list_exp_path), total=len(list_exp_path)))
