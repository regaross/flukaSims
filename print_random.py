import argparse
parser = argparse.ArgumentParser()
parser.add_argument(dest='first_seed', type=int)
args = parser.parse_args()


def get_seed(first_seed):
    '''Changes the seed to the simulation in a given input file'''
    import numpy as np
    import os
    import time

    ### Seed for random number generation from system time
    SEED = first_seed
    np.random.seed(SEED)

    sim_seed = str(np.random.randint(0,999999))

    return sim_seed

print(get_seed(args.first_seed))