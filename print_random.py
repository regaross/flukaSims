def get_seed():
    '''Changes the seed to the simulation in a given input file'''
    import numpy as np
    import os
    import time

    ### Seed for random number generation from system time
    SEED = int(time.time())
    np.random.seed(SEED)

    sim_seed = str(np.random.randint(0,999999))

    return sim_seed

print(get_seed())