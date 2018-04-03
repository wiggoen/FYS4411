import matplotlib.pyplot as plt
import numpy as np

min_block_size = 5
max_block_size = 10
n_block_samples = 5
block_ste_length = (max_block_size-min_block_size)/n_block_samples

samples = np.loadtxt("blockres.txt", unpack=True, skiprows=1)
energies = samples[1]

#for i in range(n_block_samples):
#    block_size = min_block_size + i*block_step_length

block_size = 2
block = np.array(n_block_samples)
counter = 0
block_count = 0
for i in range(block_size):
    N = len(samples[1])    
    counter +=1
    print counter 
    if (counter == block_size):
        block[block_count] = energies[(block_count)*block_size : (block_count+1)*block_size]
        block_count += 1
        counter = 0
         
