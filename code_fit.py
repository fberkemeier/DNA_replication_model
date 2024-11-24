from time import monotonic
import cProfile
from typing import Any
import argparse
import sys
import numpy as np
import scipy.optimize
import math
import os
import time

### [FUNCTION] FITTING ALGORITHM ###

def fitfunction(list, v0, st0, fit_step, maxiter, err_threshold, saveQ, info, cell_line='', chr_number=''):
    
    timel = list
    
    v = v0
    st = st0
    exp_v = np.exp(-1/v)
    x00 = np.array([(math.pi/(4*v))*i**(-2) for i in timel])
    lm = 1000 # Remove end regions for error calculation
    
    # VECTORIZED APPROACH
    
    def mse(y_true, y_pred):
        mse_value = sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)) / len(y_true)
        return mse_value
    
    def fast_roll_add(dst, src, shift):
        dst[shift:] += src[:-shift]
        dst[:shift] += src[-shift:]
    
    # Expected replication time computation (replaces bcs)
    def fp(x, L, v):
        n = len(x)
        y = np.zeros(n)
    
        last_exp_2_raw = np.zeros(n)
        last_exp_2 = np.ones(n)
        unitary = x.copy()
        for k in range(L+1):
            if k != 0:
                fast_roll_add(unitary, x, k)
                fast_roll_add(unitary, x, -k)
            exp_1_raw = last_exp_2_raw
            exp_1 = last_exp_2
            exp_2_raw = exp_1_raw + unitary / v
            exp_2 = np.exp(-exp_2_raw)
    
            # Compute the weighted sum for each j and add to the total
            y += (exp_1 - exp_2) / unitary
            
            last_exp_2_raw = exp_2_raw
            last_exp_2 = exp_2
        return y

    # Fitting iteration
    def fitf(time, lst, x0, j, fit_step):
        return x0[j] * (lst[j] / time[j])**(fit_step)

    # Alternative fitting
    def fitf0(time, lst, x0, j, fit_step):
        return x0[j]**(np.log(time[j]) / np.log(lst[j]))

    # Fitting control
    def cfit(time, lst, x0, fit_step):
        result = np.empty_like(x0)
        for j in range(len(x0)):
            fit_result = fitf(time, lst, x0, j, fit_step)
            if fit_result < 10**(-err_threshold):
                result[j] = 10**(-err_threshold)
            else:
                result[j] = fit_result
        return result
    
    xs = x00
    ys = fp(xs, len(xs)//st, v)
    new_err0 = mse(timel[lm:-lm], ys[lm:-lm])
    err = 10**10

    # Open the file to store the error values
    with open(f'data/whole-genome_mse/mse_{cell_line}_chr[{chr_number}].txt', 'a') as mse_file:
        # Write the initial error to the file before the loop
        mse_file.write(f'{new_err0:.30f}\n')

        for j in range(maxiter):
            xs0 = xs
            ys0 = ys
            xs = cfit(timel, ys, xs, fit_step)
            ys = fp(xs, len(xs)//st, v)
            
            new_err = mse(timel[lm:-lm], ys[lm:-lm])
            print(str(j+1) + '/' + str(maxiter) + ' err: ' + str('{:.30f}'.format(new_err)), end="\r")
            
            # Write the new error to the file
            mse_file.write(f'{new_err:.30f}\n')
            
            err = new_err  # Update the error with the new calculated error

    fire_rates = ['{:.30f}'.format(i) for i in xs]
    time_sim = ys
    
    if saveQ:
        with open(r'data/whole-genome_firing_rates/fire_rates_'+info+'.txt', 'w') as f:
            for rate in fire_rates:
                f.write(rate + '\n')
        np.savetxt(r'data/whole-genome_timing_simulation/time_sim_'+info+'.txt', time_sim, fmt='%.30f')
    
    return [fire_rates, time_sim]




### [FUNCTION] DATA GENERATION ###

def datagenfs(cell_line, chr_number, chrpos_min, chrpos_max, resolution, alld, dtscale, saveQ, info, sigscale=0):
    if alld:
        time_data = np.loadtxt(f'data/whole-genome_timing_data/time_data_{cell_line}_chr[{chr_number}].txt', dtype=float)
    else:
        time_data = np.loadtxt(f'data/whole-genome_timing_data/time_data_{cell_line}_chr[{chr_number}].txt', dtype=float)[chrpos_min:chrpos_max]
        np.savetxt(f"data/whole-genome_timing_data/time_data_{cell_line}_chr[{chr_number}]_{chrpos_min}-{chrpos_max}.txt", time_data, fmt='%.30f')
    return time_data



### EXAMPLE: FITTING ###

# Model parameters
cell_line = "HUVEC"
chr_number = 1
hpcQ = True # Option to run in HPC for whole-genome results
if hpcQ:
    parser = argparse.ArgumentParser()
    parser.add_argument("-cl", required=False)
    parser.add_argument("-cn", required=False)
    args = parser.parse_args()
    if len(sys.argv)>1 :
        if '-cl' in sys.argv:
            cell_line = str(args.cl)
        if '-cn' in sys.argv:
            chr_number = int(args.cn)

chrpos_min = 10000
chrpos_max = 20000
x = np.linspace(chrpos_min, chrpos_max, chrpos_max - chrpos_min)  # Chromosome positions
fork_speed = 1.4 # Fork speed
resolution = 1000 # (1 kb)
scale_factor = 6 # Scales the data
all_dataQ = True # Picks whether to fit an entire genome


# Fitting parameters
int_width = 2000
def int_widthf(time_data): return int(len(time_data)/int_width)
fit_step = 2
iterations = 100
err_threshold = 15

# Saving (Warning: replaces existing files)
saveQ = False
file_name = f'{cell_line}_chr[{chr_number}]' if all_dataQ else f'{cell_line}_chr[{chr_number}]_{chrpos_min}-{chrpos_max}'

# Single files
sing_filesQ = True
if sing_filesQ:
    # Data generation
    time_data = datagenfs(cell_line, chr_number, chrpos_min, chrpos_max, resolution, all_dataQ, scale_factor, saveQ, file_name)
    # Fitting
    fire_rates, time_sim = fitfunction(time_data, fork_speed, int_widthf(time_data), fit_step, iterations, err_threshold, saveQ, file_name,
                                      cell_line=cell_line, chr_number=chr_number)

# Multiple file fitting (long computation)
mult_fileQ = False
if mult_fileQ:
    # Whole-genome parameters
    cell_lines = ["HeLa-S3","BJ1","IMR90","HUVEC","K562","GM12878","HepG2","MCF-7"]
    chr_range = range(1,23)
    for cell_line_i in cell_lines:
        for chr_number_i in chr_range:
            print(cell_line_i+' chr '+str(chr_number_i)+'/22')
            file_name = cell_line_i+'_chr['+str(chr_number_i)+']' if all_dataQ else cell_line_i+'_chr['+str(chr_number_i)+']_'+str(chrpos_min)+'-'+str(chrpos_max)
            # Data generation
            time_data = datagenfs(cell_line_i, chr_number_i, chrpos_min, chrpos_max, resolution, all_dataQ, scale_factor, saveQ, file_name)
            # Fitting
            fire_rates, time_sim = fitfunction(time_data, fork_speed, int_widthf(time_data), fit_step, iterations, err_threshold, saveQ, file_name)
            