# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 11:41:32 2024

@author: aaron.cone
"""

#%%  

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import os
import seaborn as sns 

# Directory containing the files
directory = 'D:/Fiber Photometry Food Restrict'

# Dictionaries to store DataFrames
dfs_iso = {}
dfs_gcamp = {}

# Function to generate the dictionary key from the filename
def generate_key(filename):
    parts = filename.split('_')
    key = '_'.join(parts[2:6])  # Adjust this based on the filename structure
    return key

# Loop through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith('.CSV'):
        filepath = os.path.join(directory, filename)
        key = generate_key(filename)
        if '415_Signal' in filename:
            dfs_iso[key] = pd.read_csv(filepath)
        elif '470_Signal' in filename:
            dfs_gcamp[key] = pd.read_csv(filepath)

# Ensure that all keys match in both dictionaries
matching_keys = set(dfs_iso.keys()).intersection(set(dfs_gcamp.keys()))
non_matching_keys = set(dfs_iso.keys()).symmetric_difference(set(dfs_gcamp.keys()))

# Make DataFrames for matching keys the same length
for key in matching_keys:
    len_iso = len(dfs_iso[key])
    len_gcamp = len(dfs_gcamp[key])
    min_length = min(len_iso, len_gcamp)
    dfs_iso[key] = dfs_iso[key].iloc[:min_length]
    dfs_gcamp[key] = dfs_gcamp[key].iloc[:min_length]


#%%

# ZScore the 415 and 470 Signals and adds time to the dataframe

for mouse in dfs_iso:
    dfs_iso[mouse]['Region0G_zscore'] = stats.zscore(dfs_iso[mouse]['Region0G'])
    dfs_iso[mouse]['Time'] = pd.Series(np.arange(start = 0, stop=len(dfs_iso[mouse]['Region0G_zscore']), step=0.0331))


for mouse in dfs_gcamp:
    dfs_gcamp[mouse]['Region0G_zscore'] = stats.zscore(dfs_gcamp[mouse]['Region0G'])
    dfs_gcamp[mouse]['Time'] = pd.Series(np.arange(start = 0, stop=len(dfs_gcamp[mouse]['Region0G_zscore']), step=0.0331))



#%%
from itertools import chain


### FOR PERI-EVENT ANALYSIS ###

#Peri-event histogram for continuous values.
def contvar_peh(var_ts, var_vals, ref_ts, min_max, bin_width = False):
    
    if bin_width:
        ds_ts = np.linspace(var_ts.min(), var_ts.max(), int((var_ts.max()-var_ts.min())/bin_width))
        ds_vals = np.interp(ds_ts, var_ts, var_vals)
        rate = bin_width
    
    else:
        rate = np.diff(var_ts).mean()
        ds_ts, ds_vals = (np.array(var_ts), np.array(var_vals))       
        
    left_idx = int(min_max[0]/rate)
    right_idx = int(min_max[1]/rate)
    
    all_idx = np.searchsorted(ds_ts,ref_ts, "right")   
    all_trials = np.vstack([ds_vals[idx+left_idx:idx+right_idx] for idx in all_idx])
    
    return all_trials



######### GCAMP SIGNAL ##########



# Process all animals and store results in food_depr_470
food_depr_470 = {}

to_start = -150
to_end = 150


for mouse in dfs_gcamp:
        events = contvar_peh(var_ts = dfs_gcamp[mouse]['Time'], 
                              var_vals = dfs_gcamp[mouse]['Region0G_zscore'],
                              ref_ts = [300], 
                              min_max = (to_start,to_end), 
                              bin_width = False)
        # break
        food_depr_470[mouse]=events ## saves array for each peri-event for all animals




# creates time series for plotting, add a zero after events (e.g. num = events[0}.size])
time_peri_event = np.linspace(start = to_start, stop = to_end, num = events[0].size, retstep=0.016)

# Formatted so can process SEM and heatmap
to_line_points = chain.from_iterable(food_depr_470.values())#### TImestamped ambush events

lined_up_points = np.array(list(to_line_points))

# Formatted so can process SEM and heatmap
lined_up_points = np.mean(list(food_depr_470.values()), axis=1)


#%%

# Calculates means for all data points
points = lined_up_points.mean(axis=0)

# Calculates standard error of mean for data points
points_sem = stats.sem(lined_up_points)

# Creates dataframe to plot 
to_plot = pd.DataFrame({'Time': time_peri_event[0], 'zdFF': points})
#%%
# Make figure
fig, ax = plt.subplots(figsize=(16, 10)) # you change dimensions of plot here Width x Length

# ax = plt.gca() # needed for line below - change y axis min and max
# ax.set_ylim([-1.5, 2.5]) #change y axis min and max

# Makes line plot
ax.plot('Time', 'zdFF', data = to_plot)
ax.fill_between(to_plot['Time'], to_plot['zdFF'] - points_sem, 
                   to_plot['zdFF'] + points_sem, alpha=0.15)
ax.set_xlabel('Time (sec)')
ax.set_ylabel('Z-Score')
ax.set_title('Aggregrated 470 Signal')
ax.margins(x=0)        


fig, ax = plt.subplots(figsize=(16, 10)) # you change dimensions of plot here Width x Length

# Heatmap
sns.heatmap(lined_up_points, cbar = True, xticklabels = False, yticklabels = False)
ax.set_title('Aggregrated 470 Signal')
#%%

######### ISOBESTIC SIGNAL ##########


food_depr_415 = {}


to_start = -150
to_end = 150


for mouse in dfs_gcamp:
        events = contvar_peh(var_ts = dfs_iso[mouse]['Time'], 
                             var_vals = dfs_iso[mouse]['Region0G_zscore'],
                             ref_ts = [300], 
                             min_max = (to_start,to_end), 
                             bin_width = False)
        # break
        food_depr_415[mouse]=events ## saves array for each peri-event for all animals





# creates time series for plotting, add a zero after events (e.g. num = events[0}.size])
time_peri_event = np.linspace(start = to_start, stop = to_end, num = events[0].size, retstep=0.016)

# Formatted so can process SEM and heatmap
to_line_points = chain.from_iterable(food_depr_415.values())#### TImestamped ambush events

lined_up_points = np.array(list(to_line_points))

# Formatted so can process SEM and heatmap
lined_up_points = np.mean(list(food_depr_415.values()), axis=1)


#%%

# Calculates means for all data points
points = lined_up_points.mean(axis=0)

# Calculates standard error of mean for data points
points_sem = stats.sem(lined_up_points)

# Creates dataframe to plot 
to_plot = pd.DataFrame({'Time': time_peri_event[0], 'zdFF': points})
#%%
# Make figure
fig, ax = plt.subplots(figsize=(16, 10)) # you change dimensions of plot here Width x Length

# ax = plt.gca() # needed for line below - change y axis min and max
# ax.set_ylim([-1.5, 2.5]) #change y axis min and max

# Makes line plot
ax.plot('Time', 'zdFF', data = to_plot)
ax.fill_between(to_plot['Time'], to_plot['zdFF'] - points_sem, 
                   to_plot['zdFF'] + points_sem, alpha=0.15)
ax.set_xlabel('Time (sec)')
ax.set_ylabel('Z-Score')
ax.set_title('Aggregrated 415 Signal')
ax.margins(x=0)       

fig, ax = plt.subplots(figsize=(16, 10)) # you change dimensions of plot here Width x Length

# Plot heatmap
sns.heatmap(lined_up_points, cbar = True, xticklabels = False, yticklabels = False)
ax.set_title('Aggregrated 415 Signal')
#%%



