# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 18:39:07 2020

@author: eliza

comparing data from mono vs duo
"""
import numpy as np
import copy
import pandas as pd
import os
import glob
import argparse
from datetime import datetime
from dataclasses import dataclass

from setup_file import RunConfig, print_fig, load_data, plot_heatmap

def check_setup_same (mono_results, duo_results):
    """Assert that both simulations that are being compared are compatible"""
    assert mono_results.setup.sim_type == duo_results.setup.sim_type, "simulation types don't match"
    assert mono_results.setup.p == duo_results.setup.p, "parameters don't match"
    assert mono_results.setup.initial_values == duo_results.setup.initial_values, "initial values don't match"
    assert mono_results.setup.setup_start == duo_results.setup.setup_start, "setup start doesn't match"
    assert mono_results.setup.substrate == duo_results.setup.substrate, "substrate doesn't match"
    assert mono_results.setup.solver_args == duo_results.setup.solver_args
    
    if mono_results.setup.sim_type == "parameter_sweep_1D" or mono_results.setup.sim_type == "parameter_sweep_2D":
        assert mono_results.setup.p1_name == duo_results.setup.p1_name, "parameter 1 name doesn't match"
        assert np.all(mono_results.setup.p1_sweep == duo_results.setup.p1_sweep) == True, "parameter 1 sweep values don't match"
        if mono_results.setup.sim_type == "parameter_sweep_2D":
            assert mono_results.setup.p2_name == duo_results.setup.p2_name, "parameter 2 name doesn't match"
            assert np.all(mono_results.setup.p2_sweep == duo_results.setup.p2_sweep) == True, "parameter 2 sweep values don't match"
      
    setup=copy.deepcopy(mono_results.setup)
    setup.now = datetime.now()
    return(setup)

"""Finding the maximum results"""

@dataclass   
class MaxPoints(): 
    point: float
    index: float or tuple
    p1: float
    p2: float

def find_highest_results(results, setup):
    """Search array and identify the maximum result"""
    max_value = np.nanmax(results)
    index_max = np.where(results == max_value)
    p1_index = int(index_max[0])
    p1_max_value = setup.p1_sweep[p1_index]
    if setup.sim_type == "parameter_sweep_2D":
        p2_index = int(index_max[1])
        p2_max_value = setup.p2_sweep[p2_index]
    highest_results = MaxPoints(point=max_value, index=index_max, p1=p1_max_value, p2=p2_max_value)
    return(highest_results) 

def print_highest_results(mono_results, duo_results, setup, value, f):
    """Prints the comparisons for maximum results"""
     
    mono_max = find_highest_results(results=mono_results.results[value], setup=setup)
    f.write(f'  Monoculture max {value}: {mono_max.point} at p1: {mono_max.p1} and p2: {mono_max.p2}\n')
    
    if value in RunConfig().cell_a_only:
        duo_a_max = find_highest_results(results=duo_results.results.cell_a[value] , setup=setup)
        f.write(f'  Cell a max: {value} {duo_a_max.point} at p1: {duo_a_max.p1} and p2: {duo_a_max.p2}\n')
        highest_point_dict = {"mono_max" : mono_max, "duo_a_max" : duo_a_max}
    
    elif value in RunConfig().cell_b_only:
        duo_b_max = find_highest_results(results=duo_results.results.cell_b[value], setup=setup)
        f.write(f'  Cell b max {value}: {duo_b_max.point} at p1: {duo_b_max.p1} and p2: {duo_b_max.p2}\n')
        highest_point_dict = {"mono_max" : mono_max, "duo_b_max" : duo_b_max}
    
    else:
        duo_avg_max = find_highest_results(results=duo_results.results.avg[value], setup=setup)
        f.write(f'  Duoculture max {value}: {duo_avg_max.point} at p1: {duo_avg_max.p1} and p2: {duo_avg_max.p2}\n')

        cell_a_value = duo_results.results.cell_a[value][duo_avg_max.index]
        cell_b_value = duo_results.results.cell_b[value][duo_avg_max.index]
        f.write(f'  Cell a {value} when average is max: {cell_a_value[0]} at p1: {duo_avg_max.p1} and p2: {duo_avg_max.p2}\n')
        f.write(f'  Cell b {value} when average is max: {cell_b_value[0]} at p1: {duo_avg_max.p1} and p2: {duo_avg_max.p2}\n')
        
        duo_a_max = find_highest_results(results=duo_results.results.cell_a[value], setup=setup)
        cell_b_when_a_max = duo_results.results.cell_b[value][duo_a_max.index]
        f.write(f'  Cell a max {value}: {duo_a_max.point} at p1: {duo_a_max.p1} and p2: {duo_a_max.p2}\n')
        f.write(f'  Cell b when a is max: {cell_b_when_a_max}\n')
        
        duo_b_max = find_highest_results(results=duo_results.results.cell_b[value], setup=setup)
        cell_a_when_b_max = duo_results.results.cell_a[value][duo_b_max.index]
        f.write(f'  Cell a max {value}: {duo_b_max.point} at p1: {duo_b_max.p1} and p2: {duo_b_max.p2}')
        f.write(f'  Cell a when b is max: {cell_a_when_b_max}')
        
        highest_point_dict = {"mono_max" : mono_max, "duo_avg_max" : duo_avg_max, "duo_a_max" : duo_a_max, "duo_b_max" : duo_b_max}
    
    return(highest_point_dict)
    
def plot_comparison_graphs(setup, data, title, *args, **kwds):
    """Plots graphs"""
    plot_heatmap(data=data, 
             title= title, 
             xaxis=setup.p2_sweep, 
             xlabel=setup.p2_name, 
             yaxis=setup.p1_sweep, 
             ylabel=setup.p1_name, 
             *args, **kwds)
    figname = title
    print_fig(figname= figname, setup=setup)


"""LOAD DATA"""

#### MAKE SURE YOU ONLY HAVE 1 MONO AND 1 DUO FILE IN THE SAME FOLDER

mono_file = glob.glob('2D_paramsweep_mono*.*')
duo_file = glob.glob('2D_paramsweep_duo*.*')

mono_results=load_data(mono_file[0])
duo_results=load_data(duo_file[0])
setup = check_setup_same(mono_results, duo_results)

####SETUP TXT SAVE

now = setup.now
dt_string = now.strftime("%Y.%m.%d_%H.%M.%S")
folder_name =  "Results-" + dt_string
if not os.path.exists(folder_name):
    os.mkdir(folder_name)
dir_name = os.path.join(os.getcwd(), folder_name)
save_name = "highest results_" + dt_string + ".txt"
file_path = os.path.join(dir_name, save_name)
f = open(file_path, mode='w')


"""COMPARE RESULTS"""

############################### GROWTH RATE ################################################################ 


f.write("Growth rates")
lam_highest_results = print_highest_results(mono_results, duo_results, setup, value="lam", f=f)
f.write("\n")

####DIFFERENCE IN GROWTH RATE

#make copies so any changes don't affect original
lam_mono = copy.deepcopy(mono_results.results["lam"])
lam_duo_cell_a = copy.deepcopy(duo_results.results.cell_a["lam"])
lam_duo_cell_b = copy.deepcopy(duo_results.results.cell_b["lam"])
lam_duo_avg = copy.deepcopy(duo_results.results.avg["lam"])

#######Figure 5a
plot_comparison_graphs(setup=setup,
                       data=lam_mono, 
                       title= "Mono growth rate", 
                       cbarlabel="growth rate mono",
                       cmap="Blues")

#######Figure 5b
plot_comparison_graphs(setup=setup,
                       data=lam_duo_avg, 
                       title= "Duo growth rate", 
                       cbarlabel="growth rate duo avg",
                       cmap="Blues")


diff_lam = lam_duo_avg - lam_mono

#######Figure 5c
plot_comparison_graphs(setup=setup,
                       data=diff_lam, 
                       title= "Duo avg - mono growth",  
                       cbarlabel="growth rate duo - mono",
                       cmap="RdBu",
                       center=0)

#to do further calculations very very low values are removed to avoid extremes in the heatmaps
lam_mono[lam_mono < 1e-10] = np.NaN
lam_duo_avg[lam_duo_avg < 1e-10] = np.NaN
lam_duo_cell_a[lam_duo_cell_a < 1e-10] = np.NaN
lam_duo_cell_b[lam_duo_cell_b < 1e-10] = np.NaN


### DUO - MAX MONO - find parameter values where the duo culture outperforms the best performing monoculture

duo_more_than_mono_lam = lam_duo_avg > lam_highest_results["mono_max"].point

#graph that just shows the parameters where duo outperforms ANY monoculture
plot_comparison_graphs(setup=setup,
                       data=duo_more_than_mono_lam, 
                       title= "Where duo grows better than max mono", 
                       cbarlabel="growth rate duo - max mono",
                       cmap="Blues")

#graph for the actual differences when duo is greater than mono
duo_more_than_mono_lam_2 = np.where(lam_duo_avg > lam_highest_results["mono_max"].point, lam_duo_avg - lam_highest_results["mono_max"].point, np.nan)

#graph for just Difference between duo and max mono
plot_comparison_graphs(setup=setup,
                       data=duo_more_than_mono_lam_2, 
                       title= "lam Duo - mono MAX", 
                       cbarlabel="lam duo - max mono",
                       cmap="Blues")

lam_best = find_highest_results(results=duo_more_than_mono_lam_2, setup=setup)
f.write(f'  Maximum difference between duo and the maximum mono is: {lam_best.point} at p1: {lam_best.p1} and p2: {lam_best.p2}\n')


############################ PROTEINS A and B ###################################################################
f.write("\n")

###HIGHEST RESULTS PER CELL
f.write("A-amylase production PER CELL\n")
ea_highest_results = print_highest_results(mono_results, duo_results, setup, value="ea", f=f)
f.write("Glucoamylase production PER CELL\n")
eb_highest_results = print_highest_results(mono_results, duo_results, setup, value="eb", f=f)
f.write("\n")
# RESULTS ARE GIVEN AS PER CELL

##MULTIPLY BY POPULATIONS

mono_a_protein = (setup.p["Na"] + setup.p["Nb"] ) * mono_results.results["ea"] 
duo_a_protein = setup.p["Na"] * duo_results.results.cell_a["ea"] 
mono_b_protein = ( setup.p["Na"] + setup.p["Nb"] ) * mono_results.results["eb"] 
duo_b_protein = setup.p["Nb"] * duo_results.results.cell_b["eb"]  
both_amylases_mono = mono_a_protein + mono_b_protein
both_amylases_duo = duo_a_protein + duo_b_protein

#remove any values where growth is approx 0 - because of steady state there may be some protein left over
mono_a_protein[mono_results.results["lam"] < 1e-10] = np.NaN
duo_a_protein[duo_results.results.avg["lam"] < 1e-10] = np.NaN
mono_b_protein[mono_results.results["lam"] < 1e-10] = np.NaN
duo_b_protein[duo_results.results.avg["lam"] < 1e-10] = np.NaN
both_amylases_mono[mono_results.results["lam"] < 1e-10] = np.NaN
both_amylases_duo[duo_results.results.avg["lam"] < 1e-10] = np.NaN

diff_a_protein = duo_a_protein - mono_a_protein
diff_b_protein = duo_b_protein - mono_b_protein
diff_both_amylases = both_amylases_duo - both_amylases_mono

plot_comparison_graphs(setup=setup,
                       data=diff_a_protein, 
                       title= "Duo - mono protein a", 
                       cbarlabel="duo - mono",
                       cmap="RdBu",
                       center=0)

plot_comparison_graphs(setup=setup,
                       data=diff_b_protein, 
                       title= "Duo - mono protein b", 
                       cbarlabel="duo - mono",
                       cmap="RdBu",
                       center=0)

plot_comparison_graphs(setup=setup,
                       data=diff_both_amylases, 
                       title= "Duo - mono both amylases", 
                       cbarlabel="duo - mono",
                       cmap="RdBu",
                       center=0)

#########################MALTODEXTRINS & GLUCOSE############################################################

plot_comparison_graphs(setup=setup,
                       data=mono_results.results["s1"], 
                       title= "Mono maltodextrins", 
                       cbarlabel="s1",
                       cmap="Blues")

plot_comparison_graphs(setup=setup,
                       data=duo_results.results.ext["s1"], 
                       title= "Duo maltodextrins", 
                       cbarlabel="s1",
                       cmap="Blues")


plot_comparison_graphs(setup=setup,
                       data=mono_results.results["s"], 
                       title= "Mono glucose", 
                       cbarlabel="s",
                       cmap="Blues")

plot_comparison_graphs(setup=setup,
                       data=duo_results.results.ext["s"], 
                       title= "Duo glucose", 
                       cbarlabel="s",
                       cmap="Blues")

diff_glucose = duo_results.results.ext["s"] - mono_results.results["s"]

plot_comparison_graphs(setup=setup,
                       data=diff_glucose, 
                       title= "Duo-mono glucose", 
                       cbarlabel="duo-mono s",
                       cmap="Blues",
                       center=0)





