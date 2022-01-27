# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 10:18:04 2020
@author: Eliza
Scans two parameters
Enabled for multiple processors if available
"""

#setup
import time
import os
import pandas as pd
from functools import partial
from dataclasses import dataclass
from dask.distributed import Client, LocalCluster

from setup_file import arg_parsing, setup_default, process_2D_sweep, check_server, save_data, updateProgressBar, full_title_dict, results_to_plot, plot_heatmap, print_fig
from run_ODE import run_sim_conditional

@dataclass
class SimResultsP1P2:
    """Custom dataclass for 2D sweep sim results"""
    p1: int #value of the parameter being changed
    p2: int #value of the 2nd parameter being changed
    solution: tuple #simulation results
    
@dataclass
class SweepResults2D:
    """Dataclass containing all the setup and solutions"""
    def __init__(self, data, setup):
        self.setup= setup
        self.data= data
    
    def process_data(self):
        self.results = process_2D_sweep(data=self.data, setup=self.setup)
        
    def save_data(self):
        dataname = "2D_paramsweep_" + self.setup.ODE_type + "_" + self.setup.p1_name
        save_data(data=self, dataname=dataname, setup=self.setup)


"""Basic sweep def and distribution"""

def param_scan_2D(setup, comb_param_values):
    """Run sim with 2 specific parameter values - accepts all the inputs of the run_ODE file"""
    selected_p1 = setup.p1_name
    selected_p2 = setup.p2_name
    if not selected_p1 in dict.keys(setup.p):
        raise ValueError("Parameter 1 does not exist") 
    elif not selected_p2 in dict.keys(setup.p):
        raise ValueError("Parameter 2 does not exist") 
    else: 
        #simulate
        setup.p[selected_p1] = comb_param_values[0]
        setup.p[selected_p2] = comb_param_values[1]
        sim_results= run_sim_conditional(default_setup=setup)
        return SimResultsP1P2(p1=comb_param_values[0], p2=comb_param_values[1], solution=sim_results.solution)

def run_sweep_client(claster, setup):
    """Send sims to different processors"""
    # computation
    client = Client(claster)
    print(client)
    # partial
    f = partial(param_scan_2D, setup)
    A = client.map(f, setup.comb_param_values)
    data = client.gather(A)
    return(data)

"""Output files"""

def plot_2D_results_v_params(results, setup):
    """Plotting heatmaps vs parameter values"""
    titles = full_title_dict(setup)
    molecs_to_plot= ["ea", "eb", "s0", "s1", "s", "lam"]
    for i in molecs_to_plot:
        data = results_to_plot(i=i, setup=setup, results=results)
        if len(data) == 3:
            names = ["cell a","cell b","avg"]
            index = 0
        for datum in data:
            if len(data) == 3: 
                title = setup.ODE_type.capitalize() + " " +titles[i] + " " + names[index]
                index = index+1
            else:
                title = setup.ODE_type.capitalize() + " " +titles[i]
            plot_heatmap(data=datum, 
                         title= title, 
                         xaxis=setup.p2_sweep, 
                         xlabel=setup.p2_name, 
                         yaxis=setup.p1_sweep, 
                         ylabel=setup.p1_name, 
                         cmap = "coolwarm",
                         cbarlabel="amount of " + i)
            #save figure
            figname = title + " vs " + setup.p1_name + " and " + setup.p2_name
            print_fig(figname= figname, setup=setup) 


def param_2D_sweep_to_excel(sweep_2D_package):    
    
    now = sweep_2D_package.setup.now
    dt_string = now.strftime("%Y.%m.%d_%H.%M.%S")
    
    folder_name =  "Results-" + dt_string
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    dir_name = os.path.join(os.getcwd(), folder_name)
    excel_single_sim_name = "param_sweep_2D_" + sweep_2D_package.setup.ODE_type + "_" + dt_string + ".xlsx"
    file_path = os.path.join(dir_name, excel_single_sim_name)
    
    parameters_df = pd.DataFrame.from_dict(sweep_2D_package.setup.p, orient='index')
    initial_values_df = pd.DataFrame.from_dict(sweep_2D_package.setup.initial_values, orient='index')
    other_info = {'ODE_type': sweep_2D_package.setup.ODE_type, 
                  'substrate': sweep_2D_package.setup.substrate, 
                  'steady_state_start': sweep_2D_package.setup.setup_start,
                  'p1':sweep_2D_package.setup.p1_name,
                  'p1_sweep': sweep_2D_package.setup.p1_sweep,
                  'p2': sweep_2D_package.setup.p2_name,
                  'p2_sweep': sweep_2D_package.setup.p2_sweep}
    other_info_df = pd.DataFrame.from_dict(other_info, orient='index')
    
    writer = pd.ExcelWriter(file_path, engine='xlsxwriter')
    parameters_df.to_excel(writer, sheet_name='parameters')
    initial_values_df.to_excel(writer, sheet_name='initial states')
    other_info_df.to_excel(writer, sheet_name='input arguments')
    
    if sweep_2D_package.setup.ODE_type == "mono":
        for i in sweep_2D_package.results:
            index=[]
            columns=[]
            for p1_value in sweep_2D_package.setup.p1_sweep:
                index.append("p1= " + str(p1_value))
            for p2_value in sweep_2D_package.setup.p2_sweep:
                columns.append("p2= " + str(p2_value))         
                
            results_df = pd.DataFrame(sweep_2D_package.results[i], index=index, columns=columns)
            results_df.to_excel(writer, sheet_name=i)
    else:
        for i in sweep_2D_package.results.cell_a:
            index=[]
            columns=[]
            for p1_value in sweep_2D_package.setup.p1_sweep:
                index.append("p1= " + str(p1_value))
            for p2_value in sweep_2D_package.setup.p2_sweep:
                columns.append("p2= " + str(p2_value))         
                
            results_df = pd.DataFrame(sweep_2D_package.results.cell_a[i], index=index, columns=columns)
            results_df.to_excel(writer, sheet_name="cell a " + i)
            
        for i in sweep_2D_package.results.cell_b:
            index=[]
            columns=[]
            for p1_value in sweep_2D_package.setup.p1_sweep:
                index.append("p1= " + str(p1_value))
            for p2_value in sweep_2D_package.setup.p2_sweep:
                columns.append("p2= " + str(p2_value))         
                
            results_df = pd.DataFrame(sweep_2D_package.results.cell_b[i], index=index, columns=columns)
            results_df.to_excel(writer, sheet_name="cell b " + i)  
            
        for i in sweep_2D_package.results.avg:
            index=[]
            columns=[]
            for p1_value in sweep_2D_package.setup.p1_sweep:
                index.append("p1= " + str(p1_value))
            for p2_value in sweep_2D_package.setup.p2_sweep:
                columns.append("p2= " + str(p2_value))         
                
            results_df = pd.DataFrame(sweep_2D_package.results.avg[i], index=index, columns=columns)
            results_df.to_excel(writer, sheet_name="average " + i) 
   
    writer.save()

"""Run sweep"""
    
def parameter_sweep_2D(default_setup):
    """Changes the value of 2 parameters over a range of values, runs a simulation
    for each pair of values and returns the solutions for each"""
    
    #step 2: run sweep based on enviornment
    t = time.time()
    #on server run through the run_sweep_client function
    if check_server() == True:
            local = LocalCluster(n_workers=32)
            data=run_sweep_client(claster=local, setup=default_setup)
    #off server run with a for loop
    else:
        data = []
        i=0
        complete = len(default_setup.comb_param_values)
        #for each combination of parameters run a sweep
        for both_ps in default_setup.comb_param_values:
            simulation = param_scan_2D(setup=default_setup, comb_param_values=both_ps)
            data.append(simulation)
            #progress bar
            i = i + 1
            percent = int(i / complete * 100)
            updateProgressBar(percent)
            time.sleep(0.05)
    #print time
    elapsed = time.time() - t
    print(f'Parameter sweep took {elapsed} seconds')
    
    #step 3: package data
    sweep_2D_package = SweepResults2D(setup=default_setup, data=data)
    
    #step 4: process data
    sweep_2D_package.process_data()
    
    #save data
    sweep_2D_package.save_data()    
    param_2D_sweep_to_excel(sweep_2D_package)
    
    #plot data
    plot_2D_results_v_params(results=sweep_2D_package.results, setup=sweep_2D_package.setup)
    
    return(sweep_2D_package)

if __name__ == '__main__':
    args= arg_parsing()
    args["sim_type"] = "parameter_sweep_2D"
    default_setup=setup_default(args)
    sweep_2D_package=parameter_sweep_2D(default_setup)


    
