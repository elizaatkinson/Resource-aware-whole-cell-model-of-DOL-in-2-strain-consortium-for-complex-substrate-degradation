# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:38:17 2020
@author: Eliza
Scans a single parameter
Enabled for multiple processors if available
"""

#import functions
import time
import os
import pandas as pd
from functools import partial
from dataclasses import dataclass
from dask.distributed import Client, LocalCluster
import seaborn as sns
sns.set()
from setup_file import arg_parsing, setup_default, check_server, save_data, updateProgressBar, process_1D_sweep, results_to_plot, full_title_dict, plot_line_graph, print_fig
from run_ODE import run_sim_conditional

#custom dataclasses to save simulation data
@dataclass
class SimResultsP1:
    """Custom dataclass for 1D parameter sweep"""
    p1: int #value of the parameter
    solution: tuple #simulation results
    
class SweepResults1D:
    """Dataclass containing all the setup and solutions"""
    def __init__(self, data, setup):
        self.setup= setup
        self.data= data
    
    def process_data(self):
        self.results = process_1D_sweep(data=self.data, setup=self.setup)
        
    def save_data(self):
        dataname = "1D_paramsweep_" + self.setup.ODE_type + "_" + self.setup.p1_name
        save_data(data=self, dataname=dataname, setup=self.setup)

"""Basic sweep def and distribution"""

def param_scan_1D(setup, p1_value):
    """Run sim with specific parameter values - accepts all the inputs of the run_ODE file"""
    selected_parameter = setup.p1_name
    if selected_parameter in dict.keys(setup.p):
        #change value of parameter
        setup.p[selected_parameter] = p1_value
        #simulate
        sim_results = run_sim_conditional(default_setup=setup)
        return SimResultsP1(p1=p1_value, solution=sim_results.solution)
    else:
        raise ValueError("That parameter does not exist") 
        
def run_sweep_client(claster,setup):
    """Send sims to different processors"""
    # computation
    client = Client(claster)
    print(client)
    # partial
    f = partial(param_scan_1D, setup)
    A = client.map(f, setup.p1_sweep)
    data = client.gather(A)
    return(data)

"""Output data"""

def plot_results_v_param(results, setup):
    """Plotting results vs parameter values"""
    titles = full_title_dict(setup)
    molecs_to_plot= ["ea", "eb", "s0", "s1", "s", "lam"]
    for i in molecs_to_plot:
        y = results_to_plot(i=i, setup=setup, results=results)
        if len(y) == 3:
            labels=["cell a", "cell b", "average"]
        else:
            labels=False
        plot_line_graph(x=setup.p1_sweep, 
                        y=y,
                        title= setup.ODE_type.capitalize() + " " + titles[i],
                        xlabel= setup.p1_name,
                        ylabel= "amount of " + titles[i],
                        legendlabel=labels)
        #save figure
        figname = setup.ODE_type + " " + titles[i] + "vs" + setup.p1_name
        print_fig(figname= figname, setup=setup) 

def param_1D_sweep_to_excel(sweep_1D_package):    
    
    now = sweep_1D_package.setup.now
    dt_string = now.strftime("%Y.%m.%d_%H.%M.%S")
    folder_name = "Results-" + dt_string
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    dir_name = os.path.join(os.getcwd(), folder_name)
    excel_single_sim_name = "param_sweep_1D_" + sweep_1D_package.setup.ODE_type + "_" + dt_string + ".xlsx"
    file_path = os.path.join(dir_name, excel_single_sim_name)
    
    parameters_df = pd.DataFrame.from_dict(sweep_1D_package.setup.p, orient='index')
    initial_values_df = pd.DataFrame.from_dict(sweep_1D_package.setup.initial_values, orient='index')
    other_info = {'ODE_type': sweep_1D_package.setup.ODE_type, 
    'substrate': sweep_1D_package.setup.substrate, 
    'steady_state_start': sweep_1D_package.setup.setup_start,
    'p1':sweep_1D_package.setup.p1_name,
    'p1_sweep': sweep_1D_package.setup.p1_sweep}
    other_info_df = pd.DataFrame.from_dict(other_info, orient='index')
    
    writer = pd.ExcelWriter(file_path, engine='xlsxwriter')
    parameters_df.to_excel(writer, sheet_name='parameters')
    initial_values_df.to_excel(writer, sheet_name='initial states')
    other_info_df.to_excel(writer, sheet_name='input arguments')
    
    if sweep_1D_package.setup.ODE_type == "mono":
        for i in sweep_1D_package.results:
            results_df = pd.DataFrame(sweep_1D_package.results[i], index=sweep_1D_package.setup.p1_sweep)
            results_df.to_excel(writer, sheet_name=i)
    else:
        for i in sweep_1D_package.results.cell_a:
            results_df = pd.DataFrame(sweep_1D_package.results.cell_a[i], index=sweep_1D_package.setup.p1_sweep)
            results_df.to_excel(writer, sheet_name="cell a " + i)
        for i in sweep_1D_package.results.cell_b:
            results_df = pd.DataFrame(sweep_1D_package.results.cell_b[i], index=sweep_1D_package.setup.p1_sweep)
            results_df.to_excel(writer, sheet_name="cell b " + i) 
        for i in sweep_1D_package.results.avg: 
            results_df = pd.DataFrame(sweep_1D_package.results.avg[i], index=sweep_1D_package.setup.p1_sweep)
            results_df.to_excel(writer, sheet_name="average " + i) 
    writer.save()


"""Run sweep"""

def parameter_sweep_1D(default_setup):
    """Changes the value of a parameter over a range of values, runs a simulation
    for each value and returns the solutions"""
    
    #run sweep based on enviornment
    t = time.time()
    #on server run with the run_sweep_client function
    if check_server() == True:
        local = LocalCluster(n_workers=32)
        data=run_sweep_client(claster=local, setup=default_setup)
    #off server use a for loop
    else:
        data = []
        i=0
        complete = len(default_setup.p1_sweep)
        for p1 in default_setup.p1_sweep:
            simulation = param_scan_1D(setup=default_setup, p1_value = p1)
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
    sweep_1D_package = SweepResults1D(setup=default_setup, data=data)
    
    #step 4: process data
    sweep_1D_package.process_data()
    
    #save data
    sweep_1D_package.save_data()    
    param_1D_sweep_to_excel(sweep_1D_package)
    
    #plot data
    plot_results_v_param(results=sweep_1D_package.results, setup=sweep_1D_package.setup)

    return(sweep_1D_package)


if __name__ == '__main__':
    #get args
    args = arg_parsing()
    args["sim_type"] = "parameter_sweep_1D"
    #setup
    default_setup=setup_default(args)
    #get data
    sweep_1D_package = parameter_sweep_1D(default_setup)
