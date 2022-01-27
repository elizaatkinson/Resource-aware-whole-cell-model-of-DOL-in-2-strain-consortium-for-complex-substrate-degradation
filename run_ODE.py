# -*- coding: utf-8 -*-
"""
Created on Fri May  1 15:42:17 2020
@author: eliza
Run single simulation of co-culture system, measuring intracellular
concentrations over time.The transient dynamics of the intracellular molecules
are not biologically relevant: only the steady state values should be considered
"""

#setup
import time
import copy
import pandas as pd
import numpy as np
import os
import scipy.integrate
from collections import OrderedDict
import seaborn as sns
sns.set()
from setup_file import setup_default, arg_parsing, save_data, ProcessedResults, full_title_dict, results_to_plot, plot_line_graph, print_fig
from ODE_mono import dydt_mono
from ODE_duo import dydt_duo

class SimulationResults():
    def __init__(self, solution, setup):
        """Basic class just setup and solutions from solver"""
        self.setup= setup
        self.solution= solution
        
    def process_data(self):
        """Processes data into a more readable form"""
        self.processed_data = ProcessedResults(setup=self.setup, solution=self.solution)
        
    def store_first_run(self, first_run_results):
        """If running from steady state the results class can also record the
        first simulation results"""
        self.initial_run = first_run_results
        
    def save_data(self):
        """Saves data in a Results folder with date + time"""
        dataname=self.setup.ODE_type + "_results"
        save_data(data=self, dataname=dataname, setup=self.setup)        

"""Setup"""

def setup_on_glucose(setup):
    """Change initial values to 0 starch and a set amount of glucose"""
    #I have to use indexes here because the setup input could be a dict (if from zero) or from a list (if from ss)
    glc_setup = copy.deepcopy(setup)
    glc_setup.initial_values[-1] = setup.p["s_in"] #glucose
    glc_setup.initial_values[-3] = 0 #starch
    glc_setup.y0[-1] = setup.p["s_in"] #glucose
    glc_setup.y0[-3] = 0 #starch
    return(glc_setup)

def setup_on_starch(setup): # can be set up from default setup or ss setup
    """Change initial values to set amount of starch and 0 glucose"""
    starch_setup = copy.deepcopy(setup)
    starch_setup.initial_values[-1] = 0.0 #glucose
    starch_setup.initial_values[-3] = setup.p["s0_in"] #starch
    starch_setup.y0[-1]= 0.0
    starch_setup.y0[-3]= setup.p["s0_in"] 
    return(starch_setup)

def setup_from_ss(default_setup, first_run_solution):
    """Use the steady states from the previous glucose simulation as intitial
    values in a second simulation setup"""
    second_run_setup = copy.deepcopy(default_setup)
    if second_run_setup.ODE_type == "mono":
        second_run_setup.initial_values = OrderedDict(zip(second_run_setup.initial_values.keys(), first_run_solution.y[:,-1]))
        second_run_setup.y0= list(second_run_setup.initial_values.values())
        #the heterologous enzymes can be set to zero- since they are extracellular
        #second_run_setup.change_initial_values(new_initial_values= {"ea":0.0,"eb":0.0})
    elif second_run_setup.ODE_type == "duo":
        #not using usual dict method here because the two cells will have same molecules in different amounts
        second_run_setup.initial_values = list(first_run_solution.y[:,-1])
        second_run_setup.y0= list(first_run_solution.y[:,-1])
        #second_run_setup.y0[3] = 0.0 #ea
        #second_run_setup.y0[23]= 0.0 #eb 
    return(second_run_setup)



"""Run simulations"""

def run_on_glucose(glucose_setup):
    """Run simulation with glucose set to continuous - can add arg starch = "fedbatch" to make it non-constant"""
    #run first sim to steady state
    glucose_solution = run_sim(setup=glucose_setup, substrate="glucose") 
    #process results
    glucose_results = SimulationResults(solution=glucose_solution, setup=glucose_setup)
    return(glucose_results)

def run_on_starch(starch_setup):
    """Run simulation with starch as continuous - can add arg starch = "fedbatch" to make it non-constant"""
    #run sim
    starch_solution=run_sim(setup=starch_setup, substrate="starch") 
    #process results
    starch_results= SimulationResults(solution=starch_solution, setup=starch_setup)
    return(starch_results)

def run_sim(setup, *args, **kwds):
    """Solves the mono or duo version of the RHS ODEs and returns the solution-
    can accept arguments for glucose/starch : cont/non-cont/fed-batch: see ODE
    files for more specifics"""
    if setup.ODE_type == "mono":
        dydt_withps = lambda t,y: dydt_mono(t, y, setup.p, *args, **kwds)
    if setup.ODE_type == "duo":
        dydt_withps = lambda t,y: dydt_duo(t, y, setup.p, *args, **kwds)
    solution = scipy.integrate.solve_ivp(fun=dydt_withps, t_span=setup.solver_args.t_span, y0=setup.y0, method=setup.solver_args.method, rtol=setup.solver_args.rtol, atol=setup.solver_args.atol)
    return(solution)

def run_sim_conditional(default_setup):
    """Run sim for different experimental setups- default sim is from zero
    on starch. Arguments can change to run sim from zero on glucose, or
    first run on glucose with or without heterologous protein then second run
    from steady state on starch"""
    
    #just give default y0 values for now
    default_setup.y0_setup() 
    
    #### glucose from zero
    if default_setup.substrate == "glucose" and default_setup.setup_start == "zero":
        glucose_setup=setup_on_glucose(setup=default_setup)
        sim_results = run_on_glucose(glucose_setup=glucose_setup) 
        
    #### starch from zero
    elif default_setup.substrate == "starch" and default_setup.setup_start == "zero":
        starch_setup = setup_on_starch(setup=default_setup)
        sim_results = run_on_starch(starch_setup=starch_setup)
        
    #### running from steadys state
    elif default_setup.setup_start != "zero":
        #if running from steady state the first sim is set up on glucose
        glucose_setup=setup_on_glucose(setup=default_setup)
        
        # run on glucose first - without population to just get steady state values for single cell
        first_sim_on_glc = run_on_glucose(glucose_setup=glucose_setup)
            
        second_run_setup=setup_from_ss(default_setup=default_setup, first_run_solution=first_sim_on_glc.solution)
        
        #### running on glucose from steady state
        if default_setup.substrate == "glucose":
            sim_results = run_on_glucose(glucose_setup=second_run_setup)
        
        #### running on starch from steady state
        else: 
            ss_setup_on_starch = setup_on_starch(setup=second_run_setup)
            sim_results= run_on_starch(starch_setup=ss_setup_on_starch)
        
        sim_results.store_first_run(first_run_results=first_sim_on_glc)
    
    else:
        raise ValueError("Argument not recognised- check -s or -ss")
    
    return(sim_results)

"""Output data"""

def plot_molecs_v_time(results, setup):
    """Plotting intracellular molecules against time"""
    titles = full_title_dict(setup)
    #this is the only part of the code that needs to be changed manually:
        #determines which values are plotted
    molecs_to_plot= ["s0", "s1", "s", "ea", "eb", "lam"]
    for i in molecs_to_plot:
        y = results_to_plot(i=i, setup=setup, results=results.processed_data.y)
        if len(y) == 3:
            labels=["cell a", "cell b", "average"]
        else:
            labels=False
        ax= plot_line_graph(x=results.solution.t, 
                    y=y,
                    title= setup.ODE_type.capitalize() + " " + titles[i],
                    xlabel="time (mins)",
                    ylabel= "amount of " + titles[i], 
                    legendlabel=labels)
        #plot steady state
        plot_steady_state(setup=setup, i=i, ax=ax, results=results.processed_data)
        #save figure
        figname = setup.ODE_type + " " + titles[i]
        print_fig(figname= figname , setup=setup)  

def plot_steady_state(setup, i, ax, results):
    """Give steady state values on the graph of molecs/time"""
    if setup.ODE_type == "mono":
        if i == "lam":
            steadystate_value=results.ss[i]
            doubling_time= np.log(2)/(results.ss[i])
            steadystate_label = "Steady state: lambda= " + str(steadystate_value) + "\n doubling time= " + str(doubling_time)
            ax.text(0.98,0.02, steadystate_label, 
                    horizontalalignment="right", 
                    verticalalignment="bottom", 
                    transform=ax.transAxes)
        #for other molecule just label steady state
        else:
            steadystate_value=results.ss[i]
            steadystate_label= "Steady state= " + str(steadystate_value)
            ax.text(0.98,0.02, 
                    steadystate_label, 
                    horizontalalignment="right", 
                    verticalalignment="bottom", 
                    transform=ax.transAxes)
    if setup.ODE_type == "duo":
        if i == "lam":
            doubling_time_cell_a= np.log(2)/(results.ss.cell_a["lam"])
            doubling_time_cell_b= np.log(2)/(results.ss.cell_b["lam"])
            doubling_time_avg = (doubling_time_cell_a + doubling_time_cell_b) /2
            #just to double check the averages work out the same
            doubling_time_avg2 = np.log(2)/(results.ss.avg["lam"])
            print(f'doubling time from avg doubling time: (( ln2 / {results.ss.cell_a["lam"]} ) +  ( ln2 / {results.ss.cell_b["lam"]} )) / 2 = {doubling_time_avg}')
            print(f'doubling time from avg lam: ln2 / {results.ss.avg["lam"]} = {doubling_time_avg2}')
            steadystate_label = "Steady state doubling time= \ncell a:" + str(doubling_time_cell_a) + "\ncell b: " +str(doubling_time_cell_b) + "\naverage:" +str(doubling_time_avg)
            ax.text(0.98,0.02, 
                    steadystate_label, 
                    horizontalalignment="right", 
                    verticalalignment="bottom", 
                    transform=ax.transAxes)
        else:
            steadystate_value = results_to_plot(i, setup=setup, results=results.ss)
            if len(steadystate_value) == 3:
                steadystate_label = "Steady state= \ncell a:" + str(steadystate_value[0]) + "\ncell b: " +str(steadystate_value[1]) + "\naverage:" +str(steadystate_value[2])  
            else:    
                steadystate_label= "Average steady state= " + str(steadystate_value[0])
            ax.text(0.98,0.02, 
                    steadystate_label, 
                    horizontalalignment="right", 
                    verticalalignment="bottom", 
                    transform=ax.transAxes)

def single_sim_data_to_excel(sim_results):
        
    now = sim_results.setup.now
    dt_string = now.strftime("%Y.%m.%d_%H.%M.%S")
    
    folder_name =  "Results-" + dt_string
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    dir_name = os.path.join(os.getcwd(), folder_name)
    excel_single_sim_name = "runODE_" + sim_results.setup.ODE_type + "_" + dt_string + ".xlsx"
    file_path = os.path.join(dir_name, excel_single_sim_name)
    
    parameters_df = pd.DataFrame.from_dict(sim_results.setup.p, orient='index')
    y0_df = pd.DataFrame(sim_results.setup.y0)
    other_info = {'ODE_type': sim_results.setup.ODE_type, 
                  'substrate': sim_results.setup.substrate, 
                  'steady_state_start': sim_results.setup.setup_start}
    other_info_df = pd.DataFrame.from_dict(other_info, orient='index')
    
    writer = pd.ExcelWriter(file_path, engine='xlsxwriter')
    parameters_df.to_excel(writer, sheet_name='parameters')
    y0_df.to_excel(writer, sheet_name='y0')
    other_info_df.to_excel(writer, sheet_name='input arguments')

    if sim_results.setup.ODE_type == "mono":    
        y_values_df = pd.DataFrame.from_dict(sim_results.processed_data.y)
        ss_values_df = pd.DataFrame.from_dict(sim_results.processed_data.ss, orient='index')
        y_values_df.to_excel(writer, sheet_name="y solutions")
        ss_values_df.to_excel(writer, sheet_name="steady state results")
    else:
        cell_a_y_values_df = pd.DataFrame.from_dict(sim_results.processed_data.y.cell_a)
        cell_a_ss_values_df = pd.DataFrame.from_dict(sim_results.processed_data.ss.cell_a, orient='index')
        cell_a_y_values_df.to_excel(writer, sheet_name="cell a y solutions")
        cell_a_ss_values_df.to_excel(writer, sheet_name="cell a steady state results")
        
        cell_b_y_values_df = pd.DataFrame.from_dict(sim_results.processed_data.y.cell_b)
        cell_b_ss_values_df = pd.DataFrame.from_dict(sim_results.processed_data.ss.cell_b, orient='index')
        cell_b_y_values_df.to_excel(writer, sheet_name="cell b y solutions")
        cell_b_ss_values_df.to_excel(writer, sheet_name="cell b steady state results")
        
        ext_y_values_df = pd.DataFrame.from_dict(sim_results.processed_data.y.ext)
        ext_ss_values_df = pd.DataFrame.from_dict(sim_results.processed_data.ss.ext, orient='index')
        ext_y_values_df.to_excel(writer, sheet_name="ext y solutions")
        ext_ss_values_df.to_excel(writer, sheet_name="ext steady state results")
        
        avg_y_values_df = pd.DataFrame.from_dict(sim_results.processed_data.y.avg)
        avg_ss_values_df = pd.DataFrame.from_dict(sim_results.processed_data.ss.avg, orient='index')
        avg_y_values_df.to_excel(writer, sheet_name="avg y solutions")
        avg_ss_values_df.to_excel(writer, sheet_name="avg steady state results")
   
    writer.save()
    
    
"""Run single ODE simulation"""
    
def run_ODE(default_setup):
    "Master run_ODE function that can be called in master_file"
    t = time.time() 
    #run sim with different args
    sim_results= run_sim_conditional(default_setup=default_setup)
    #process data
    sim_results.process_data()
    #save data
    sim_results.save_data()  
    #single_sim_data_to_excel(sim_results)
    elapsed = time.time() - t
    print(f'Simulation took {elapsed} seconds')
    #save important results to excel
    single_sim_data_to_excel(sim_results)
    #plot results
    plot_molecs_v_time(results=sim_results, setup=sim_results.setup)
    #return the processed results
    return(sim_results)
    
if __name__ == '__main__':
    #get start point and substrate from argparse
    args = arg_parsing()
    args["sim_type"] = "run_ODE"
    default_setup=setup_default(args)
    sim_results=run_ODE(default_setup)

