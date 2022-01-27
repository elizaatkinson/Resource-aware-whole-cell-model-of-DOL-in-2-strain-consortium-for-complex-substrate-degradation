# -*- coding: utf-8 -*-
"""
Created on Sun May 24 15:01:51 2020

@author: eliza

Setup file for DOL model
"""

import pandas as pd
from dataclasses import dataclass
from collections import OrderedDict
from datetime import datetime
import numpy as np
import os
import re
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
import argparse

class RunConfig():
    def __init__(self):
        """"Values to check/change to ensure sims are running correctly.
        Error messages within this file draw the "correct" values from this 
        class- if more species or parameters are added to the sim at any point
        these values may need to be changed manually."""
        self.no_of_params = 49                  #number of parameters
        self.no_of_initial_values = 27          #number of unique species
        self.cell_a_size = 21                #number of species in cell a (duo only)
        self.cell_b_size = 21                #number of species in cell b (duo only)
        self.ext_size = 3                    #number of extracellular species (duo only)
        self.init_vals_in_duo = (self.cell_a_size 
                            + self.cell_b_size 
                            + self.ext_size)     #number of species in duo full
        #all species names
        self.mono_keys = ["ma", "rma", "ea",
                          "mb", "rmb", "eb",
                          "mq", "rmq", "q",
                          "mm", "rmm", "em",
                          "mt", "rmt", "et",
                          "mr", "rmr", "r",
                          "si", "a", 
                          "mc", "rmc", "ec",
                          "z",
                          "s0", "s1", "s"]  #shorthand names expected
        self.full_names= ["enzyme A mRNA",
                          "ea mRNA:ribosome complex",
                          "enzyme A (protein)",
                          "enzyme B mRNA",
                          "eb mRNA:ribosome complex",
                          "enzyme B (protein)",
                          "housekeeping mRNA",
                          "housekeeping mRNA:ribosome complex",
                          "housekeeping proteins",
                          "metabolic mRNA",
                          "metabolic mRNA:ribosome complex",
                          "metabolic enzymes",
                          "transport mRNA",
                          "transport mRNA:ribsome complex",
                          "transport proteins",
                          "ribosomal mRNA",
                          "ribosomal mRNA:ribosome complex",
                          "free ribosomes",
                          "imported substrate",
                          "available energy",
                          "enzyme C mRNA",
                          "ec mRNA:ribosome complex",
                          "enzyme C (protein)",
                          "product",
                          "starch",
                          "maltodextrins",
                          "glucose",
                          "growth rate"] #full names
        #which species are only found in one or the other
        self.cell_a_only = ["ma", "rma", "ea"]
        self.cell_b_only = ["mb", "rmb", "eb"]
        self.ext_only = ["s0", "s1", "s"]
        #the name of the server profile
        self.server_profile = '/home/ea2915'
        #name of the input file
        self.input_csv = "DOL_model_input_file.csv"

"""
External inputs
"""

def arg_parsing():
    """Argument parsing for simulations.
    The master file will either run a new sim or read data from old sims
    All sims require an ODE_type (-t) mono or duo argument- this will cause 
    the mono or duo versions of the RHS ODEs to be run.
    1D or 2D sweeps require parameter inputs: -p1/2n provides the name, -p1/2s
    provides the values of the sweep."""
    ap = argparse.ArgumentParser()
    #can input a pre-run simulation and check values
    ap.add_argument("-read", "--data_file", default="none", help="read an existing data file")
    #compare data from 2 sims(this is not actually implemented yet)
    ap.add_argument("-compare", "-comp_files", default="none")    
    #sim type - can just run one sim, or either a 1D or 2D parameter sweep
    ap.add_argument("-run", "--sim_type", default="none", help="run_ODE, parameter_sweep_1D or parameter_sweep_2D")
    #ODE type - monoculture or duo-culture
    ap.add_argument("-t", "--ODE_type", default="none", help="mono or duo or both")
    #run from zero or steady state: zero is default
    ap.add_argument("-ss", "--steady_state", default="zero", help="start from steady state on glucose")
    #may want to run on glucose or starch: default is starch
    ap.add_argument("-s", "--substrate", default="starch", help="starch or glucose")
    #require sweep arguments depending on sweep type
    ap.add_argument("-p1n", "--p1_name", default="none", help="name of the 1st parameter")
    ap.add_argument("-p1s", "--p1_sweep", default="none", help="parameter 1 sweep values start:end:steps format: \d+:\d+:\d+")
    ap.add_argument("-p2n", "--p2_name", default="none", help="name of the 2nd parameter")
    ap.add_argument("-p2s", "--p2_sweep", default="none", help="parameter 2 sweep valies start:end:steps format: \d+:\d+:\d+")
    args = vars(ap.parse_args())
    return(args) 

@dataclass
class SolverArgs:
    "Class for the arguments for the ODE solver"
    t_span: tuple
    method: str
    rtol: float
    atol: float

@dataclass
class Inputfile:
    "Class for all the inputs from the csv file"
    p: dict
    initial_values: OrderedDict or list
    solver_args: SolverArgs

def CSV_reader(dataname):
    """Reads the csv file as a pandas dataframe- collates the parameters,
    initial values and arguments for the solver into a custom class"""
    with open(dataname) as csv_file:
        df = pd.read_csv(csv_file)
        p_names = df["parameter shorthand name"].astype(str)
        p_values = df["parameter value"].astype(float)
        p = dict(zip(p_names, p_values))
        if "nan" in p:
            del p["nan"]
        initial_names = df["initial states shorthand name"].astype(str)
        initial_data = df["initial states value"].astype(float)
        initial_values = OrderedDict(zip(initial_names, initial_data))
        if "nan" in initial_values:
            del initial_values["nan"]
        t_span = (0, float(df["time period end"][0]))
        method = str(df["method"][0])
        rtol = float(df["rtol"][0])
        atol = float(df["atol"][0])            
        solver_args=SolverArgs(t_span=t_span, method=method,rtol=rtol,atol=atol)
    return(Inputfile(p=p, initial_values=initial_values, solver_args=solver_args))

"""Setup for simulations"""

def setup_default(args):
    """The default setup is generated from the CSV file and the arguments
    provided by argparse"""
    #get parameters, initial values, time span and method from input file
    input_file = CSV_reader(RunConfig().input_csv)
    #put in standard setup format
    setup= Setup(args=args, input_file=input_file)
    return(setup)

@dataclass
class Setup:
    """Base level setup for every simulation.
    The input_file should contain the parameters, initial values and args for 
    intergation solver."""
    def __init__(self, input_file, args):
        #record what kind of sim
        self.sim_type = args["sim_type"]
        assert self.sim_type == "run_ODE" or self.sim_type == "parameter_sweep_1D" or self.sim_type == "parameter_sweep_2D"
        #get ODE type from args
        self.ODE_type = args["ODE_type"]
        assert self.ODE_type == "mono" or self.ODE_type == "duo" or self.ODE_type == "both", "ODE type should be mono or duo or both"
        #determines if running from zero or steady state
        self.setup_start = args["steady_state"]
        #determines if second (or only) sim is run on glucose or starch
        self.substrate = args["substrate"]
        #get parameters from csv
        self.p= input_file.p 
        assert len(self.p) == RunConfig().no_of_params, "incorrect number of parameters"
        #get initial values from csv
        self.initial_values = input_file.initial_values
        assert len(self.initial_values) == RunConfig().no_of_initial_values, "incorrect number of initial values"
        assert list(self.initial_values.keys()) == RunConfig().mono_keys, "incorrect keys for initial values"
        #tspan, method, rtol, atol
        self.solver_args= input_file.solver_args
        #get date & time of simulation
        self.now = datetime.now()
        #if doing a sweep get the args for p1/p2
        if args["sim_type"] == "parameter_sweep_1D" or args["sim_type"] == "parameter_sweep_2D":
            setup_sweep(self, args)
        
    def change_parameters(self, new_parameters):
        """Changes any parameters after the initial setup"""
        for new_p, p_value in new_parameters.items():
            if new_p in self.p:
                self.p[new_p] = p_value
            else:
                raise ValueError("parameter not found")

    def change_initial_values(self, new_initial_values):
        """Changes any of the initial values after the initial setup"""
        for new_init, i_value in new_initial_values.items():
            if new_init in self.initial_values:
                self.initial_values[new_init] = i_value
            else:
                raise ValueError("initial value not found") 
                
    def change_solver_args(self, new_t_span="none", new_method="none", new_rtol="none", new_atol="none"):
        """Changes t_span, method or tolerances for the ODE solver"""
        #change t_span
        if new_t_span != "none":
            self.solver_args.t_span = new_t_span
        elif new_method != "none":
            self.solver_args.method = new_method
        elif new_rtol != "none":
            self.solver_args.rtol = new_rtol      
        elif new_atol != "none":
            self.solver_args.atol = new_atol   
        else: 
            raise ValueError("must provide a value to change")
        
    def y0_setup(self):
        """Sets up y0 values from the initial values dictionary.
        In a mono case the initial values will just be taken from the setup
        class and converted into a list for the solver.
        In a duo case the initial values must be arraged into the order of 
        expected values from both cells and extracellular environment"""
        if self.ODE_type == "mono":
            self.y0 = list(self.initial_values.values())
            assert len(self.y0) == RunConfig().no_of_initial_values, "unexpected number of y0 values for mono"
        elif self.ODE_type == "duo":
            self.y0 = duo_reorder(setup=self, mono_order=list(self.initial_values.values()))
            assert len(self.y0) == RunConfig().init_vals_in_duo, "unexpected number of y0 values for duo"

  
def duo_reorder(setup, mono_order):
    """Rearranges a list of values or keys from the values provided in the
    setup function, into the required number of values in the duo form 
    RHS ODE.
    Some molecules are only found in one of the two cells or in the 
    extracellular environment- these molecules are defined in the RunConfig function
    at the beginning of this file, and will need to be changed manually if
    the RHS ODEs change"""
    mono_names = RunConfig().mono_keys
    cell_a=[]
    cell_b=[]
    ext=[]
    for i in mono_names:
        #looks at the RunConfig file- which tells the code which shorthand names belong to each compartment
        #appends value to correct compartment
        if i in RunConfig().cell_a_only:
            cell_a.append(mono_order[mono_names.index(i)])
        elif i in RunConfig().cell_b_only:
            cell_b.append(mono_order[mono_names.index(i)])
        elif i in RunConfig().ext_only:
            ext.append(mono_order[mono_names.index(i)])
        else:
            cell_a.append(mono_order[mono_names.index(i)])
            cell_b.append(mono_order[mono_names.index(i)])
    #double check that all the lengths are correct
    assert len(cell_a) == RunConfig().cell_a_size
    assert len(cell_b) == RunConfig().cell_b_size
    assert len(ext) == RunConfig().ext_size
    return(cell_a + cell_b + ext)

def setup_sweep(self, args):
    """Processes parameter-sweep-specific args from argparse.
    The parameter names are taken as strings- will flag up an error if the
    str does not exist in the dict keys for the parameters.
    The sweeps are accepted as strings in the form start:end:steps and converted
    into floats. Some error messages are flagged if numbers are not logical"""
    #for both 1d and 2d sweeps:
    #set parameter 1 name from argparse, and check name exists
    self.p1_name = args["p1_name"]
    assert self.p1_name in self.p, "-p1n parameter does not exist" 
    #check p1 sweep was given in correct format
    assert re.match('\d*\.?\d*:\d*\.?\d*:\d*\.?\d*', args["p1_sweep"]), "-p1s incorrect format"
    #convert into 3 inputs
    p1_start, p1_end, p1_steps = args["p1_sweep"].split(":")
    #check values are logical
    assert float(p1_start) < float(p1_end), "end point must be after start point"
    assert float(p1_steps) < float(p1_end), "steps must be smaller than total range"
    #store p1_sweep in tuple
    self.p1_sweep = np.arange(float(p1_start), float(p1_end), float(p1_steps))
    
    #for only 2d sweeps:
    if args["sim_type"] == "parameter_sweep_2D":
        #set name, check it exists
        self.p2_name = args["p2_name"]
        assert self.p2_name in self.p, "-p2n parameter does not exist"
        #check sweep in correct format
        assert re.match('\d*\.?\d*:\d*\.?\d*:\d*\.?\d*', args["p2_sweep"]), "-p2s incorrect format"     
        #convert into 3 imputs, check values are logical
        p2_start, p2_end, p2_steps = args["p2_sweep"].split(":")
        assert float(p2_start) < float(p2_end), "end point must be after start point"
        assert float(p2_steps) < float(p2_end), "steps must be smaller than total range"
        #store p2_sweep in tuple
        self.p2_sweep = np.arange(float(p2_start), float(p2_end), float(p2_steps))
        #make a list of every combination of p1 and p2 value
        self.comb_param_values = list(itertools.product(self.p1_sweep,self.p2_sweep))


"""
Process single simulations
"""

class ProcessedResults():
    def __init__(self, setup, solution):
        """Processing solutions for run_ODE.
        For ease of plotting the t solution is stored in results alongside
        the y solutions.
        The y solutions are assigned to the correct shorthand names for their
        respective molecules.
        Steady state is calculated and stored"""
        #assign simulation to molecule
        self.y=assign_names(setup=setup, results_values=solution.y)
        #calculate growth rate
        get_growth_rate(results=self.y, setup=setup)
        #create separate section for the steadt states
        self.ss = get_steady_state(solution=solution, setup=setup)  
        
def results_dict(results_keys, results_values):
    """Zip function for keys and values, with error flagging
    This is just the basic zip function, but for all uses in this code, I want
    to ensure the keys=values, so this function just adds in some error flagging"""
    #check that keys and values are equal- else the zip function won't work
    if len(results_keys) != len(results_values):
        if len(results_keys) > len(results_values):
            raise ValueError("more keys than solutions")
        elif len(results_keys) < len(results_values):
            raise ValueError("more solutions than keys")
    else:
        #zips keys to values to make Ordered Dict
        return(OrderedDict(zip(results_keys, results_values)))
    
def assign_names(setup, results_values):
    """Creates dictionaries with the shorthand names for each molecule and the
    values corresponding to that molecule- may be solutions to a single sim,
    the steady state values, the parameter solutions etc."""
    if setup.ODE_type == "mono":
        results = results_dict(results_keys=RunConfig().mono_keys, results_values=results_values)
    elif setup.ODE_type == "duo":
        results = Duo_results_dict(setup=setup, values=results_values)
    else:
        raise ValueError("ODE_type must be either mono or duo")
    get_growth_rate(results=results, setup=setup)
    return(results)

class Duo_assigned():
    """Divides keys or values into the subregions of a consortia system.
    The output of a duo form simulation will give cell a, b and extracellular
    environment in one list. If we wish to assign these solutions to the correct
    molecules in a dictionary they need to be subdivided, since there will
    be repeating molecules in the two cells"""
    def __init__(self, duo_order):
        #duo_order = the values spit out by the solver which come in
        #cell a then cell b then extracellular
        cell_a_size= RunConfig().cell_a_size
        cell_b_size= RunConfig().cell_b_size
        self.cell_a= duo_order[:cell_a_size]
        self.cell_b= duo_order[cell_a_size:(cell_a_size+cell_b_size)]
        self.ext= duo_order[(cell_a_size+cell_b_size):]

class Duo_results_dict():
    """Results class for duo simulations
    Assigns shorthand names to each of the results, subdivided based on the cell
    or extracellular environment"""
    def __init__(self, setup, values):
        #use monoculture shorthand names and reorder so they are repeated for cell b
        y_names = duo_reorder(setup=setup, mono_order=RunConfig().mono_keys)
        #assign shorthand names to each compartment
        y_keys = Duo_assigned(duo_order=y_names)
        #assign solutions to compartments
        y_values = Duo_assigned(duo_order=values)
        #zip up each compartment as its own dictionary
        self.cell_a = results_dict(results_keys=y_keys.cell_a, results_values=y_values.cell_a)
        self.cell_b = results_dict(results_keys=y_keys.cell_b, results_values=y_values.cell_b)
        self.ext = results_dict(results_keys=y_keys.ext, results_values=y_values.ext)
        #get average of the two cells
        self.avg = get_avg(self=self, setup=setup)
        
"""Calculations"""

def get_growth_rate(results, setup):
    """Calculate growth rates- uses the Weisse equation for growth to calculate 
    growth rate at each point from each result array"""
    if setup.ODE_type == "mono":
        results["lam"] = (results["rmq"] + results["rmr"] + results["rmt"] + results["rmm"] + results["rma"] + results["rmb"] + results["rmc"])*(setup.p["gmax"]*results["a"]/(setup.p["Kgamma"] + results["a"])/setup.p["M"])
    elif setup.ODE_type == "duo":
        results.cell_a["lam"] = (results.cell_a["rmq"] + results.cell_a["rmr"] + results.cell_a["rmt"] + results.cell_a["rmm"] + results.cell_a["rma"] + results.cell_a["rmc"])*(setup.p["gmax"]*results.cell_a["a"]/(setup.p["Kgamma"] + results.cell_a["a"])/setup.p["M"])
        results.cell_b["lam"] = (results.cell_b["rmq"] + results.cell_b["rmr"] + results.cell_b["rmt"] + results.cell_b["rmm"] + results.cell_b["rmb"] + results.cell_b["rmc"])*(setup.p["gmax"]*results.cell_b["a"]/(setup.p["Kgamma"] + results.cell_b["a"])/setup.p["M"])
        results.avg["lam"] = (results.cell_a["lam"] + results.cell_b["lam"]) / 2
    else:
        raise ValueError("ODE_type must be either mono or duo")
    return(results)

def get_steady_state(solution, setup):
    """Get all steady state values (for use in further sims)- creates a dictionary
    (or multiple for duo) with correct assigned molecule"""
    ss_values = []
    for i in solution.y:
        ss_values.append(i[-1])
    if setup.ODE_type == "mono":
        steady_state = results_dict(results_keys=RunConfig().mono_keys, results_values=ss_values)
        get_growth_rate(results=steady_state, setup=setup)
    elif setup.ODE_type == "duo":
        steady_state = Duo_results_dict(setup=setup, values=ss_values)
        get_growth_rate(results=steady_state, setup=setup)
    else:
        raise ValueError("ODE_type must be either mono or duo")
    return(steady_state)

def get_avg(self, setup):
    """Getting averages of cells in duo system"""
    avg_names=[]
    avg_values=[]
    for name in RunConfig().mono_keys:
        #if the shorthand name is found in both a & b- then an average can be found
        if name in self.cell_a and name in self.cell_b:
            avg_names.append(name)
    for i in avg_names:
        avg = (self.cell_a[i] + self.cell_b[i]) / 2
        avg_values.append(avg)
    assert len(avg_names) == len(avg_values), "unexpected average values"
    average_dict = results_dict(results_keys=avg_names, results_values=avg_values)
    return(average_dict)

"""
Process sweeps
"""

def process_1D_sweep(data, setup):
    """Process 1D sweeps"""
    #initialise arrays
    index_sols = np.arange(0, len(data[0].solution.y), 1)
    sol_arrays= []
    for i in index_sols:
        array=np.empty((len(setup.p1_sweep)))
        array[:]=np.NaN
        sol_arrays.append(array)
    #populate arrays with results
    for datum in data:
        norm_p1= abs(datum.p1 - setup.p1_sweep)
        idx_p1=np.argmin(norm_p1)
        for i in index_sols:
            sol_arrays[i][idx_p1]=datum.solution.y[i][-1]
    #assign arrays to correct molecules
    sweep_results=assign_names(setup=setup, results_values=sol_arrays)
    return(sweep_results)

def process_2D_sweep(data, setup):
    """Process 2D sweeps"""
    #initialise arrays
    index_sols = np.arange(0, len(data[0].solution.y), 1)
    sol_arrays= []
    for i in index_sols:
        array=np.empty((len(setup.p1_sweep),len(setup.p2_sweep)))
        array[:]=np.NaN
        sol_arrays.append(array)
    #populate arrays with results
    for datum in data:
        norm_p1= abs(datum.p1 - setup.p1_sweep)
        idx_p1=np.argmin(norm_p1)
        norm_p2= abs(datum.p2 - setup.p2_sweep)
        idx_p2=np.argmin(norm_p2)
        for i in index_sols:
            sol_arrays[i][idx_p1][idx_p2]=datum.solution.y[i][-1]
    #assign arrays to correct molecules
    sweep_results=assign_names(setup=setup, results_values=sol_arrays)
    return(sweep_results)       
    
"""
Plotting figures
"""

def full_title_dict(setup):
    """Set up titles for plotting""" 
    #get shorthand names
    dict_keys = RunConfig().mono_keys
    dict_keys.append("lam")
    #assign each molecule a full title
    dict_values = RunConfig().full_names
    titles = results_dict(results_keys=dict_keys, results_values=dict_values)
    return(titles)

def results_to_plot(i, setup, results):
    """For duo culture, results are stored in different subclasses. If a molecule is found in both cells it will plot both 
    cells and the average- on line graphs this will be plotted on the same graph,
    heatmaps will be 3 separate graphs"""
    #for mono cells only one dictionary
    if setup.ODE_type == "mono":
        y=[results[i]]
    if setup.ODE_type == "duo":
        if i in results.avg:
            y=[results.cell_a[i], 
            results.cell_b[i],
            results.avg[i]]
        elif i in results.cell_a:
            y=[results.cell_a[i]]
        elif i in results.cell_b:
            y=[results.cell_b[i]]
        elif i in results.ext:
            y=[results.ext[i]]
        else:
            raise ValueError ("Results not found")
    return(y)

def plot_line_graph(x, y, title, xlabel, ylabel, legendlabel=False):
    """Plot line graph- basic matplotlib setup"""
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for ys in y:
        ax.plot(x, ys)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if legendlabel != False:
        ax.legend(legendlabel)
    return(ax)

#tick marks for heatmaps
def tick_marks(xaxis, yaxis):
    """Making x and y axis tick marks for heatmaps- seaborn's automatic
    labelling leaves something to be desired"""
    x_axis= []
    y_axis= []
    if len(xaxis) > 20 or len(yaxis) > 20:
        
        major_ticks_x = 10
        major_ticks_y = 10
        for i in xaxis:
            if i % major_ticks_x == 0:
                #if it is the desired step value it is added to the list
                x_axis.append(str(int(i)))
            else:
                #if it is not the desired step value it is made blank
                x_axis.append("")
        for i in yaxis:
            if i % major_ticks_y == 0:
                y_axis.append(str(int(i)))
            else:
                y_axis.append("")
                
    else:
        for x in xaxis:
            x_axis.append(str(int(x)))
        for y in yaxis:
            y_axis.append(str(int(y)))
    return(x_axis, y_axis)

  
def plot_heatmap(data, title, xaxis, yaxis, xlabel, ylabel, cbarlabel="none", *args, **kwds):
    """Plot heatmap using seaborn"""
    tickmarks=tick_marks(xaxis, yaxis)
    xticklabels=tickmarks[0]
    yticklabels=tickmarks[1]
    sns.set()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    if cbarlabel == "none":
        ax = sns.heatmap(data= data, xticklabels=xticklabels, yticklabels=yticklabels, *args, **kwds) 
    else:
        ax = sns.heatmap(data= data, xticklabels=xticklabels, yticklabels=yticklabels, cbar_kws={'label': cbarlabel}, *args, **kwds) 
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
        
"""
#Supporting functions
"""

def check_server():
    """check if the working directory is in the server or not 
    (determines how many processors can be used)"""
    #my windows laptop and the server have different "HOME"s,
    #so it tries both
    try:
        home=os.environ['HOME']
        if home==RunConfig().server_profile:
            return(True)
        else:
            return(False)
    except KeyError:
        home=os.environ['HOMEPATH']
        if home==RunConfig().server_profile:
            return(True)
        else:
            return(False)

 
def save_data(data, dataname, setup):
    """Function for saving data as a pickle file- will save as the custom
    results class relevant to the type of simulation. Saves in a Results file
    with the date and time"""
    #get date and time
    now = setup.now
    dt_string = now.strftime("%Y.%m.%d_%H.%M.%S")
    #check if folder for date already exists
    folder_name =  "Results-" + dt_string
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    dir_name = os.path.join(os.getcwd(), folder_name)
    file_name = dataname + "_" + dt_string + ".pickle"
    file_path = os.path.join(dir_name, file_name)
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)


def load_data(dataname, *args, **kwds):
    """Function for loading data from a pickle file in the current
    working directory"""
        #get date and time
    dir_name = os.getcwd()
    file_path = os.path.join(dir_name, dataname)
    with open(file_path, 'rb') as f:
        return(pickle.load(f))
    
def read_file(dataname):
    """Opens a pickle file and reads out the important features- can
    be used to check the values used in specific Results folders"""
    data=load_data(dataname)
    print(f'sim_type == {data.setup.sim_type}')
    print(f'ODE_type == {data.setup.ODE_type}')
    print(f'parameter values == {data.setup.p}')
    print(f'initial values == {data.setup.initial_values}')
    print(f'substrate == {data.setup.substrate}')
    print(f'steady state start? == {data.setup.setup_start}')
    return(data)

def print_fig(figname, setup):
    """Saving figures as png with high resolution"""
    #get date and time
    now = setup.now
    dt_string = now.strftime("%Y.%m.%d_%H.%M.%S")
    #check if folder for date already exists
    folder_name =  "Results-" + dt_string
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    #save
    dir_name = os.path.join(os.getcwd(), folder_name)
    file_name = figname + "_" + dt_string + ".png"
    file_path= os.path.join(dir_name, file_name)
    plt.savefig(file_path, dpi=1200)    #dpi sets the resolution
    plt.show()
   
#progress bar - function found on the internet
barWidth = 50 
def updateProgressBar(value):
    line = '\r%s%%[%s]' % ( str(value).rjust(3), '-' * int ((float(value)/100) * barWidth))  
    print(line, end='')
    sys.stdout.flush()

