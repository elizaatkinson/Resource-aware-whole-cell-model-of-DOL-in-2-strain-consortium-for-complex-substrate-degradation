# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 14:35:59 2020

@author: eliza

Master file for all simulation types
"""

from setup_file import read_file, arg_parsing, setup_default
from run_ODE import run_ODE, SimulationResults, plot_molecs_v_time
from parameter_sweep_1D import parameter_sweep_1D, SweepResults1D, plot_results_v_param
from parameter_sweep_2D import parameter_sweep_2D, SweepResults2D, plot_2D_results_v_params

if __name__ == '__main__':
    
    #get arguments to determine what to do
    args = arg_parsing()
    
    
    #if no pickle file is provided- the code will run a simulation
    if args["data_file"] == "none":
        
        ##convert arguments into setup information for simulation
        setup = setup_default(args)
        
    
        #run a single simulation
        if args["sim_type"] == "run_ODE":
            if args["ODE_type"] == "mono" or args["ODE_type"] == "duo":
                sim_results = run_ODE(setup)
            elif args["ODE_type"] == "both":
                setup.ODE_type = "mono"
                mono_results= run_ODE(setup)
                setup.ODE_type  = "duo"
                duo_results = run_ODE(setup)
            else:
                raise ValueError ("no such ODE_type")
           
        #run a 1D parameter sweep        
        elif args["sim_type"]== "parameter_sweep_1D":
            assert args["p1_name"] != "none" and args["p1_sweep"] != "none", "parameter sweep requires -p1n and -p1s argument"
            if args["ODE_type"] == "mono" or args["ODE_type"] == "duo":
                sim_results = parameter_sweep_1D(setup)
            elif args["ODE_type"] == "both":
                setup.ODE_type = "mono"
                mono_results= parameter_sweep_1D(setup)
                setup.ODE_type  = "duo"
                duo_results = parameter_sweep_1D(setup)
            else:
                raise ValueError ("no such ODE_type")
        
        #run a 2D parameter sweep
        elif args["sim_type"]== "parameter_sweep_2D":
            assert args["p1_name"] != "none" and args["p1_sweep"] != "none", "parameter sweep requires -p1n and -p1s argument"
            assert args["p2_name"] != "none" and args["p2_sweep"] != "none", "parameter sweep requires -p2n and -p2s argument"
            if args["ODE_type"] == "mono" or args["ODE_type"] == "duo":
                sim_results = parameter_sweep_2D(setup)
            elif args["ODE_type"] == "both":
                setup.ODE_type = "mono"
                mono_results= parameter_sweep_2D(setup)
                setup.ODE_type = "duo"
                duo_results = parameter_sweep_2D(setup)    
            else:
                raise ValueError ("no such ODE_type")
            
        else:
            raise ValueError("invalid simulation type- ensure -run argument is valid")
       
    #if a pickle file is provided the code will read it and plot the default graphs for the sim type        
    else:
        
        sim_results = read_file(dataname=args["data_file"])
        if sim_results.setup.sim_type == "run_ODE":
            plot_molecs_v_time(results=sim_results, setup=sim_results.setup)
        elif sim_results.setup.sim_type == "parameter_sweep_1D":
            plot_results_v_param(results=sim_results.results, setup=sim_results.setup)
        elif sim_results.setup.sim_type == "parameter_sweep_2D":
            plot_2D_results_v_params(results=sim_results.results, setup=sim_results.setup)
        else:
            raise ValueError ("no sim type provided")