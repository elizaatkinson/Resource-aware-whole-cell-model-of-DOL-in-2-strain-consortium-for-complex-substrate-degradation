# DOL_model

A host aware model to compare starch degradation in a monoculture system and a consortium (duo-culture) system. 

## Biological background 
Starch degradation requires a-amylase and glucoamylase, which can be co-expressed in a single cell or each expressed in a different cell, forming a 2 strain consortium. This model is used to investigate the difference in metabolic burden, and the difference in resource allocation, experienced by cells in the monoculture and in the consortium.  

## Getting started

### Understanding the files

**master_file.py:** The main driver for any simaultion type: can read existing data files or run original simulations. This file gives access to all the features of the model. 

**ODE files:** ODE_mono.py contains all the RHS ODEs for each molecule in the monoculture system. ODE_duo.py contains the RHS ODEs for the duo-culture system. If reactions need to be added to the model they can be added here.

**setup_file.py:** This contains all the functions which will be shared between the different simulation types. If reactions need to be added to the model, the RunConfig() class within the setup file will need to be altered to account for the new reactions. 

**Simulation files:** run_ODE.py, parameter_sweep_1D.py and parameter_sweep_2D.py contain functions specific to the different simulation types.

**compare.py:** This will run comparisons between the same simulation done for monoculture and duoculture and produce graphs for the differences.

### Creating a CSV file input

The DOL_model_input_file.csv file contains the default values for the parameters, initial_values and arguments for the ODE solver. These values can be changed manually within this file. 
Do not change the file name, or any of the full or shorthand names, or it will not be read correctly. If new parameters or initial values need to be added, ensure the correct changes are made in the ODE files and setup_file.py first.  

### Input arguments

#### -run (choose which simulation to run)

To run a simulation from scratch use the -run argument. The run argument accepts three inputs:
* run_ODE : runs a single simulation measuring the change in intracellular molecule over time
* parameter_sweep_1D: varies a single parameter: running a simulation to steady state for each value of the parameter
* parameter_sweep_2D: varies 2 parameters, running a simulation to steady state for each pair of values

#### -t (choose which type to run)

The model can run the same parameters for a monoculture or a consortium (duoculture) case. The type argument accepts:
* mono : runs the monoculture case where both cells co-express both amylases
* duo : runs the duoculture case where one cell expresses a-amylase and the other expresses glucoamylase
* both : will run mono followed by duo - not advised for large parameter sweeps as it requires high memory)

#### -s (choose which substrate to provide)

Cells can be provided, in a chemostat, with glucose or starch as the primary substrate. Substrate argument accepts two inputs:
* glucose : glucose is supplied as the substrate
* starch : starch is supplied as the substrate

#### -ss (choose whether to run from zero or from steady state)

If this argument is not provided the model will run using the initial states given in the input CSV file. Steady state argument accepts one input:
* w : the model will first run a simulation to steady state on glucose, then it will initialise a run on starch with those steady states

## Run a single simulation to steady state

run_ODE runs a single simulation of molecules over time with the parameters, initial values and solver arguments provided in the csv file. It requires an ODE_type argument (-t) to run either the mono RHS ODEs or the duo RHS ODEs, or both. The biologically-meaningful result from this simulation is the steady state value which represents the amount of that molecule present in the average cell in exponential growth. 

To run the mono version of the RHS ODEs on glucose:

```
python3 master_file.py -run run_ODE -t mono -s glucose
```
This will output a new Results folder in the current working directory. This Results folder will contain a pickle file and an excel file, which both store the results of the simulation, and png files of any plotted graphs.  

To run the duo version of the RHS ODEs:

```
python3 master_file.py -run run_ODE -t duo -s glucose
```

To automatically run both mono and duo versions:
```
python3 master_file.py -run run_ODE -t both -s glucose
```

***Running from steady state*** : (steady state (ss) argument)

In some cases running from zero values will not support the growth of the cell on starch. Instead the cell should be growth on glucose to steady state, and then those steady state values used in a second simulation for growth on starch. The run_ODE code can provide two different steady state starting points: one without expression of heterologous protein (wo), and the other with (w) the desired level of heterologous protein expression. The run_ODE file will run the first simulation from the initial values provided to steady state on glucose, then run a second simulation on starch, providing the mono_results file as before. 

To run with heterologous protein in the first simulation on glucose, then run the second simulation with the same heterologous protein expression on starch use the -ss w argument: 
```
python3 master_file.py -run run_ODE -t mono -s starch -ss w
```

The default value for -ss is "zero"- which will just run a single simulation from the intitial values provided. 


## Parameter sweeps
Parameter sweeps return the steady state values achieved by the cell across a range of values for a single, or two prarameters. To find the correct name to input for each parameter, see the "parameter shorthand name" column in the CSV file. 

**1D_parameter_sweep** : run sweeps of a single parameter. 

To run the mono RHS ODEs with a sweep of enzyme a transcription rate (parameter w_ea) with a range 0-100, increasing by 10 each step:
```
python3 master_file.py -run parameter_sweep_1D -t mono -p1n w_ea -p1s 0:100:10
```


***2D_parameter_sweep*** : run sweeps of 2 parameters

To run the mono RHS ODEs, on glucose, with a sweep of enzyme a transcription rate across a range 0-100 and enzyme b transcription rate across a range 0-200, with steps of 20 for both
```
python3 master_file.py -run parameter_sweep_2D -t mono -s glucose -p1n w_ea -p1s 0:100:20 -p2n w_eb -p2s 0:200:20
```

***Options:*** Both sweep files will also accept the steady state and substrate arguments described in the "Run a single simulation to steady state" section above.

To run the duo RHS ODEs, on starch, from the glucose steady state, with a sweep of enzyme a transcription rate across a range 0-100 and enzyme b transcription rate across a range 0-200, with steps of 20 for both
```
python3 master_file.py -run parameter_sweep_2D -t duo -s starch -ss w -p1n w_ea -p1s 0:100:20 -p2n w_eb -p2s 0:200:20
```


## Read an existing pickle file
Outputs of simulations are stored as pickle files. To read the key infomation about a pre-run simulation and plot the results use the -read argument:
```
python3 master_file.py -read INSERT-FILENAME-HERE.pickle
```
Ensure the file you want to read is in your current working directory.


## Reproducing Figures

Figure 3a & 3b
```
python3 master_file.py -run run_ODE -t both -s glucose
```

Figure 3c
```
python3 master_file.py -run parameter_sweep_2D -t mono -s glucose -p1n w_ea -p1s 0:50:1 -p2n w_eb -p2s 0:50:1
```
```
python3 master_file.py -run parameter_sweep_2D -t duo -s glucose -p1n w_ea -p1s 0:50:1 -p2n w_eb -p2s 0:50:1
```

Figure 5a
```
python3 master_file.py -run parameter_sweep_2D -t mono -s starch -ss w -p1n w_ea -p1s 0:100:1 -p2n w_eb -p2s 0:100:1
```

Figure 5b
```
python3 master_file.py -run parameter_sweep_2D -t duo -s starch -ss w -p1n w_ea -p1s 0:100:1 -p2n w_eb -p2s 0:100:1
```

Figure 5c
Move the pickle files from 5a and 5b into the same folder as the compare.py file and then run:
```
python3 compare.py
```

