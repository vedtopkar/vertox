# VERTOX  
  
**Ve**d's **R**NA **To**olbo**x**  
  
This is a repo of scripts that I have written for manipulating primary sequence and secondary structure information.  
  
The initial consolidation of these scripts into this package is targeted towards generating libraries for mutate-rescue functional experiments.  
  
The long-term hope is to make this a general swiss-army-knife package that can easily be installed and used for RNA informatics.  


## Current Functionality
- A set of `secondary structure` functions that allow for quick parsing of dot-bracket strings into base-pairs and helices
- A `score_structure_perturbations` that scores variants of a WT sequence for the disruption of a certain helix and/or the recovery of the rest of the structure.
- A **parallelized** `stochastic_design` module that generates large numbers of helix variants for scoring.
- A `mutate_rescue` module that allows for easy and rapid design of mutate-rescue libraries
- A `linearfold` module allows for easy usage of that package withing Python


## Installation
In the future, I plan to distribute via pip. For now, there are three steps:

1. `git clone` this repo
2. Make sure you have dependencies installed: `pip3 install -r requirements.txt`
3. Set up the Linearfold RNA folding engine.
    1. Clone and compile the [linearfold repo](https://github.com/LinearFold/LinearFold) in the location of your choice.
    2. Copy the `path.example.py` file from this repo to a new file called `path.py` in the `vertox` subdirectory. (i.e. from the top of this repo run `cp path.example.py vertox/path.`)
    3. Specify the location of the `linearfold` executable by changing the value of `PATH_LINEARFOLD_EXECUTABLE` in your new `path.py` file.
    

## Usage
TODO 
  
  
## TODOs:  
- Add SLAP-seq enrichment data analysis and visualization
