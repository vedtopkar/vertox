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
  
  
## TODOs:  
- Add SLAP-seq enrichment data analysis and visualization