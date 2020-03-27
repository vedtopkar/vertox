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
In the future, I plan to distribute via pip. For now, there are a few manual steps:

1. `git clone` this repo
2. Make sure you have dependencies installed: `pip3 install -r requirements.txt`
3. Set up the Linearfold RNA folding engine.
    1. Clone and compile the [linearfold repo](https://github.com/LinearFold/LinearFold) in the location of your choice.
    2. Copy the `path.example.py` file from this repo to a new file called `path.py` in the `vertox` subdirectory. (i.e. from the top of this repo run `cp path.example.py vertox/path.`)
    3. Specify the location of the `linearfold` executable by changing the value of `PATH_LINEARFOLD_EXECUTABLE` in your new `path.py` file.
    

## Usage
Generation of a mutate-rescue library can be done by setting up a `.cfg` python configuration file with the fields as indicated in the example below. Note that for now, these are all required fields!

### Data Input
- `master_sequence`: Sequence of the overall region being scanned
- `master_structure`: Secondary structure of the master sequence
- `regions`: A list of indices specifying which specific sub-sequences of the `master_sequence` should be used for mutate-rescue design. Each set of indices should be on its own line and indented, and the start:stop indices should be separated by a colon (`:`). These indices are ZERO INDEXED!

## Library Parameters
- `five_prime_adapter`: Sequence to prepend to each final designed sequence (e.g. a cloning or priming overhang)
- `three_prime_adapter`: Sequence to append to each final designed sequence (e.g. a cloning or priming overhang)
- `total_iterations`: Number of random structures to try during stochastic stem disruption and rescue design.
- `min_helix_size`: Minimum size of helix (number of consecutive base pairs) that should be considered for mutate-rescue design.
- `n_stochastic_results`: How many of the highest ranking designs to keep for each stochastic design run. (1 means that only the top design is used).

## IO Parameters
- `experiment_name`: Name of experiment
- `output_folder`: Path to folder to put output files

```
[Data Input]
master_sequence=TTCTAATACGACTCACTATAGGCCAAAGGCGTCGAGTAGACGCCAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCGTCAGCGAGTAGCTGACAAAAAGAAACAACAACAACAAC
master_structure=...........................((((((.....))))))...........((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...)).............((((((.....))))))......................

# regions are ZERO-INDEXED and are inclusive of the first and last index
regions= 27:43
         216:240

[Library Parameters]
five_prime_adapter=GGGGCCCCTTTTAAAA
three_prime_adapter=AAAATTTTCCCGGGG
total_iterations=100
min_helix_size=3
n_stochastic_results=1

[IO Parameters]
experiment_name=P4P6
output_folder=test/output

```

To run the program:
```
python3 /path/to/vertox.py /path/to/config.cfg
```

The program outputs two `.xlsx` files. One with the sequence, structure, and design metrics (divided into one sheet per region indicated). The other simply has all of the designed sequences plus the appropriate 5' and 3' adapters. This spreadsheet is designed to be easily inputted into the IDT oPool upload form for quick ordering.

  
## TODOs:  
- Add SLAP-seq enrichment data analysis and visualization
- Add more comprehensive/detailed logging
