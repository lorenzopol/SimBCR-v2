# SimBCR
SimBCR is a python project for simulating recombinant receptor 3D structure in a fully autonomous fashion.

> This is still an experimental build, once this project will be over a full documentation and contribution pipeline will be released.

## Current capability
1. Simulate variable regions (check literature \[1]) and join them with constant regions
2. Dynamically find CDR3 region
3. Autonomous protein folding (check literature \[2])
4. Autonomous local pdb rendering via OpenGL-based (check literature \[3]) rendering
5. Free camera movement around your 3d rendered pdb file
6. Autonomous web based pdb rendering via py3dmol and jupyter-notebook
7. Two different ways (rainbow and CDR) to color your antibody once folding is done 

## Installation
> WARNING: these are the minimum steps needed to run the CURRENT default behaviour of SimBCR. No other functionality is being supported 

1. In order to get all the necessary python packages run
`pip install -r requirements.txt`

2. You will need to install R and the immuneSIM package. Please refer to their [installation section](https://github.com/GreiffLab/immuneSIM). Assure yourself of adding R to your system PATH variables.
3. Finally, clone this repo or download and extract the raw files. 

## Usage
### Default 
Run main.py (from the root directory of the project) with 

`python main.py`

If everything runs smoothly, a window should open with your folded heavy chain (only the FAB portion).
Please keep in mind that the default settings simulate a human IgG heavy chain.
### Custom inputs
Right now simBCR allows you to customize your input in two different ways: simple mode and developer mode.
#### Custom inputs: simple mode
In order to run a custom simulation in the easiest way possible, provide the following inputs when launching main.py

`python main.py  --number_of_seqs 100 --species hs --receptor tr --chain b --name_repertoire my_custom_sim` 

This will allow you to simulate a human TCR-B chain and save it as a "name_repertoire".csv file.
The default way of running SimBCR translates to:

`python main.py  --number_of_seqs 10 --species hs --receptor ig --chain h --name_repertoire hs-igh-sim --renderer local`
#### Custom inputs: developer mode
If you really want to dig deeper into customization you can edit the `InputParser` instance in `main.py` by adding any of the available keys present in `InputParser.default_main_R_args` (check `InputParser.py`) as a kwarg. For example, replacing the current creation of `InputParser` with `
InputParser(args, smh_mode = "naive", smh_prob = 100/350)` will cause the simulation of a variable region where each 350 nucleotides 100 will get mutated in a random way.
Keep in mind that the key (smh_mode here, no quotes!) need to be present in `InputParser.default_main_R_args` and the value ("naive" here) should be treated as a raw R input.
Using custom inputs in this way is NOT recommended and please refer to the [original documentation](https://immunesim.readthedocs.io/en/latest/parameters.html) if you need more information. 


## Develop timeline from 24/12/23 on:
- [x] Simulate variable regions with this [R project](https://github.com/GreiffLab/immuneSIM)
- [x] Previously using [SWISS-MODEL](https://swissmodel.expasy.org/interactive) in order to get the 3d structure from the aminoacid sequence. Switched to ESMatlas API (even if it does not support more than 400AA prediction. Consider hybrid approach)
- [x] Dreaming of a custom OpenGL-based 3D renderer for custom display of folded protein. Currently, I am using py3dmol with hacky jupiter-notebook trick (check show_3d.py and render_pdb.ipynb)
- [ ] AlphaFold support? 
- [ ] CDR1 and CDR2 are still find in a hard coded, literature-based method. Is good enough? 


# Literature and resource references
1. [ImmuneSIM article](https://academic.oup.com/bioinformatics/article/36/11/3594/5802461) and its [GitHub repo](https://github.com/GreiffLab/immuneSIM)
2. [ESMatlas API](https://esmatlas.com/about#api)
3. [Base Script for OpenGL renderer](https://github.com/StanislavPetrovV/3D-Graphics-Engine) from StanislavPetrovV