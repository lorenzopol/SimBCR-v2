# SimBCR
Python script for simulating recombinant receptor 3D structure in a fully autonomous fashion.

> This is still an experimental build, once this project will be over a full documentation and contribution pipeline will be released.

## Current capability
1. Simulate variable regions (check literature \[1]) and join them with constant regions
2. Dynamically find CDR3 region
3. Autonomous protein folding (check literature \[2])
4. Autonomous pdb rendering via py3dmol and jupyter-notebook
5. Two different ways (rainbow and CDR) to color your antibody once folding is done 

## Installation
> WARNING: these are the minimum steps needed to run the CURRENT default behaviour of SimBCR. No other function is being supported 

1. In order to get all the necessary python packages run
`pip install -r requirements.txt`

2. You will need to install R and the immuneSIM package. Please refer to their [installation section](https://github.com/GreiffLab/immuneSIM).
3. Finally, clone this repo or download and extract the raw files. 

## Usage
Run main.py (from the root directory of the project) with 

`python main.py`

If everything runs smoothly, a web page should open with your folded heavy chain (only the FAB portion) colored by CDR (red is CDR1, green is CDR2 and blue is CDR3).
Please keep in mind that the default settings simulate a human IgG heavy chain.

## Develop timeline from 04/12/23 on:
- [x] Simulate variable regions with this [R project](https://github.com/GreiffLab/immuneSIM)
- [x] AlphaFold. Previously using [SWISS-MODEL](https://swissmodel.expasy.org/interactive) in order to get the 3d structure from the aminoacid sequence. Switched to ESMatlas API (even if it does not support more than 400AA prediction. Consider hybrid approach)
- [ ] CDR1 and CDR2 are still find in a hard coded, literature-based method. Is good enough? 
- [ ] Dreaming of a custom OpenGL-based 3D renderer for custom display of folded protein. Currently using py3dmol with hacky jupiternotebook trick (check show_3d.py and render_pdb.ipynb)


# Literature and resource references
1. https://academic.oup.com/bioinformatics/article/36/11/3594/5802461
2. https://esmatlas.com/about#api