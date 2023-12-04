# SimBCR
Python script for simulating recombinant receptor 3D structure in a fully autonomous fashion.

> This is still an experimental build, once this project will be over a full documentation and contribution pipeline will be released.

## Current capability
1. Simulate variable regions (check literature \[1]) and join them with constant regions
2. Dynamically find CDR3 region
3. Autonomous protein folding (check literature \[2])
4. Autonomous pdb rendering via py3dmol and jupyter-notebook
5. Two different ways (rainbow and CDR) to color your antibody once folding is done 

## Develop timeline from 04/12/23 on:
- [x] Simulate variable regions with this [R project](https://github.com/GreiffLab/immuneSIM)
- [x] AlphaFold. Previously using [SWISS-MODEL](https://swissmodel.expasy.org/interactive) in order to get the 3d structure from the aminoacid sequence. Switched to ESMatlas API (even if it does not support more than 400AA prediction. Consider hybrid approach)
- [ ] CDR1 and CDR2 are still find in a hard coded, literature-based method. Is good enough? 
- [ ] Dreaming of a custom OpenGL-based 3D renderer for custom display of folded protein. Currently using py3dmol with hacky jupiternotebook trick (check show_3d.py and render_pdb.ipynb)


# Literature and resource references
1. https://academic.oup.com/bioinformatics/article/36/11/3594/5802461
2. https://esmatlas.com/about#api