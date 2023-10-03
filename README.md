# Code for *ZeoSyn: A Comprehensive Zeolite Synthesis Dataset Enabling Machine-learning Rationalization of Hydrothermal Parameters*

Elton Pan,† Soonhyoung Kwon,‡ Zach Jensen,† Manuel Moliner,¶ Yuriy Roman,‡ and Elsa Olivetti∗,†

† Department of Materials Science and Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts 02139, United States

‡ Department of Chemical Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts 02139, United States

¶ Instituto de Tecnolog ́ıa Qu ́ımica, Universitat Politecnica de Valencia-Consejo Superior de
Investigaciones Cient ́ıficas, 46022 Valencia, Spain

**ZeoSyn** is a dataset of **23,925 unique zeolite hydrothermal synthesis routes**, encompassing over 200 zeolite topologies and 900 organic structure-directing agents (OSDAs).
Each unique synthesis route consists of a **comprehensive set of key synthesis variables**:
1. Gel compositions (molar ratios between heteroatoms, mineralizing agents, and water)
2. Reaction conditions (crystallization temperature and time)
3. Organic structure-directing agent
4. Resultant zeolite product

### Overview of ZeoSyn dataset
![Alt text](/figures/overview.png "overview")

### Common zeolite frameworks in the ZeoSyn dataset (by pores size)
![Alt text](/figures/zeo_distribution_by_zeotype_pore.png "frameworks")

### Common organic structure-directing agents in the ZeoSyn dataset
![Alt text](/figures/osda_hierarchy.png "osda")

### SHAP analysis revealing the most important synthesis parameters favoring the formation of specific zeolite frameworks and composite building units
![Alt text](/figures/SHAP_zeolite_cbu.png "shap")

## Setup and installation

Run the following terminal commands 

1. Clone repo to local directory

```bash
  git clone https://github.com/eltonpan/ZeoSyn_dataset.git
```

2. Set up and activate conda environment
```bash
  cd ZeoSyn_dataset
```
```bash
  conda env create -f env.yml
```
```bash
  conda activate zeosyn
```

3. Add conda environment to Jupyter notebook
```bash
  conda install -c anaconda ipykernel
```
```bash
  python -m ipykernel install --user --name=zeosyn
```

4. Open jupyter notebooks
```bash
  jupyter notebook <notebook_name>.ipynb
```


make sure the `zeosyn` is the environment under dropdown menu `Kernel` > `Change kernel`

All visualizations, model training and SHAP analysis in the paper can be reproduced by running the code in the following: 

* `visualization.ipynb`: For visualizations of the ZeoSyn dataset
* `classifier.ipynb`: SHAP of zeolite phase predictor model


The data can be found in `/datasets` directory:

* `ZEOSYN.xlsx`: ZeoSyn dataset
* `osda_descriptors.csv`: Descriptors of organic structure-directing agents
* `zeolite_descriptors.csv`: Descriptors of zeolite frameworks

# Cite
If you use this dataset or code, please cite this paper:
```
<INSERT BIBTEX HERE>
```