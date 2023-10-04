# Code for *ZeoSyn: A Comprehensive Zeolite Synthesis Dataset Enabling Machine-learning Rationalization of Hydrothermal Parameters*

Elton Pan,† Soonhyoung Kwon,‡ Zach Jensen,† Manuel Moliner,¶ Yuriy Roman,‡ and Elsa Olivetti∗,†

† Department of Materials Science and Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts 02139, United States

‡ Department of Chemical Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts 02139, United States

¶ Instituto de Tecnolog ́ıa Qu ́ımica, Universitat Politecnica de Valencia-Consejo Superior de
Investigaciones Cient ́ıficas, 46022 Valencia, Spain

### Quick links:
1. [Setup and installation](#setup-and-installation)
2. [Demo notebook](/demo.ipynb)
3. [ZeoSyn dataset](/dataset/ZEOSYN.xlsx)

### To-do:
- [ ] Test conda installation on Mac and Windows
- [ ] Check OSDA and zeolite descriptors have any redundant/data
- [ ] Add Bibtex

**[ZeoSyn](/dataset/ZEOSYN.xlsx)** is a dataset of **23,925 unique zeolite hydrothermal synthesis routes**, encompassing over 200 zeolite topologies and 900 organic structure-directing agents (OSDAs).
Each unique synthesis route consists of a **comprehensive set of key synthesis variables**:
1. [Gel compositions](#overview-of-zeosyn-dataset) (molar ratios between heteroatoms, mineralizing agents, and water)
2. [Reaction conditions](#overview-of-zeosyn-dataset) (crystallization temperature and time)
3. [Organic structure-directing agent](#common-organic-structure-directing-agents-in-the-zeosyn-dataset) (SMILES)
4. [Resultant zeolite product](#common-zeolite-frameworks-in-the-zeosyn-dataset-by-pores-size) (3-letter IZA code)

### Overview of ZeoSyn dataset
![Alt text](/figures/overview.png "overview")
**(a)** Example of a zeolite synthesis route in the
dataset, consisting of the gel composition, inorganic precursors, reaction conditions, organic
structure-directing agent (OSDA), and the resultant zeolite framework. Paper metadata of
the scientific paper containing the synthesis route is also provided. **(b)** Frequency of elements
present in the dataset. The values correspond to the log number of synthetic routes with a
specific element. **(c)** Total number of synthesis routes of small, medium, large, and extra-large pore zeolites extracted from literature across time in the dataset. Distributions of key
gel composition variables in the dataset, including ratio between **(d)** heteroatoms, and **(e)**
mineralizing agents, metal cations and OSDA ratios (T = ∑
i ni where ni is the amount of
the ith heteroatom present in synthesis).

### Common zeolite frameworks in the ZeoSyn dataset (by pores size)
![Alt text](/figures/zeo_distribution_by_zeotype_pore.png "frameworks")
Zeolite frameworks can be divided into different categories based on their maximum ring
size. ZeoSyn contains 5250, 5494, 5769, and 716 synthesis routes for small (8MR), medium
(10MR), large (12MR), and extra-large pore (>12MR) zeolites, respectively.

### Common organic structure-directing agents in the ZeoSyn dataset
![Alt text](/figures/osda_hierarchy.png "osda")
**(a)** Hierarchical clustering of the top 50 most frequent OSDAs in the dataset,
labled with the main classes of molecular structures. Splits are obtained through agglomer-
ative hierarchical clustering of OSDA Morgan fingerprints. Each OSDA is colored by its
molecular volume (orange), and median largest included sphere of zeolites formed by the
OSDA (purple). The concomitant intensities of the colors show a positive correlation between the two properties. **(b)** Positive correlation between zeolite largest included sphere vs.
OSDA volume. Red points refer to high asphericity, which account for outliers **(c)** Positive
correlation between zeolite ring size vs. OSDA volume.

### SHAP analysis revealing the most important synthesis parameters favoring the formation of specific zeolite frameworks and composite building units
![Alt text](/figures/SHAP_zeolite_cbu.png "shap")
**(a)** Framework-level SHAP analysis revealing the top 10 (out of 43) most important
synthesis parameters favoring the formation of specific frameworks. Each framework belongs
to 1 out of 3 types of synthesis based on its top synthesis parameters: 1) Gel-dominated synthesis (CAN, KFI) where most top parameters are inorganic-related, 2) OSDA-dominated
synthesis (ISV, ITE) where most top parameters are OSDA-related, and 3) balanced syn-
thesis (IWW, RUT) where even attribution is given to inorganic and OSDA parameters.
Every point is an individual synthesis colored by the value of synthesis parameter (orange
and blue colors indicate high and low values, respectively). **(b)** CBU-level SHAP analysis
of large CBUs showing OSDA parameters favoring their formation.

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

## Getting started
A demo of loading and accessing key components of the ZeoSyn data can be found in `demo.ipynb`

## Reproducibility
The data can be found in `/datasets` directory:

* **[`ZEOSYN.xlsx`](/dataset/ZEOSYN.xlsx): ZeoSyn dataset**
* [`osda_descriptors.csv`](/dataset/osda_descriptors.csv): Descriptors of organic structure-directing agents
* [`zeolite_descriptors.csv`](/dataset/zeolite_descriptors.csv): Descriptors of zeolite frameworks

All visualizations, model training and SHAP analysis in the paper can be reproduced by running the code in the following: 

* [`visualization.ipynb`](/visualization.ipynbvisu): For visualizations of the ZeoSyn dataset
* [`classifier.ipynb`](/classifier.ipynb): SHAP of zeolite phase predictor model

(OPTIONAL) Computation of SHAP values takes a while (~2 hours to run on 32 CPU cores). To avoid computation of SHAP, you can choose to download and load the precomputed SHAP values:
1. Download `shap_values.pkl` from [here](https://figshare.com/s/5519f7668ff2f631f47f)
2. Place `shap_values.pkl` in `shap/` directory


## Cite
If you use this dataset or code, please cite this paper:
```
<INSERT BIBTEX HERE>
```