# Paper title

This is the official repository of

[**Title**](url) **(Published in Nature Materials, 2025)**

Authors....

<p align="center">
  <img src="/figures/TOC.png" width="600"/> 
</p>

<p align="center">
  <a href="/LICENSE">
      <img alt="license" src="https://img.shields.io/badge/license-MIT-green" />
  </a>

  <a href="https://colab.research.google.com/drive/1pIdzgTtcXFj7JGqIAyhQLl41j4ksE11E?usp=sharing">
      <img alt="demo" src="https://colab.research.google.com/assets/colab-badge.svg" />
  </a>
</p>



**ZeoSyn** is a large **zeolite synthesis dataset** comprising **23,961 zeolite synthesis routes**, 233 zeolite topologies and 921 organic structure-directing agents (OSDAs).
Each unique synthesis route consists of a **comprehensive set of key synthesis parameters**:
1. Gel compositions (molar ratios between heteroatoms, mineralizing agents, and water)
2. Reaction conditions (crystallization/aging temperature and time)
3. Organic structure-directing agent (SMILES)
4. Resultant zeolite product (3-letter IZA code)

## 1) Quick demo

We highly encourage you to check out our **[Demo notebook](/demo.ipynb)** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1pIdzgTtcXFj7JGqIAyhQLl41j4ksE11E?usp=sharing) for a gentle introduction (**< 3 min** 🎉) on the key components of dataset + SHAP for frameworks and building units.

**⚠️Note:** We strongly recommend the **Chrome** browser for the Google Colab notebooks


## 2) The ZeoSyn dataset

### A) Overview
![Alt text](/figures/figure1.jpg "overview")
**(a)** Example of a zeolite synthesis route in the
dataset, consisting of the gel composition, inorganic precursors, reaction conditions, organic
structure-directing agent (OSDA), and the resultant zeolite framework. Paper metadata of
the scientific paper containing the synthesis route is also provided. **(b)** Frequency of elements
present in the dataset. The values correspond to the log number of synthetic routes with a
specific element. **(c)** Total number of synthesis routes of small, medium, large, and extra-large pore zeolites extracted from literature across time in the dataset. Distributions of key gel composition variables in the dataset, including ratio between **(d)** heteroatoms, and **(e)**
mineralizing agents, metal cations and OSDA ratios (T = ∑i ni where ni is the amount of the ith heteroatom present in synthesis).

### B) Zeolite frameworks
![Alt text](/figures/zeo_distribution_by_zeotype_pore.png "frameworks")
Zeolite frameworks can be divided into different categories based on their maximum ring
size. ZeoSyn contains 5250, 5494, 5769, and 716 synthesis routes for small (8MR), medium
(10MR), large (12MR), and extra-large pore (>12MR) zeolites, respectively.

### C) Organic structure-directing agents
![Alt text](/figures/osda_hierarchy.jpg "osda")
**(a)** Hierarchical clustering of the top 50 most frequent OSDAs in the dataset,
labled with the main classes of molecular structures. Splits are obtained through agglomer-
ative hierarchical clustering of OSDA Morgan fingerprints. Each OSDA is colored by its
molecular volume (orange), and median largest included sphere of zeolites formed by the
OSDA (purple). The concomitant intensities of the colors show a positive correlation between the two properties. **(b)** Positive correlation between zeolite largest included sphere vs.
OSDA volume. Red points refer to high asphericity, which account for outliers **(c)** Positive
correlation between zeolite ring size vs. OSDA volume.

### D) SHAP analysis reveals key zeolite structure-synthesis relationships
![Alt text](/figures/figure4.jpg "shap")
**(a)** Framework-level SHAP analysis revealing the top 10 (out of 43) most important
synthesis parameters favoring the formation of specific frameworks. Each framework belongs
to 1 out of 3 types of synthesis based on its top synthesis parameters: 1) Gel-dominated synthesis (CAN, KFI) where most top parameters are inorganic-related, 2) OSDA-dominated
synthesis (ISV, ITE) where most top parameters are OSDA-related, and 3) balanced syn-
thesis (IWW, RUT) where even attribution is given to inorganic and OSDA parameters.
Every point is an individual synthesis colored by the value of synthesis parameter (orange
and blue colors indicate high and low values, respectively). **(b)** CBU-level SHAP analysis
of large CBUs showing OSDA parameters favoring their formation.

## 3) Setup and installation

The code in this repo has been tested on a Linux machine running Python 3.8.8

Run the following terminal commands 

1. Clone repo to local directory

```bash
  git clone https://github.com/eltonpan/zeosyn_dataset.git
```

2. Set up and activate conda environment
```bash
  cd zeosyn_dataset
```
```bash
  conda env create -f env/env.yml
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

## 4) Code reproducibility

All data required to reproduce results in the paper can be found in `datasets/` directory:

* [`ZEOSYN.xlsx`](/dataset/ZEOSYN.xlsx): ZeoSyn dataset
* [`osda_descriptors.csv`](/dataset/osda_descriptors.csv): Descriptors of organic structure-directing agents
* [`zeolite_descriptors.csv`](/dataset/zeolite_descriptors.csv): Descriptors of zeolite frameworks

All visualizations, model training and SHAP analysis in the paper can be reproduced by running the code in the following: 
* [visualization.ipynb](/visualization.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1yUN0HfLVGgvQfThnRSFJjWw1Xwyncg0s?usp=sharing) (~10 min, in-depth visualization of dataset)
* [classifier.ipynb](/classifier.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1vKvWxqcP0Cs4CxCqDcL1TQwp8EJYw1PZ?usp=sharing) (~15 min, deeper dive into zeolite classifier + SHAP for frameworks, building units, competing phases and intergrowths)

**If you are not using Colab:** Computation of SHAP values takes a while (~2 hours to run on 32 CPU cores). To avoid computation of SHAP, you can choose to download and load the precomputed SHAP values:
1. Download `shap_values.pkl` from [here](https://figshare.com/s/5519f7668ff2f631f47f)
2. Place `shap_values.pkl` in `shap/` directory
3. Make sure the following block in `classifier.ipynb` is uncommented
  ```python 
    with open('shap/shap_values.pkl', 'rb') as handle:
        shap_values = pickle.load(handle)
  ``` 


## 5) Cite
If you use this dataset or code, please cite this paper:
```
```

## 6) Contact
If you have any questions, please contact us at [elsao@mit.edu](mailto:elsao@mit.edu).

## To-do:
- [ ] Test conda installation on non-Linux systems
- [x] Check OSDA and zeolite descriptors have any redundant data
- [x] Add Bibtex
- [X] Update Bibtex (after issue is out)
- [x] Add Colab notebook option
- [x] Add contact info