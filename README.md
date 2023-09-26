# Code for *ZeoSyn: A Comprehensive Zeolite Synthesis Dataset Enabling Machine-learning Rationalization of Hydrothermal Parameters*

Elton Pan,† Soonhyoung Kwon,‡ Zach Jensen,† Manuel Moliner,¶ Yuriy Roman,‡ and Elsa Olivetti∗,†

† Department of Materials Science and Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts 02139, United States

‡ Department of Chemical Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts 02139, United States

¶ Instituto de Tecnolog ́ıa Qu ́ımica, Universitat Politecnica de Valencia-Consejo Superior de
Investigaciones Cient ́ıficas, 46022 Valencia, Spain

**ZeoSyn** is a comprehensive dataset of **23,925 unique zeolite hydrothermal synthesis routes**, encom-
passing over 200 zeolite topologies and 900 organic structure-directing agents (OSDAs).
Each unique synthesis route consists of a **comprehensive set of key synthesis variables**: gel compositions (molar ratios between heteroatoms, mineralizing agents, and water), reaction conditions (crystallization temperature and time), organic structure-directing agent, and resultant zeolite product.

![Alt text](/figures/overview.png "overview")

![Alt text](/figures/osda_hierarchy.png "osda")

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