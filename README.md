# Code for *ZeoSyn: A Comprehensive Zeolite Synthesis Dataset Enabling Machine-learning Rationalization of Hydrothermal Parameters*

Elton Pan,† Soonhyoung Kwon,‡ Zach Jensen,† Manuel Moliner,¶ Yuriy Roman,‡ and Elsa Olivetti∗,†

† Department of Materials Science and Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts 02139, United States

‡ Department of Chemical Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts 02139, United States

¶ Instituto de Tecnolog ́ıa Qu ́ımica, Universitat Politecnica de Valencia-Consejo Superior de
Investigaciones Cient ́ıficas, 46022 Valencia, Spain

Zeolites, crystalline aluminosilicate materials with well-defined porous structures,
have emerged as versatile materials with applications in various fields, including catal-
ysis, gas separation, and ion exchange. Hydrothermal synthesis is the most widely used
method for zeolite production, offering control over composition, crystallinity, and pore
size. However, the intricate interplay of synthesis parameters and their impact on ze-
olite properties necessitates a comprehensive understanding to optimize the synthesis
process. Hitherto, publicly available zeolite synthesis databases only contain a subset
of key parameters. 

We present ZeoSyn, a comprehensive dataset of 23,925 unique zeolite hydrothermal synthesis routes, encom-
passing over 200 zeolite topologies and 900 organic structure-directing agents (OSDAs).
Each unique synthesis route consists of a comprehensive set of key synthesis variables:
1) gel compositions (molar ratios between heteroatoms, mineralizing agents, and wa-
ter) 2) reaction conditions (crystallization temperature and time), 3) OSDA, and 4)
resultant zeolite product. Extracted from over 50 years of zeolite scientific literature
using a semi-automatic approach leveraging a natural language processing pipeline, the
dataset enables a holistic analysis of zeolite synthesis parameters and their influence on
the final zeolite product. 

![Alt text](/figures/overview.png "overview")

![Alt text](/figures/osda_hierarchy.png "osda")

![Alt text](/figures/SHAP_zeolite_cbu.png "shap")