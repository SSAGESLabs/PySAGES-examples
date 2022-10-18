
# PySAGES-examples

Examples are divided in two main categories: `classic examples` and `research examples`.

The repo is organized as follows:

```
PySAGES-examples
│   README.md
│
└───backend01
│   │
│   └───classic
│   │   │
│   │   └───system01
│   │   │   │
│   │   │   └───method01
│   │   │   │   │   system01_method01.ipynb
│   │   │   │   │   system01_method01.md
│   │   │   │   │   system01_method01.py
│   │   │   │
│   │   │   └───method02
│   │   │   │   │   system01_method02.ipynb
│   │   │   │   │   system01_method02.md
│   │   │   │   │   system01_method02.py
│   │   │   │
│   │   │   └───method03
│   │   │       │   ...
│   │   │
│   │   └───system02
│   │       │   ...
│   │
│   └───research
│       │
│       └───system01
│       │   │
│       │   └───method01
│       │   │   │   system01_method01.ipynb
│       │   │   │   system01_method01.md
│       │   │   │   system01_method01.py
│       │   │
│       │   └───method02
│       │       │   ...
│       │  
│       └───system02
│           │   ...
│   
└───backend02
│   │   ...
...
```
where `backend?? = ['ase', 'hoomd-blue', 'openmm']`; 
`system?? = ['butane', 'alaninedp', 'nacl', etc.]`;
`method?? = ['abf', 'cff', etc.]`.


Moreover, the folder [inputs/](./inputs/) contains some auxiliary files
such as starting coordinates.


## List of available examples

| system     | method       | CV                 | script                                             | notebook |
|------------|--------------|--------------------|----------------------------------------------------|----------| 
| butane     | ANN          | dihedral angle     | [script.py](./hoomd-blue/classic/butane/ann/butane_ANN.py) | [![butane_ANN](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/hoomd-blue/classic/butane/ann/butane_ANN.ipynb) |
| butane     | FUNN         | dihedral angle     | [script.py](./hoomd-blue/classic/butane/funn/butane_FUNN.py) | [![butane_FUNN](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/hoomd-blue/classic/butane/funn/butane_FUNN.ipynb) |
| butane     | CFF          | dihedral angle     | [script.py](./hoomd-blue/classic/butane/cff/butane_CFF.py) | [![butane_CFF](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/hoomd-blue/classic/butane/cff/butane_CFF.ipynb) |
| butane     | SpectralABF  | dihedral angle     | [script.py](./hoomd-blue/classic/butane/spectral_abf/butane_SpectralABF.py) | [![butane_SpectralABF](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/hoomd-blue/classic/butane/spectral_abf/butane_SpectralABF.ipynb) |
| alanine dp | Metadynamics | 2 x dihedral angle | [script.py](./openmm/classic/alaninedipeptide/metad/adp_Metadynamics.py) | [![adp_Metadynamics](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/openmm/classic/alaninedipeptide/metad/adp_Metadynamics.ipynb) |
| alanine dp | SpectralABF  | 2 x dihedral angle | [script.py](./openmm/classic/alaninedipeptide/spectral_abf/adp_SpectralABF.py) | [![adp_SpectralABF](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/openmm/classic/alaninedipeptide/spectral_abf/adp_SpectralABF.ipynb) |

*ANN*  = Artificial Neural Network sampling

*FUNN* = adpative Force-biasing sampling Using Neural Networks, or FUNN-ABF [[JCP **2018**, 148, 134108]](https://doi.org/10.1063/1.5020733)

*CFF*  = Combined Force-Frequency sampling [[JCTC **2020**, 16, 1448−1455]](https://doi.org/10.1021/acs.jctc.9b00883)

*SpectralABF* = Spectral Adaptive Biasing Force [[arXiv:2202.01876, **2022**]](https://arxiv.org/abs/2202.01876) 
