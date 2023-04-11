
# PySAGES-examples

Repository containing an extensive set of examples for enhanced sampling simulations using [PySAGES](https://arxiv.org/abs/2301.04835).\
The repo containing the PySAGES code can be found [here](https://github.com/SSAGESLabs/PySAGES).


## Repo organization

Examples are divided in two main categories: `classic examples` and `research examples`.\
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
| butane     | Metadynamics | dihedral angle     | [script.py](./hoomd-blue/classic/butane/metad/butane_Metadynamics.py) | --- |
| butane     | SpectralABF  | dihedral angle     | [script.py](./hoomd-blue/classic/butane/spectral_abf/butane_SpectralABF.py) | [![butane_SpectralABF](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/hoomd-blue/classic/butane/spectral_abf/butane_SpectralABF.ipynb) |
| alanine dp | ABF          | 2 x dihedral angle | [script.py](./openmm/classic/alaninedipeptide/abf/adp_ABF.py) | --- |
| alanine dp | CFF          | 2 x dihedral angle | [script.py](./openmm/classic/alaninedipeptide/cff/adp_CFF.py) | --- |
| alanine dp | Metadynamics | 2 x dihedral angle | [script.py](./openmm/classic/alaninedipeptide/metad/adp_Metadynamics.py) | [![adp_Metadynamics](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/openmm/classic/alaninedipeptide/metad/adp_Metadynamics.ipynb) |
| alanine dp | SpectralABF  | 2 x dihedral angle | [script.py](./openmm/classic/alaninedipeptide/spectral_abf/adp_SpectralABF.py) | [![adp_SpectralABF](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/openmm/classic/alaninedipeptide/spectral_abf/adp_SpectralABF.ipynb) |
| alanine dp | ANN          | 2 x dihedral angle | [script.py](./openmm/classic/alaninedipeptide/ann/adp_ANN.py) | --- |
| alanine dp | unbiased     | 2 x dihedral angle | [script.py](./openmm/classic/alaninedipeptide/unbiased/adp_unbiased.py) | --- |
| NaCl       | Metadynamics | distance           | [script.py](./openmm/classic/NaCl/metad/nacl_Metadynamics.py) | [![nacl_Metadynamics](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SSAGESLabs/PySAGES-examples/blob/main/openmm/classic/NaCl/metad/nacl_Metadynamics.ipynb) |

*ANN*  = Artificial Neural Network sampling; 
*FUNN* = adpative Force-biasing sampling Using Neural Networks, or FUNN-ABF;
*CFF*  = Combined Force-Frequency sampling;
*SpectralABF* = Spectral Adaptive Biasing Force. 

See [the PySAGES manuscript](https://arxiv.org/abs/2301.04835) for more info on the methods and their citations.

