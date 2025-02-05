Cross-product Penalized Component Analysis (XCAN) v0.1

coded by: Jose Camacho Paez (josecamacho@ugr.es)
          Evrim Acar Ataman (evrim.acarataman@gmail.com)

last modification: 2/Feb/25


Please, see the following reference for a deep explanation on the method:

- J. Camacho, E. Acar, M. Rasmussen, R. Bro. Cross-product Penalized Component Analysis (XCAN), Submitted to Chemometrics and Intelligent Laboratory Systems, 2019.


The contents of the folder are:

- Examples: scripts with examples in the paper.

- xcan.m: main XCAN routine

- crossprod.m: routine to compute cross-product matrices

- rest of funtions are performing low-level computation.

The XCAN repo needs the Poblano Toolbox, so please remember to pull with the option: --recursive

git clone --recursive https://github.com/josecamachop/X-CAN.git

Or simply download the Poblano repository, version 1.2 (from https://github.com/sandialabs/poblano_toolbox/releases/tag/v1.2), in the subfolder poblano_toolbox.

To start, write in the Matlab command line:

# help xcan 
