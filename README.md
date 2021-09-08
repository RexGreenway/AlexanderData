Alexander Data
==============

A library for the calculation and understanding of the Alexander Data.

Created in progress toward completion of COMP702 - Dissertation
Module, MA Data Science and Artificial Intelligence @ University
of Liverpool.

This library was written for the purpose of calculating the 
Alexander Data, an invariant for textiles. To achieve this goal
this package allows for the creation of Braid objects - 
specifically the Braid Kernel, a representation of repeating
textile structures. From this object the user can calculate the
entire Alexander Data, or each of the internel stages piece-meal.

Features
--------
- Braid, and Braid Kernel, class object representation.
- Reduced Burau, Alexander Data, and Alexander Data calculation.
- Briad and Kernel visulaisation functionality.

Installation
------------
This package is dependent on some other common Python data science
libraries:

- Numpy
- Sympy
- Matplotlib (For Visualisation Purposes)

One reccomended way of running the Alexander Data module is to run
it using Anaconda, or in a conda environment, as it provides all
above requirements.

Once the dependencies are accounted for you can install the library
very simply using Pip:
```
pip install git+https://github.com/RexGreenway/polylatlib
```

