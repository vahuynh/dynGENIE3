# dynGENIE3
Semi-parametric approach for the inference of gene regulatory networks from time series of expression data

The dynGENIE3 method is described in the following paper (available [here](https://www.nature.com/articles/s41598-018-21715-0)):
```
dynGENIE3: dynamical GENIE3 for the inference of gene networks from time series expression data
Huynh-Thu, V. A. and Geurts, P.
Scientific Reports, 8:3384, 2018.
```

Three implementations of dynGENIE3 are available: Python, MATLAB and R. Each folder contains a PDF file with a step-by-step tutorial showing how to run the code.

Note: All the results presented in the paper were generated using the Python implementation.

dynGENIE3 is based on regression trees. To learn these trees, the Python implementation uses the [scikit-learn](http://scikit-learn.org/) library, and the MATLAB and R implementations are respectively MATLAB and R wrappers of a C code written by [Pierre Geurts](http://www.montefiore.ulg.ac.be/~geurts/). 