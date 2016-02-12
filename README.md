HEigs
=====

Haskell interface to ARPACK for large sparse eigenvalue problems.

ARPACK is a Fortran code for computing a few eigenpairs associated with large
sparse linear systems. This package wraps a subset of ARPACK's functionality and
attempts to deliver something similar to scipy or MATLAB's eigs functions.

To solve an eigensystem `A x = lambda x` the user needs to define an
`ArpackLinearOp`


    type ArpackLinearOp = (SV.IOVector CDouble -> SV.IOVector CDouble -> IO ())


This operator should overwrite the second vector with the matrix times the first
vector. To compute the eigenvalues call

    eigs :: ArpackLinearOp -> ProblemDim -> Which -> NumEV -> Tolerance -> 
            MaxIter -> IO (Bool, [(Complex Double, V.Vector (Complex Double))])
        
Where the parameters are:

- `type ProblemDim = Int` : The size of the linear system.
- `data Which = LM | SM` : Which eigenvalues to compute, those with largest
  magnitude (`LM`) or smallest magnitude (`SM`)
- `type NumEV = Int` : The number of eigenpairs to compute.
- `type Tolerance = Double` : The error tolerance. Setting to 0 uses the machine
  eps.
- `type MaxIter = Int` : The maximum number of Arnoldi iterations.

This function returns `(False, [])` if ARPACK was unsuccessful. Otherwise the
second element in the tuple contains a list of the computed NumEV eigenvalues
and eigenvectors.

Notes
======

Currently building with:

- ARPACK: http://forge.scilab.org/index.php/p/arpack-ng/source/tree/3.1.5/
- OpenBlas:
  https://github.com/xianyi/OpenBLAS/tree/7b8604ea29f7a7c258b3d7faf1f81eef750280f6


History
=======

0.0.1 First working interface to ARPACK. Can solve Ax = lambda x for largest and
      smallest magnitude eigenproblems.
      
Future Plans
============

- Support additional values for Which.
- Write a second wrapper for ARPACK's symmetric driver.

I usually only write what I need for today. If you need to access other portions of ARPACKs functionality
please let me know and I'll do my best to support it.

Chris Miller
cmiller730@gmail.com
