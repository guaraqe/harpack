HArpack
=======

This is a fork of [HEigs](https://github.com/cmiller730/HEigs/) with
similar functionality but a fairly different interface.

It currently supports symmetrical and non-symmetrical eigenvalues
problems of the form `Mx = λx`

Its fundamental data type is the `SparseMatrix`:

    data SparseMatrix =
      SparseMatrix { dim ::Int
                   , indexes :: [((Int,Int),Double)] }

The eigenpair problem is represented by the `Problem` type:

    data Problem = Symmetric SWhich
                 | NonSymmetric NWhich

where `SWhich` and `NWhich` represents the eigenvalues to be found.
  - The first character refers to 'S'ymmetric or 'N'on-Symmetric;
  - The second refers to 'L'argest or 'S'mallest;
  - The third refers to 'M'agnitude, 'A'lgebraic, 'R'eal, 'I'maginary;
  - The option 'BE' means 'BE'tween, for central eigenvalues.

The data types are:

    data SWhich = SLM | SSM | SLA | SSA | SBE
    data NWhich = NLM | NSM | NLR | NSR | NLI | NSI
    
The result is given as:

    data Result = RError
                | RReal [(Double,Vector Double)]
                | RComplex [(Complex Double, Vector (Complex Double))]

This decision can certainly be improved, I accept suggestions. Finally, the `eigs` function does what we expect:

    eigs :: SparseMatrix
         -> Problem
         -> HowMany
         -> Result


Notes
=====

Currently building with:

- ARPACK: http://forge.scilab.org/index.php/p/arpack-ng/source/tree/3.1.5/
- OpenBlas:
  https://github.com/xianyi/OpenBLAS/tree/7b8604ea29f7a7c258b3d7faf1f81eef750280f6
