# SymCOOSolverInterface

This package provides a simple linear solver interface for factorization of symmetric matrices in COO format, where only the lower triangle is stored. It is made especially for optimization solvers using [NLPModels.jl](https://github.com/JuliaSmoothOptimizers/NLPModels.jl).

## Usage

Given lower triangular coordinates of a symmetric matrix in COO format `(cols, rows, vals)` one can create a structure ready for factorization using either (if installed) MA_57 from [`HSL.jl`](https://github.com/JuliaSmoothOptimizers/HSL.jl)
```
MA57Struct(n, rows, cols, vals)
```
or [`LDLFactorizations.jl`](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl)
```
LDLFactorizationStruct(n, rows, cols, vals)
```
Independently of the chosen structure, this package export the functions `factorize!`, `solve!`, `success`, `isposdef` and `num_neg_eig`.
