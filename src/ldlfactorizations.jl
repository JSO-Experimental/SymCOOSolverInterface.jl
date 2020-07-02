using LDLFactorizations, LinearAlgebra, SparseArrays

export LDLFactorizationStruct

# LDLFactorizations
mutable struct LDLFactorizationStruct <: SymCOOSolver
  ndim :: Int
  rows :: Vector{Int}
  cols :: Vector{Int}
  vals :: Vector
  factor
end

function LDLFactorizationStruct(ndim, rows, cols, vals)
  LDLFactorizationStruct(ndim, rows, cols, vals, nothing)
end

function factorize!(M :: LDLFactorizationStruct)
  try
    A = sparse(M.cols, M.rows, M.vals, M.ndim, M.ndim)
    M.factor = ldl(A, upper=true)
  catch ex
    M.factor = nothing
  end
end

function solve!(x, M :: LDLFactorizationStruct, b)
  ldiv!(x, M.factor, b)
end

function success(M :: LDLFactorizationStruct)
  M.factor != nothing
end

function isposdef(M :: LDLFactorizationStruct)
  success(M) && count(M.factor.d .≤ 1e-14) == 0
end

function num_neg_eig(M :: LDLFactorizationStruct)
  count(M.factor.d .≤ -1e-14)
end