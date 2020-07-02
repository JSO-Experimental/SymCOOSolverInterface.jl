using HSL

export MA57Struct

if isdefined(HSL, :libhsl_ma57)
  mutable struct MA57Struct <: SymCOOSolver
    factor :: Ma57
  end

  function MA57Struct(N, rows, cols, vals)
    MA57Struct(ma57_coord(N, rows, cols, vals))
  end
else
  function MA57Struct(N, rows, cols, vals)
    error("MA57 not installed. See HSL.jl's README")
  end
end

function factorize!(M :: MA57Struct)
  ma57_factorize(M.factor)
end

function solve!(x, M :: MA57Struct, b)
  x .= b
  ma57_solve!(M.factor, x)
end

function success(M :: MA57Struct)
  M.factor.info.info[1] == 0
end

function isposdef(M :: MA57Struct)
  M.factor.info.num_negative_eigs == 0
end

function num_neg_eig(M :: MA57Struct)
  M.factor.info.num_negative_eigs
end