using SymCOOSolverInterface

using HSL, LinearAlgebra, Random, SparseArrays, Test

function tests()
  if isdefined(HSL, :libhsl_ma57)
    test_solver(MA57Struct)
  end
  test_solver(LDLFactorizationStruct)
end

function test_solver(LinearSolver)
  Random.seed!(0)
  nvar, ncon = 5, 3
  Q = spdiagm(nvar, nvar, 0 => 2 * ones(nvar), -1 => -ones(nvar-1), 1 => -ones(nvar-1))
  Qt = tril(Q)
  A = spdiagm(ncon, nvar, (i => ones(ncon) for i = (0:nvar-ncon))...)

  @testset "Positive definite systems" begin
    rows, cols, vals = findnz(Qt)
    M = LinearSolver(nvar, rows, cols, vals)
    sol = rand(nvar)
    rhs = Q * sol
    x = zeros(nvar)

    factorize!(M)
    solve!(x, M, rhs)
    @test norm(x - sol) ≤ 1e-8 * norm(sol)
    @test success(M)
    @test isposdef(M)
  end

  @testset "Indefinite systems" begin
    λ = sort(eigen(Matrix(Q)).values)
    α = (λ[1] + λ[2]) / 2
    rows, cols, vals = findnz(Qt - α * I)
    M = LinearSolver(nvar, rows, cols, vals)
    sol = rand(nvar)
    rhs = Q * sol - α * sol
    x = zeros(nvar)

    factorize!(M)
    solve!(x, M, rhs)
    @test norm(x - sol) ≤ 1e-8 * norm(sol)
    @test success(M)
    @test !isposdef(M)
    @test num_neg_eig(M) == 1
  end

  @testset "Failed factorization" begin
    rows = [1, 5, 3]
    cols = [1, 3, 5]
    vals = ones(3)
    M = LinearSolver(nvar, rows, cols, vals)

    factorize!(M)
    @test !success(M)
    @test !isposdef(M)
  end

  @testset "SQD system" begin
    rows, cols, vals = findnz([[Qt; A] zeros(nvar + ncon, ncon)])
    M = LinearSolver(nvar + ncon, rows, cols, vals)
    sol = rand(nvar + ncon)
    rhs = [Q A'; A zeros(ncon, ncon)] * sol
    x = zeros(nvar + ncon)

    factorize!(M)
    solve!(x, M, rhs)
    @test norm(x - sol) ≤ 1e-8 * norm(sol)
    @test success(M)
    @test !isposdef(M)
    @test num_neg_eig(M) == ncon
  end
end

tests()