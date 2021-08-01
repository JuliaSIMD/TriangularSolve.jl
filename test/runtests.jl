using TriangularSolve, LinearAlgebra
using Test

function test_solve(::Type{T}) where {T}
  for n ∈ 1:(T === Float32 ? 100 : 200)
    @show n
    for m ∈ max(1,n-10):n+10
      A = rand(T, m, n); res = similar(A);
      B = rand(T, n, n) + I
      @test TriangularSolve.rdiv!(res, A, UpperTriangular(B)) * UpperTriangular(B) ≈ A
      @test TriangularSolve.rdiv!(res, A, UnitUpperTriangular(B)) * UnitUpperTriangular(B) ≈ A
      @test TriangularSolve.rdiv!(res, A, UpperTriangular(B), Val(false)) * UpperTriangular(B) ≈ A
      @test TriangularSolve.rdiv!(res, A, UnitUpperTriangular(B), Val(false)) * UnitUpperTriangular(B) ≈ A
      A = rand(T, n, m); res = similar(A);
      @test LowerTriangular(B) * TriangularSolve.ldiv!(res, LowerTriangular(B), A) ≈ A
      @test UnitLowerTriangular(B) * TriangularSolve.ldiv!(res, UnitLowerTriangular(B), A) ≈ A
      @test LowerTriangular(B) * TriangularSolve.ldiv!(res, LowerTriangular(B), A, Val(false)) ≈ A
      @test UnitLowerTriangular(B) * TriangularSolve.ldiv!(res, UnitLowerTriangular(B), A, Val(false)) ≈ A
    end
  end
end

@testset "TriangularSolve.jl" begin
  @time @testset "Float64" begin
    test_solve(Float64)
  end
  @time @testset "Float32" begin
    test_solve(Float32)
  end
end
