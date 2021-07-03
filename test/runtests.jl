using TriangularSolve, LinearAlgebra
using Test

function test_solve(::Type{T}) where {T}
  for n ∈ 1:100
    for m ∈ max(1,n-10):n+10
      A = rand(T, m, n); res = similar(A);
      B = rand(T, n, n);      
      @test TriangularSolve.rdiv!(res, A, UpperTriangular(B)) ≈ A / UpperTriangular(B)
      @test TriangularSolve.rdiv!(res, A, UnitUpperTriangular(B)) ≈ A / UnitUpperTriangular(B)
      A = rand(T, n, m); res = similar(A);
      @test TriangularSolve.ldiv!(res, LowerTriangular(B), A) ≈ LowerTriangular(B) \ A
      @test TriangularSolve.ldiv!(res, UnitLowerTriangular(B), A) ≈ UnitLowerTriangular(B) \ A
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
