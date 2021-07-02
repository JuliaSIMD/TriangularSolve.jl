using TriangularSolve
using Test

@testset "TriangularSolve.jl" begin

  @time @testset "Float64" begin
    for n ∈ 1:100
      for m ∈ max(1,n-10):n+10
        A = rand(m, n);
        U = TriangularSolve.UpperTriangular(rand(n, n));
        Bref = A / U
        @test TriangularSolve.rdiv!(similar(Bref), A, U) ≈ Bref
      end
    end
  end
  @time @testset "Float32" begin
    for n ∈ 1:100
      for m ∈ max(1,n-10):n+10
        A = rand(Float32, m, n);
        U = TriangularSolve.UpperTriangular(rand(Float32, n, n));
        Bref = A / U
        @test TriangularSolve.rdiv!(similar(Bref), A, U) ≈ Bref
      end
    end
  end
end
