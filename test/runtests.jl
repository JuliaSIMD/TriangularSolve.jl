using TriangularSolve, LinearAlgebra
using Test

function check_box_for_nans(A, M, N)
  # blocks start at 17, and are MxN
  @test all(isnan, @view(A[1:16, :]))
  @test all(isnan, @view(A[17+M:end, :]))
  @test all(isnan, @view(A[17:16+M, 1:16]))
  @test all(isnan, @view(A[17:16+M, 17+N:end]))
end

function test_solve(::Type{T}) where {T}
  maxN = (T === Float32 ? 100 : 200)
  maxM = maxN + 10
  AA = fill(T(NaN), maxM + 32, maxM + 32)
  RR = fill(T(NaN), maxM + 32, maxM + 32)
  BB = fill(T(NaN), maxN + 32, maxN + 32)
  for n ∈ 1:maxN
    @show n
    for m ∈ max(1, n - 10):n+10
      A = @view AA[17:16+m, 17:16+n]
      res = @view RR[17:16+m, 17:16+n]
      B = @view BB[17:16+n, 17:16+n]

      A .= rand.(T)
      B .= rand.(T)
      @view(B[diagind(B)]) .+= one(T)

      for C in (
        UpperTriangular(B),
        UnitUpperTriangular(B),
        LowerTriangular(B),
        UnitLowerTriangular(B)
      )
        @test TriangularSolve.rdiv!(res, A, C) * C ≈ A
        check_box_for_nans(RR, m, n)
        @test TriangularSolve.rdiv!(res, A, C, Val(false)) * C ≈ A
        check_box_for_nans(RR, m, n)
      end

      res .= NaN
      A .= NaN

      A = @view AA[17:16+n, 17:16+m]
      res = @view RR[17:16+n, 17:16+m]
      A .= rand.(T)

      for C in (
        UpperTriangular(B),
        UnitUpperTriangular(B),
        LowerTriangular(B),
        UnitLowerTriangular(B)
      )
        @test C * TriangularSolve.ldiv!(res, C, A) ≈ A
        check_box_for_nans(RR, n, m)
        @test C * TriangularSolve.ldiv!(res, C, A, Val(false)) ≈ A
        check_box_for_nans(RR, n, m)
      end

      res .= NaN
      A .= NaN
      B .= NaN
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

using Aqua
Aqua.test_all(TriangularSolve; ambiguities = false)
@test isempty(Test.detect_ambiguities(TriangularSolve))
