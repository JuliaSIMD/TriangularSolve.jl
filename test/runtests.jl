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

      @test TriangularSolve.rdiv!(res, A, UpperTriangular(B)) *
            UpperTriangular(B) ≈ A
      @test TriangularSolve.rdiv!(res, A, UnitUpperTriangular(B)) *
            UnitUpperTriangular(B) ≈ A
      @test TriangularSolve.rdiv!(res, A, UpperTriangular(B), Val(false)) *
            UpperTriangular(B) ≈ A
      @test TriangularSolve.rdiv!(res, A, UnitUpperTriangular(B), Val(false)) *
            UnitUpperTriangular(B) ≈ A

      @test TriangularSolve.rdiv!(res, A, LowerTriangular(B)) *
            LowerTriangular(B) ≈ A
      @test TriangularSolve.rdiv!(res, A, UnitLowerTriangular(B)) *
            UnitLowerTriangular(B) ≈ A
      @test TriangularSolve.rdiv!(res, A, LowerTriangular(B), Val(false)) *
            LowerTriangular(B) ≈ A
      @test TriangularSolve.rdiv!(res, A, UnitLowerTriangular(B), Val(false)) *
            UnitLowerTriangular(B) ≈ A

      check_box_for_nans(RR, m, n)
      res .= NaN
      A .= NaN

      A = @view AA[17:16+n, 17:16+m]
      res = @view RR[17:16+n, 17:16+m]
      A .= rand.(T)

      @test LowerTriangular(B) *
            TriangularSolve.ldiv!(res, LowerTriangular(B), A) ≈ A
      @test UnitLowerTriangular(B) *
            TriangularSolve.ldiv!(res, UnitLowerTriangular(B), A) ≈ A
      @test LowerTriangular(B) *
            TriangularSolve.ldiv!(res, LowerTriangular(B), A, Val(false)) ≈ A
      @test UnitLowerTriangular(B) *
            TriangularSolve.ldiv!(res, UnitLowerTriangular(B), A, Val(false)) ≈
            A

      @test UpperTriangular(B) *
            TriangularSolve.ldiv!(res, UpperTriangular(B), A) ≈ A
      @test UnitUpperTriangular(B) *
            TriangularSolve.ldiv!(res, UnitUpperTriangular(B), A) ≈ A
      @test UpperTriangular(B) *
            TriangularSolve.ldiv!(res, UpperTriangular(B), A, Val(false)) ≈ A
      @test UnitUpperTriangular(B) *
            TriangularSolve.ldiv!(res, UnitUpperTriangular(B), A, Val(false)) ≈
            A
      check_box_for_nans(RR, n, m)
      res .= NaN
      A .= NaN
      B .= NaN
    end
  end
end

function test_packed_lu(::Type{T}) where {T}
  for n ∈ (8, 16, 33, 64, 129, 256, 512), nrhs ∈ (2, 5, 8)
    A = rand(T, n, n) + T(n) * I
    # only the triangle indicated by the wrapper may be read; the other
    # triangle holds the other LU factor
    F = lu!(copy(A)).factors
    B = rand(T, n, nrhs)
    for thread ∈ (Val(false), Val(true))
      X = TriangularSolve.ldiv!(UpperTriangular(F), copy(B), thread)
      @test X ≈ Matrix(UpperTriangular(F)) \ B
      Y = TriangularSolve.ldiv!(UnitLowerTriangular(F), copy(B), thread)
      @test Y ≈ Matrix(UnitLowerTriangular(F)) \ B
    end
    b = rand(T, n)
    for thread ∈ (Val(false), Val(true))
      xv = TriangularSolve.ldiv!(UpperTriangular(F), copy(b), thread)
      @test xv ≈ Matrix(UpperTriangular(F)) \ b
      yv = TriangularSolve.ldiv!(UnitLowerTriangular(F), copy(b), thread)
      @test yv ≈ Matrix(UnitLowerTriangular(F)) \ b
    end
    C = rand(T, nrhs, n)
    @test TriangularSolve.rdiv!(copy(C), UnitLowerTriangular(F)) ≈
          C / Matrix(UnitLowerTriangular(F))
  end
end

@testset "TriangularSolve.jl" begin
  @time @testset "Float64" begin
    test_solve(Float64)
  end
  @time @testset "Float32" begin
    test_solve(Float32)
  end
  @time @testset "packed LU factors" begin
    test_packed_lu(Float64)
    test_packed_lu(Float32)
  end
  @testset "native upper-ldiv / lower-rdiv dispatch" begin
    # these must hit TriangularSolve's own kernels, not the LinearAlgebra
    # catch-all fallbacks
    for T ∈ (Float32, Float64)
      @test which(
        TriangularSolve.ldiv!,
        (UpperTriangular{T,Matrix{T}}, Matrix{T}, Val{false})
      ).sig <: Tuple{
        typeof(TriangularSolve.ldiv!),
        UpperTriangular,
        AbstractMatrix,
        Val
      }
      @test which(
        TriangularSolve.rdiv!,
        (Matrix{T}, LowerTriangular{T,Matrix{T}}, Val{false})
      ).sig <: Tuple{
        typeof(TriangularSolve.rdiv!),
        AbstractMatrix,
        LowerTriangular,
        Val
      }
    end
  end
  @testset "vector right-hand sides" begin
    # every size runs the naive sweeps (the BLAS deferral is gone); sizes
    # cover SIMD-width remainders and the formerly-deferred n > 128 range
    for T ∈ (Float64, Float32),
      n ∈ (1, 2, 5, 8, 16, 33, 64, 127, 128, 129, 200, 500, 1201)

      P = rand(T, n, n) + T(n) * I
      # unit solves ignore the diagonal, so scale the strict triangle to
      # keep solution growth bounded at large n
      Pu = rand(T, n, n) ./ T(2n) + I
      b = rand(T, n)
      for wrap ∈ (
          UpperTriangular,
          UnitUpperTriangular,
          LowerTriangular,
          UnitLowerTriangular
        ),
        thread ∈ (Val(false), Val(true))

        unit = wrap === UnitUpperTriangular || wrap === UnitLowerTriangular
        U = wrap(unit ? Pu : P)
        x = TriangularSolve.ldiv!(U, copy(b), thread)
        @test x ≈ Matrix(U) \ b rtol = sqrt(eps(T)) * n
        c = similar(b)
        @test TriangularSolve.ldiv!(c, U, copy(b), thread) ≈ x
        xa = copy(b)
        @test TriangularSolve.ldiv!(xa, U, xa, thread) ≈ x
      end
    end
    F = lu!(rand(100, 100) + 100I).factors
    bf = rand(100)
    Uf = UpperTriangular(F)
    @test TriangularSolve.ldiv!(Uf, copy(bf), Val(false)) ≈ Matrix(Uf) \ bf
    xf = copy(bf)
    cf = similar(bf)
    for thread ∈ (Val(false), Val(true))
      TriangularSolve.ldiv!(Uf, xf, thread)
      xf .= bf
      @test iszero(@allocated TriangularSolve.ldiv!(Uf, xf, thread))
      TriangularSolve.ldiv!(cf, Uf, xf, thread)
      @test iszero(@allocated TriangularSolve.ldiv!(cf, Uf, xf, thread))
    end
  end
  @testset "non-contiguous strided vector solves" begin
    for T ∈ (Float64, Float32)
      n = 37
      M = rand(T, 2n, 2n) + T(2n) * I
      Pv = @view M[1:2:2n, 1:2:2n]
      v = rand(T, 2n)
      bv = @view v[1:2:2n]
      for wrap ∈ (
        UpperTriangular,
        UnitUpperTriangular,
        LowerTriangular,
        UnitLowerTriangular
      )
        U = wrap(Pv)
        xref = Matrix(U) \ Vector(bv)
        @test TriangularSolve.ldiv!(U, copy(bv)) ≈ xref rtol = sqrt(eps(T)) * n
        cv = @view similar(v)[1:2:2n]
        @test TriangularSolve.ldiv!(cv, U, bv) ≈ xref rtol = sqrt(eps(T)) * n
      end
    end
  end
  @testset "inner-product sweep kernels" begin
    # these forms are selected when the parent's rows are the contiguous
    # direction; no Base strided type has that layout, so test them directly
    for T ∈ (Float64, Float32), n ∈ (1, 7, 40), UNIT ∈ (false, true)
      P = rand(T, n, n) + T(n) * I
      b = rand(T, n)
      lo = UNIT ? UnitLowerTriangular(P) : LowerTriangular(P)
      up = UNIT ? UnitUpperTriangular(P) : UpperTriangular(P)
      x = copy(b)
      TriangularSolve._naive_vsolve_fwd_dot!(x, P, Val(UNIT))
      @test x ≈ Matrix(lo) \ b rtol = sqrt(eps(T)) * n
      x = copy(b)
      TriangularSolve._naive_vsolve_bwd_dot!(x, P, Val(UNIT))
      @test x ≈ Matrix(up) \ b rtol = sqrt(eps(T)) * n
    end
  end
  @testset "dimension mismatch throws" begin
    @test_throws DimensionMismatch TriangularSolve.rdiv!(
      rand(4, 8),
      LowerTriangular(rand(6, 6) + 6I),
      Val(false)
    )
    @test_throws DimensionMismatch TriangularSolve.ldiv!(
      UpperTriangular(rand(6, 6) + 6I),
      rand(8, 3),
      Val(false)
    )
    @test_throws DimensionMismatch TriangularSolve.ldiv!(
      rand(4, 3),
      UpperTriangular(rand(6, 6) + 6I),
      rand(6, 3),
      Val(false)
    )
    for wrap ∈ (LowerTriangular, UnitUpperTriangular)
      U = wrap(rand(6, 6) + 6I)
      @test_throws DimensionMismatch TriangularSolve.ldiv!(U, rand(8))
      @test_throws DimensionMismatch TriangularSolve.ldiv!(rand(6), U, rand(8))
      @test_throws DimensionMismatch TriangularSolve.ldiv!(rand(4), U, rand(6))
    end
  end
  @testset "non-strided inputs keep the LinearAlgebra fallback" begin
    # a Bidiagonal parent is an AbstractMatrix{Float64} but not strided; it
    # must keep hitting the LinearAlgebra catch-all rather than the SIMD path
    dv = rand(8) .+ 8
    ev = rand(7)
    Bd = Bidiagonal(dv, ev, :U)
    B = rand(8, 4)
    @test TriangularSolve.ldiv!(UpperTriangular(Bd), copy(B)) ≈
          Matrix(UpperTriangular(Bd)) \ B
    bd = rand(8)
    @test TriangularSolve.ldiv!(UpperTriangular(Bd), copy(bd)) ≈
          Matrix(UpperTriangular(Bd)) \ bd
    Pb = big.(rand(8, 8)) + 8I
    bb = big.(rand(8))
    @test TriangularSolve.ldiv!(LowerTriangular(Pb), copy(bb)) ≈
          Matrix(LowerTriangular(Pb)) \ bb
  end
  @testset "allocations" begin
    n = 200
    U = UpperTriangular(rand(n, n) + n * I)
    B = rand(n, 8)
    X = copy(B)
    for thread ∈ (Val(false), Val(true))
      TriangularSolve.ldiv!(U, X, thread)
      X .= B
      @test iszero(@allocated TriangularSolve.ldiv!(U, X, thread))
    end
    A = rand(8, n)
    L = LowerTriangular(parent(U))
    TriangularSolve.rdiv!(A, L, Val(false))
    @test iszero(@allocated TriangularSolve.rdiv!(A, L, Val(false)))
  end
end

using Aqua
Aqua.test_all(TriangularSolve; ambiguities = false)
@test isempty(Test.detect_ambiguities(TriangularSolve))
