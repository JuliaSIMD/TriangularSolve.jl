module TriangularSolve

using VectorizationBase, LinearAlgebra #LoopVectorization
using VectorizationBase: vfnmadd_fast, AbstractStridedPointer, AbstractMask, stridedpointer_preserve, zero_offsets, gesp, StridedPointer
using CloseOpenIntervals: CloseOpen, SafeCloseOpen

@generated function solve_AU(A::VecUnroll{Nm1}, spu::AbstractStridedPointer, noff, ::Val{UNIT}) where {Nm1, UNIT}
  A_n_expr = UNIT ? :nothing : :(A_n = Base.FastMath.div_fast(A_n, U_n_n))
  N = Nm1 + 1
  quote
    $(Expr(:meta,:inline))
    Ad = VectorizationBase.data(A)
    Base.Cartesian.@nexprs $N n -> begin
      A_n = Ad[n]
      Base.Cartesian.@nexprs $(UNIT ? :(n-1) : :n) m -> begin
        U_m_n = vload(spu, (noff+(m-1),noff+(n-1)))
      end
    end
    Base.Cartesian.@nexprs $N n -> begin
      Base.Cartesian.@nexprs n-1 k -> begin
        A_n = Base.FastMath.sub_fast(A_n, Base.FastMath.mul_fast(A_k, U_k_n))
      end
      $A_n_expr
    end
    VecUnroll(Base.Cartesian.@ntuple $N A)
  end
end

# @generated function nmuladd(A::VecUnroll{Nm1},B::AbstractStridedPointer,C::VecUnroll{Nm1}) where {Nm1}
#   N = Nm1 + 1
#   quote
#     $(Expr(:meta,:inline))
#     Ad = VectorizationBase.data(A);
#     Cd = VectorizationBase.data(C);
#     bp = stridedpointer(B)
#     Base.Cartesian.@nexprs $N n -> C_n = Cd[n]
#     Base.Cartesian.@nexprs $N k -> begin
#       A_k = Ad[k]
#       Base.Cartesian.@nexprs $N n -> begin
#         C_n = Base.FastMath.sub_fast(C_n, Base.FastMath.mul_fast(A_k, vload(B, (k-1,n-1))))
#       end
#     end
#     VecUnroll(Base.Cartesian.@ntuple $N C)
#   end
# end

# @inline function solve_Wx3W(A11::V, A12::V, A13::V, U::AbstractMatrix, ::StaticInt{W}) where {V<:VecUnroll,W}
#   WS = StaticInt{W}()

#   U11 = view(U,StaticInt(1):WS,StaticInt(1):WS)
#   A11 = solve_AU(A11, U11)

#   U12 = view(U,StaticInt(1):WS, StaticInt(1)+WS:WS*StaticInt(2))
#   A12 = nmuladd(A11, U12, A12)
#   U22 = view(U,StaticInt(1)+WS:WS*StaticInt(2),StaticInt(1)+WS:WS*StaticInt(2))
#   A12 = solve_AU(A12, U22)
  
#   U13 = view(U,StaticInt(1):WS, StaticInt(1)+WS*StaticInt(2):WS*StaticInt(3))
#   A13 = nmuladd(A11, U13, A13)
#   U23 = view(U,StaticInt(1)+WS:WS*StaticInt(2),StaticInt(1)+WS*StaticInt(2):WS*StaticInt(3))
#   A13 = nmuladd(A12, U23, A13)
#   U33 = view(U,StaticInt(1)+WS*StaticInt(2):WS*StaticInt(3),StaticInt(1)+WS*StaticInt(2):WS*StaticInt(3))
#   A13 = solve_AU(A13, U33)

#   return A11, A12, A13
# end

# @inline function solve_Wx3W!(ap::AbstractStridedPointer{T}, bp::AbstractStridedPointer{T}, U, rowoffset, coloffset) where {T}
#   WS = VectorizationBase.pick_vector_width(T)
#   W = Int(WS)  
#   A11 = vload(bp, Unroll{2,1,W,1,W,0x0000000000000000,1}((rowoffset,coloffset)))
#   A12 = vload(bp, Unroll{2,1,W,1,W,0x0000000000000000,1}((rowoffset,coloffset+WS)))
#   A13 = vload(bp, Unroll{2,1,W,1,W,0x0000000000000000,1}((rowoffset,coloffset+WS+WS)))

#   A11, A12, A13 = solve_Wx3W(A11, A12, A13, U, WS)
  
#   vstore!(ap, A11, Unroll{2,1,W,1,W,0x0000000000000000,1}((rowoffset,coloffset)))
#   vstore!(ap, A12, Unroll{2,1,W,1,W,0x0000000000000000,1}((rowoffset,coloffset+WS)))
#   vstore!(ap, A13, Unroll{2,1,W,1,W,0x0000000000000000,1}((rowoffset,coloffset+WS+WS)))
# end
# @inline function solve_Wx3W!(ap::AbstractStridedPointer{T}, bp::AbstractStridedPointer{T}, U, rowoffset, coloffset, m::VectorizationBase.AbstractMask) where {T}
#   WS = VectorizationBase.pick_vector_width(T)
#   W = Int(WS)
#   A11 = vload(bp, Unroll{2,1,W,1,W,0xffffffffffffffff,1}((rowoffset,coloffset)), m)
#   A12 = vload(bp, Unroll{2,1,W,1,W,0xffffffffffffffff,1}((rowoffset,coloffset+WS)), m)
#   A13 = vload(bp, Unroll{2,1,W,1,W,0xffffffffffffffff,1}((rowoffset,coloffset+WS+WS)), m)

#   A11, A12, A13 = solve_Wx3W(A11, A12, A13, U, WS)

#   vstore!(ap, A11, Unroll{2,1,W,1,W,0xffffffffffffffff,1}((rowoffset,coloffset)), m)
#   vstore!(ap, A12, Unroll{2,1,W,1,W,0xffffffffffffffff,1}((rowoffset,coloffset+WS)), m)
#   vstore!(ap, A13, Unroll{2,1,W,1,W,0xffffffffffffffff,1}((rowoffset,coloffset+WS+WS)), m)
# end

# solve_3Wx3W!(A,B,U::UpperTriangular) = solve_3Wx3W!(A,B,parent(U))
# function solve_3Wx3W!(A::AbstractMatrix{T},B,U) where {T}
#   W = VectorizationBase.pick_vector_width(T)
#   ap = stridedpointer(A);
#   bp = stridedpointer(B);
#   solve_Wx3W!(ap, bp, U, StaticInt(1),         StaticInt(1))
#   solve_Wx3W!(ap, bp, U, StaticInt(1) + W,     StaticInt(1))
#   solve_Wx3W!(ap, bp, U, StaticInt(1) + W + W, StaticInt(1))
# end

@inline maybestore!(p, v, i) = vstore!(p, v, i)
@inline maybestore!(::Nothing, v, i) = nothing

@inline maybestore!(p, v, i, m) = vstore!(p, v, i, m)
@inline maybestore!(::Nothing, v, i, m) = nothing

@inline function store_small_kern!(spa, sp, v, spu, i, n, mask, ::Val{true})
  vstore!(spa, v, i, mask)
  vstore!(sp, v, i, mask)
end
@inline store_small_kern!(spa, ::Nothing, v, spu, i, n, mask, ::Val{true}) = vstore!(spa, v, i, mask)

@inline function store_small_kern!(spa, sp, v, spu, i, n, mask, ::Val{false})
  x = v / vload(spu, (n,n))
  vstore!(spa, x, i, mask)
  vstore!(sp, x, i, mask)
end
@inline store_small_kern!(spa, ::Nothing, v, spu, i, n, mask, ::Val{false}) = vstore!(spa, v / vload(spu, (n,n)), i, mask)

@inline function store_small_kern!(spa, sp, v, spu, i, n, ::Val{true})
  vstore!(spa, v, i)
  vstore!(sp, v, i)
end
@inline store_small_kern!(spa, ::Nothing, v, spu, i, n, ::Val{true}) = vstore!(spa, v, i)

@inline function store_small_kern!(spa, sp, v, spu, i, n, ::Val{false})
  x = v / vload(spu, (n,n))
  vstore!(spa, x, i)
  vstore!(sp, x, i)
end
@inline store_small_kern!(spa, ::Nothing, v, spu, i, n, ::Val{false}) = vstore!(spa, v / vload(spu, (n,n)), i)

function BdivU_small_kern!(spa::AbstractStridedPointer{T}, sp, spb::AbstractStridedPointer{T}, spu::AbstractStridedPointer{T}, N, mask, ::Val{UNIT}) where {T,UNIT}
  W = VectorizationBase.pick_vector_width(T)
  for n ∈ CloseOpen(N)
    Amn = vload(spb, (MM(W, StaticInt(0)),n), mask)
    for k ∈ SafeCloseOpen(n)
      Amn = vfnmadd_fast(vload(spa, (MM(W, StaticInt(0)),k), mask), vload(spu, (k,n)), Amn)
    end
    store_small_kern!(spa, sp, Amn, spu, (MM(W, StaticInt(0)),n), n, mask, Val{UNIT}())
  end
end
function BdivU_small_kern_u3!(spa::AbstractStridedPointer{T}, sp, spb::AbstractStridedPointer{T}, spu::AbstractStridedPointer{T}, N, ::Val{UNIT}) where {T,UNIT}
  W = Int(VectorizationBase.pick_vector_width(T))
  for n ∈ CloseOpen(N)
    Amn = vload(spb, Unroll{1,W,3,1,W,0x0000000000000000,1}((StaticInt(0),n)))
    for k ∈ SafeCloseOpen(n)
      Amk = vload(spa, Unroll{1,W,3,1,W,0x0000000000000000,1}((StaticInt(0),k)))
      Amn = vfnmadd_fast(Amk, vload(spu, (k,n)), Amn)
    end
    store_small_kern!(spa, sp, Amn, spu, Unroll{1,W,3,1,W,0x0000000000000000,1}((StaticInt(0),n)), n, Val{UNIT}())
  end
end
# function BdivU_small!(A::AbstractMatrix{T}, B::AbstractMatrix{T}, U::AbstractMatrix{T}) where {T}
#   W = VectorizationBase.pick_vector_width(T)
#   M, N = size(A)
#   m = 0
#   spa = stridedpointer(A)
#   spb = stridedpointer(B)
#   spu = stridedpointer(U)
#   while m < M
#     ml = m+1
#     mu = m+W
#     maskiter = mu > M
#     mask = maskiter ? VectorizationBase.mask(W, M) : VectorizationBase.max_mask(W)
#     for n ∈ 1:N
#       Amn = vload(spb, (MM(W, ml),n), mask)
#       for k ∈ 1:n-1
#         Amn = vfnmadd_fast(vload(spa, (MM(W, ml),k), mask), vload(spu, (k,n)), Amn)
#       end
#       vstore!(spa, Amn / vload(spu, (n,n)), (MM(W, ml),n), mask)
#     end    
#     m = mu
#   end
#   # @inbounds @fastmath for m ∈ 1:M
#   #   for n ∈ 1:N
#   #     Amn = B[m,n]
#   #     for k ∈ 1:n-1
#   #       Amn -= A[m,k]*U[k,n]
#   #     end
#   #     A[m,n] = Amn / U[n,n]
#   #   end
#   # end
# end
# function nmuladd!(C,A,B,D)
#   @turbo for n ∈ axes(C,2), m ∈ axes(C,1)
#     Cmn = D[m,n]
#     for k ∈ axes(B,1)
#       Cmn -= A[m,k]*B[k,n]
#     end
#     C[m,n] = Cmn
#   end
# end

@generated function rdiv_solve_W_u3!(spc, spb, spa, spu, n, ::StaticInt{W}, ::Val{UNIT}) where {W, UNIT}
  quote
    # $(Expr(:meta,:inline))
    # here, we just want to load the vectors
    C11 = VectorizationBase.data(vload(spa, Unroll{2,1,$W,1,$W,0x0000000000000000,1}(Unroll{1,$W,3,1,$W,0x0000000000000000,1}((StaticInt(0),n)))))
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    for nk ∈ SafeCloseOpen(n) # nmuladd
      A11 = vload(spc, Unroll{1,$W,3,1,$W,0x0000000000000000,1}((StaticInt(0),nk)))
      Base.Cartesian.@nexprs $W c -> C11_c = vfnmadd_fast(A11, vload(spu, (nk,n+(c-1))), C11_c)
    end
    C11vu = solve_AU(VecUnroll((Base.Cartesian.@ntuple $W C11)), spu, n, Val{$UNIT}())
    i = Unroll{2,1,$W,1,$W,0x0000000000000000,1}(Unroll{1,$W,3,1,$W,0x0000000000000000,1}((StaticInt(0),n)))
    vstore!(spc, C11vu, i)
    maybestore!(spb, C11vu, i)
  end
end
@generated function rdiv_solve_W_u1!(spc, spb, spa, spu, n, dontstorec, mask::AbstractMask{W}, ::Val{UNIT}) where {W, UNIT}
  quote
    # $(Expr(:meta,:inline))
    # here, we just want to load the vectors
    C11 = VectorizationBase.data(vload(spa, Unroll{2,1,$W,1,$W,0xffffffffffffffff,1}((StaticInt(0),n)), mask))
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    for nk ∈ SafeCloseOpen(n) # nmuladd
      A11 = vload(spc, (MM{$W}(StaticInt(0)),nk), mask)
      Base.Cartesian.@nexprs $W c -> C11_c =  vfnmadd_fast(A11, vload(spu, (nk,n+(c-1))), C11_c)
    end
    C11 = VecUnroll((Base.Cartesian.@ntuple $W C11))
    C11 = solve_AU(C11, spu, n, Val{$UNIT}())
    i = Unroll{2,1,$W,1,$W,0xffffffffffffffff,1}((StaticInt(0),n))
    dontstorec ≢ false && vstore!(spc, C11, i, mask)
    maybestore!(spb, C11, i, mask)
  end
end

function rdiv_U!(spc::AbstractStridedPointer{T}, spa, spu, M, N, ::StaticInt{1}, ::Val{UNIT}) where {T,UNIT}
  WS = pick_vector_width(T)
  W = Int(WS)
  W3 = StaticInt(3)*WS
  if VectorizationBase.register_count() > WS*StaticInt(3)
    # We have enough registers to use the 3W functions.
    # Md, Mr = VectorizationBase.vdivrem(M, W3)
    M3 = M
    Nd, Nr = VectorizationBase.vdivrem(N, W3)
  else
    # We do not have enough registers.
    M3 = 0
    Nd = 0;
    # Mr = M;
    Nr = N;
  end
  # Mrd8, Mrr8 = VectorizationBase.vdivrem(Mr, WS)
  Nrd8, Nrr8 = VectorizationBase.vdivrem(Nr, WS)
  Nrd8main = Nrd8+Nd*3
  m = 0
  while m < M3 - W3 + 1
    n = 0
    if Nrr8 > 0
      n = Nrr8
      BdivU_small_kern_u3!(spc, nothing, spa, spu, Nrr8, Val(UNIT))
      # BdivU_small!(view(C,Mrange,1:n),view(A,Mrange,1:n),view(U,1:n,1:n))
    end
    for i ∈ 1:Nrd8main
      rdiv_solve_W_u3!(spc, nothing, spa, spu, n, WS, Val(UNIT))
      n += W
    end
    m += W3
    spa = gesp(spa, (W3,StaticInt(0)))
    spc = gesp(spc, (W3,StaticInt(0)))
  end
  finalmask = VectorizationBase.mask(WS, M)
  while m < M
    ubm = m+W
    nomaskiter = ubm < M
    mask = nomaskiter ? VectorizationBase.max_mask(WS) : finalmask
    n = Nrr8
    if n > 0
      BdivU_small_kern!(spc, nothing, spa, spu, n, mask, Val(UNIT))
    end
    for i ∈ 1:Nrd8main
      rdiv_solve_W_u1!(spc, nothing, spa, spu, n, nothing, mask, Val(UNIT))
      n += W
    end
    spa = gesp(spa, (WS,StaticInt(0)))
    spc = gesp(spc, (WS,StaticInt(0)))
    m = ubm
  end
  nothing
end

const LDIVBUFFERS = Vector{UInt8}[]
function lubuffer(::Val{T}, N) where {T}
  buff = LDIVBUFFERS[Threads.threadid()]
  RS3 = StaticInt{3}()*VectorizationBase.register_size()
  L = RS3*N
  L > length(buff) && resize!(buff, L)
  ptr = Base.unsafe_convert(Ptr{T}, buff)
  StridedPointer{T,2,1,0,(1,2)}(ptr, (VectorizationBase.static_sizeof(T), RS3), (StaticInt(0),StaticInt(0)))
end
rdiv!(A, U::UpperTriangular) = (div_dispatch!(A, A, parent(U), Val(false)); return A)
rdiv!(C, A, U::UpperTriangular) = (div_dispatch!(C, A, parent(U), Val(false)); return C)
rdiv!(A, U::UnitUpperTriangular) = (div_dispatch!(A, A, parent(U), Val(true)); return A)
rdiv!(C, A, U::UnitUpperTriangular) = (div_dispatch!(C, A, parent(U), Val(true)); return C)

ldiv!(U::LowerTriangular, A) = (div_dispatch!(A', A', parent(U)', Val(false)); return A)
ldiv!(C, U::LowerTriangular, A) = (div_dispatch!(C', A', parent(U)', Val(false)); return C)
ldiv!(U::UnitLowerTriangular, A) = (div_dispatch!(A', A', parent(U)', Val(true)); return A)
ldiv!(C, U::UnitLowerTriangular, A) = (div_dispatch!(C', A', parent(U)', Val(true)); return C)

function div_dispatch!(C, A, U, ::Val{UNIT}) where {UNIT}
  M, N = size(A)
  ((N == 0) | (M == 0)) && return C
  _spa, spap = stridedpointer_preserve(A)
  _spc, spcp = stridedpointer_preserve(C)
  _spu, spup = stridedpointer_preserve(U)
  spa = zero_offsets(_spa)
  spc = zero_offsets(_spc)
  spu = zero_offsets(_spu)
  GC.@preserve spap spcp spup begin
    # if N ≥ 32
    rdiv_U!(spc, spa, spu, M, N, VectorizationBase.contiguous_axis(A), Val(UNIT))
    # else
      # rdiv_U!(spc, spa, spu, M, N, StaticInt(1), Val(UNIT))
    # end
  end
end

function rdiv_U!(spc::AbstractStridedPointer{T}, spa, spu, M, N, ::StaticInt, ::Val{UNIT}) where {T,UNIT}
  WS = pick_vector_width(T)
  W = Int(WS)
  W3 = StaticInt(3)*WS
  if VectorizationBase.register_count() > WS*StaticInt(3)
    # We have enough registers to use the 3W functions.
    M3 = M
    # Md, Mr = VectorizationBase.vdivrem(M, W3)
    Nd, Nr = VectorizationBase.vdivrem(N, W3)
  else
    # We do not have enough registers.
    Nd = 0;
    M3 = 0
    Nr = N;
  end
  # Mrd8, Mrr8 = VectorizationBase.vdivrem(Mr, WS)
  Nrd8, Nrr8 = VectorizationBase.vdivrem(Nr, WS)
  Nrd8main = Nrd8+Nd*3
  # if Nrd8 == 0
  #   spb = spa
  # else
  spb = lubuffer(Val(T), N)
  # @show M3, M, Nrd8main, Nrr8
  # end
  m = 0
  while m < M3 - W3 + 1
    n = Nrr8
    if n > 0
      BdivU_small_kern_u3!(spb, spc, spa, spu, Nrr8, Val(UNIT))
      # BdivU_small!(view(C,Mrange,1:n),view(A,Mrange,1:n),view(U,1:n,1:n))
    end
    for i ∈ 1:Nrd8main
      rdiv_solve_W_u3!(spb, spc, spa, spu, n, WS, Val(UNIT))
      n += W
    end
    m += W3
    spa = gesp(spa, (W3,StaticInt(0)))
    spc = gesp(spc, (W3,StaticInt(0)))
  end
  finalmask = VectorizationBase.mask(WS, M)
  while m < M
    ubm = m+W
    nomaskiter = ubm < M
    mask = nomaskiter ? VectorizationBase.max_mask(WS) : finalmask
    n = Nrr8
    if n > 0
      BdivU_small_kern!(spb, spc, spa, spu, n, mask, Val(UNIT))
    end
    for i ∈ 1:Nrd8main
      # @show C, n
      rdiv_solve_W_u1!(spb, spc, spa, spu, n, i ≠ Nrd8main, mask, Val(UNIT))
      n += W
    end
    spa = gesp(spa, (WS,StaticInt(0)))
    spc = gesp(spc, (WS,StaticInt(0)))
    m = ubm
  end
  nothing
end


function __init__()
  nthread = Threads.nthreads()
  resize!(LDIVBUFFERS, nthread)
  for i ∈ 1:nthread
    LDIVBUFFERS[i] = Vector{UInt8}(undef, 3VectorizationBase.register_size()*128)
  end
end

end
