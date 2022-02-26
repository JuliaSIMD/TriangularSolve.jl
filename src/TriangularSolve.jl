module TriangularSolve

using LayoutPointers: stridedpointer_preserve, StrideIndex
using VectorizationBase, LinearAlgebra #LoopVectorization
using VectorizationBase: vfnmadd_fast, AbstractStridedPointer, AbstractMask, zero_offsets, gesp, StridedPointer
using CloseOpenIntervals: CloseOpen, SafeCloseOpen
using Static
using IfElse: ifelse
using LoopVectorization
using Polyester

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
#   A11 = vload(bp, Unroll{2,1,W,1,W,zero(UInt),1}((rowoffset,coloffset)))
#   A12 = vload(bp, Unroll{2,1,W,1,W,zero(UInt),1}((rowoffset,coloffset+WS)))
#   A13 = vload(bp, Unroll{2,1,W,1,W,zero(UInt),1}((rowoffset,coloffset+WS+WS)))

#   A11, A12, A13 = solve_Wx3W(A11, A12, A13, U, WS)
  
#   vstore!(ap, A11, Unroll{2,1,W,1,W,zero(UInt),1}((rowoffset,coloffset)))
#   vstore!(ap, A12, Unroll{2,1,W,1,W,zero(UInt),1}((rowoffset,coloffset+WS)))
#   vstore!(ap, A13, Unroll{2,1,W,1,W,zero(UInt),1}((rowoffset,coloffset+WS+WS)))
# end
# @inline function solve_Wx3W!(ap::AbstractStridedPointer{T}, bp::AbstractStridedPointer{T}, U, rowoffset, coloffset, m::VectorizationBase.AbstractMask) where {T}
#   WS = VectorizationBase.pick_vector_width(T)
#   W = Int(WS)
#   A11 = vload(bp, Unroll{2,1,W,1,W,(-1%UInt),1}((rowoffset,coloffset)), m)
#   A12 = vload(bp, Unroll{2,1,W,1,W,(-1%UInt),1}((rowoffset,coloffset+WS)), m)
#   A13 = vload(bp, Unroll{2,1,W,1,W,(-1%UInt),1}((rowoffset,coloffset+WS+WS)), m)

#   A11, A12, A13 = solve_Wx3W(A11, A12, A13, U, WS)

#   vstore!(ap, A11, Unroll{2,1,W,1,W,(-1%UInt),1}((rowoffset,coloffset)), m)
#   vstore!(ap, A12, Unroll{2,1,W,1,W,(-1%UInt),1}((rowoffset,coloffset+WS)), m)
#   vstore!(ap, A13, Unroll{2,1,W,1,W,(-1%UInt),1}((rowoffset,coloffset+WS+WS)), m)
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

@inline function BdivU_small_kern!(spa::AbstractStridedPointer{T}, sp, spb::AbstractStridedPointer{T}, spu::AbstractStridedPointer{T}, N, mask::AbstractMask{W}, ::Val{UNIT}) where {T,UNIT,W}
  # W = VectorizationBase.pick_vector_width(T)
  for n ∈ CloseOpen(N)
    Amn = vload(spb, (MM{W}(StaticInt(0)),n), mask)
    for k ∈ SafeCloseOpen(n)
      Amn = vfnmadd_fast(vload(spa, (MM{W}(StaticInt(0)),k), mask), vload(spu, (k,n)), Amn)
    end
    store_small_kern!(spa, sp, Amn, spu, (MM{W}(StaticInt(0)),n), n, mask, Val{UNIT}())
  end
end
@inline function BdivU_small_kern_u!(spa::AbstractStridedPointer{T}, sp, spb::AbstractStridedPointer{T}, spu::AbstractStridedPointer{T}, N, ::StaticInt{U}, ::Val{UNIT}) where {T,U,UNIT}
  W = Int(VectorizationBase.pick_vector_width(T))
  for n ∈ CloseOpen(N)
    Amn = vload(spb, Unroll{1,W,U,1,W,zero(UInt),1}((StaticInt(0),n)))
    for k ∈ SafeCloseOpen(n)
      Amk = vload(spa, Unroll{1,W,U,1,W,zero(UInt),1}((StaticInt(0),k)))
      Amn = vfnmadd_fast(Amk, vload(spu, (k,n)), Amn)
    end
    store_small_kern!(spa, sp, Amn, spu, Unroll{1,W,U,1,W,zero(UInt),1}((StaticInt(0),n)), n, Val{UNIT}())
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

@generated function rdiv_solve_W_u!(spc, spb, spa, spu, n, ::StaticInt{W}, ::StaticInt{U}, ::Val{UNIT}) where {W, U, UNIT}
  quote
    $(Expr(:meta,:inline))
    # here, we just want to load the vectors
    C11 = VectorizationBase.data(vload(spa, Unroll{2,1,$W,1,$W,zero(UInt),1}(Unroll{1,$W,$U,1,$W,zero(UInt),1}((StaticInt(0),n)))))
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    for nk ∈ SafeCloseOpen(n) # nmuladd
      A11 = vload(spc, Unroll{1,$W,$U,1,$W,zero(UInt),1}((StaticInt(0),nk)))
      Base.Cartesian.@nexprs $W c -> C11_c = vfnmadd_fast(A11, vload(spu, (nk,n+(c-1))), C11_c)
    end
    C11vu = solve_AU(VecUnroll((Base.Cartesian.@ntuple $W C11)), spu, n, Val{$UNIT}())
    i = Unroll{2,1,$W,1,$W,zero(UInt),1}(Unroll{1,$W,$U,1,$W,zero(UInt),1}((StaticInt(0),n)))
    vstore!(spc, C11vu, i)
    maybestore!(spb, C11vu, i)
  end
end
@generated function rdiv_solve_W!(spc, spb, spa, spu, n, storec::B, mask::AbstractMask{W}, ::Val{UNIT}) where {W, UNIT, B}
  storecexpr = if (B <: Bool)
    :(storec && vstore!(spc, C11, i, mask))
  else
    :(vstore!(spc, C11, i, mask))
  end
  quote
    $(Expr(:meta,:inline))
    # here, we just want to load the vectors
    C11 = VectorizationBase.data(vload(spa, Unroll{2,1,$W,1,$W,(-1%UInt),1}((StaticInt(0),n)), mask))
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    for nk ∈ SafeCloseOpen(n) # nmuladd
      A11 = vload(spc, (MM{$W}(StaticInt(0)),nk), mask)
      Base.Cartesian.@nexprs $W c -> C11_c =  vfnmadd_fast(A11, vload(spu, (nk,n+(c-1))), C11_c)
    end
    C11 = VecUnroll((Base.Cartesian.@ntuple $W C11))
    C11 = solve_AU(C11, spu, n, Val{$UNIT}())
    i = Unroll{2,1,$W,1,$W,(-1%UInt),1}((StaticInt(0),n))
    $storecexpr
    maybestore!(spb, C11, i, mask)
  end
end

@inline function rdiv_U!(spc::AbstractStridedPointer{T}, spa::AbstractStridedPointer, spu::AbstractStridedPointer, M, N, ::StaticInt{1}, ::Val{UNIT}) where {T,UNIT}
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF*WS
  MU = UF > 1 ? M : 0
  Nd, Nr = VectorizationBase.vdivrem(N, WS)
  m = 0
  while m < MU - WU + 1
    n = Nr
    if n > 0
      BdivU_small_kern_u!(spc, nothing, spa, spu, n, UF, Val(UNIT))
    end
    for i ∈ 1:Nd
      rdiv_solve_W_u!(spc, nothing, spa, spu, n, WS, UF, Val(UNIT))
      n += W
    end
    m += WU
    spa = gesp(spa, (WU,StaticInt(0)))
    spc = gesp(spc, (WU,StaticInt(0)))
  end
  finalmask = VectorizationBase.mask(WS, M)
  while m < M
    ubm = m+W
    nomaskiter = ubm < M
    mask = nomaskiter ? VectorizationBase.max_mask(WS) : finalmask
    n = Nr
    if n > 0
      BdivU_small_kern!(spc, nothing, spa, spu, n, mask, Val(UNIT))
    end
    for _ ∈ 1:Nd
      rdiv_solve_W!(spc, nothing, spa, spu, n, nothing, mask, Val(UNIT))
      n += W
    end
    spa = gesp(spa, (WS,StaticInt(0)))
    spc = gesp(spc, (WS,StaticInt(0)))
    m = ubm
  end
  nothing
end

const LDIVBUFFERS = Vector{UInt8}[]
@inline function lubuffer(::Val{T}, ::StaticInt{UF}, N) where {T, UF}
  buff = LDIVBUFFERS[Threads.threadid()]
  RSUF = StaticInt{UF}()*VectorizationBase.register_size()
  L = RSUF*N
  L > length(buff) && resize!(buff, L)
  ptr = Base.unsafe_convert(Ptr{T}, buff)
  si = StrideIndex{2,(1,2),1}((VectorizationBase.static_sizeof(T), RSUF), (StaticInt(0),StaticInt(0)))
  stridedpointer(ptr, si, StaticInt{0}())
end
_canonicalize(x) = signed(x)
_canonicalize(::StaticInt{N}) where {N} = StaticInt{N}()
function div_dispatch!(C::AbstractMatrix{T}, A, U, ::Val{UNIT}, ::Val{THREAD}) where {UNIT,T,THREAD}
  _M, _N = size(A)
  M = _canonicalize(_M)
  N = _canonicalize(_N)
  ((N == 0) | (M == 0)) && return nothing
  _spa, spap = stridedpointer_preserve(A)
  _spc, spcp = stridedpointer_preserve(C)
  _spu, spup = stridedpointer_preserve(U)
  spa = zero_offsets(_spa)
  spc = zero_offsets(_spc)
  spu = zero_offsets(_spu)
  GC.@preserve spap spcp spup begin
    mtb = m_thread_block_size(M, N, Val(T))
    if THREAD && (VectorizationBase.num_threads() > 1)
      (M > mtb) && return multithread_rdiv!(spc, spa, spu, M, N, mtb, Val(UNIT), VectorizationBase.contiguous_axis(A))
    elseif N > block_size(Val(T))
      return rdiv_block_MandN!(spc, spa, spu, M, N, Val(UNIT), VectorizationBase.contiguous_axis(A))
    end
    return rdiv_U!(spc, spa, spu, M, N, VectorizationBase.contiguous_axis(A), Val(UNIT))
  end
end

function rdiv!(A::AbstractMatrix{T}, U::UpperTriangular{T}, ::Val{THREAD} = Val(true)) where {T<:Union{Float32,Float64},THREAD}
  div_dispatch!(A, A, parent(U), Val(false), Val(THREAD))
  return A
end
function rdiv!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, U::UpperTriangular{T}, ::Val{THREAD} = Val(true)) where {T<:Union{Float32,Float64},THREAD}
  div_dispatch!(C, A, parent(U), Val(false), Val(THREAD))
  return C
end
function rdiv!(A::AbstractMatrix{T}, U::UnitUpperTriangular{T}, ::Val{THREAD} = Val(true)) where {T<:Union{Float32,Float64},THREAD}
  div_dispatch!(A, A, parent(U), Val(true), Val(THREAD))
  return A
end
function rdiv!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, U::UnitUpperTriangular{T}, ::Val{THREAD} = Val(true)) where {T<:Union{Float32,Float64},THREAD}
  div_dispatch!(C, A, parent(U), Val(true), Val(THREAD))
  return C
end
function ldiv!(U::LowerTriangular{T}, A::AbstractMatrix{T}, ::Val{THREAD} = Val(true)) where {T<:Union{Float32,Float64},THREAD}
  div_dispatch!(transpose(A), transpose(A), transpose(parent(U)), Val(false), Val(THREAD))
  return A
end
function ldiv!(C::AbstractMatrix{T}, U::LowerTriangular{T}, A::AbstractMatrix{T}, ::Val{THREAD} = Val(true)) where {T<:Union{Float32,Float64},THREAD}
  div_dispatch!(transpose(C), transpose(A), transpose(parent(U)), Val(false), Val(THREAD))
  return C
end
function ldiv!(U::UnitLowerTriangular{T}, A::AbstractMatrix{T}, ::Val{THREAD} = Val(true)) where {T<:Union{Float32,Float64},THREAD}
  div_dispatch!(transpose(A), transpose(A), transpose(parent(U)), Val(true), Val(THREAD))
  return A
end
function ldiv!(C::AbstractMatrix{T}, U::UnitLowerTriangular{T}, A::AbstractMatrix{T}, ::Val{THREAD} = Val(true)) where {T<:Union{Float32,Float64},THREAD}
  div_dispatch!(transpose(C), transpose(A), transpose(parent(U)), Val(true), Val(THREAD))
  return C
end

ldiv!(A, B, ::Val = Val(true)) = LinearAlgebra.ldiv!(A, B)
ldiv!(Y, A, B, ::Val = Val(true)) = LinearAlgebra.ldiv!(Y, A, B)
rdiv!(A, B, ::Val = Val(true)) = LinearAlgebra.rdiv!(A, B)

function block_size(::Val{T}) where {T}
  elements_l2 = (VectorizationBase.cache_size(StaticInt(2))*StaticInt(19)) ÷ (VectorizationBase.static_sizeof(T)*StaticInt(60))
  Static.floortostaticint(sqrt(elements_l2))
end

function nmuladd!(C, A, U, M, K, N)
  @turbo for n ∈ CloseOpen(N), m ∈ CloseOpen(M)
    Cmn = A[m,n]
    for k ∈ CloseOpen(K)
      Cmn -= C[m,k]*U[k,n]
    end
    C[m,K+n] = Cmn
  end
end

function rdiv_block_N!(
  spc::AbstractStridedPointer{T}, spa, spu, M, N, ::Val{UNIT}, ::StaticInt{X}, Bsize = nothing
) where {T,UNIT,X}
  spa_rdiv = spa
  spc_base = spc
  n = 0
  W = VectorizationBase.pick_vector_width(T)
  B_normalized = Bsize === nothing ? VectorizationBase.vcld(N, VectorizationBase.vcld(N, block_size(Val(T)))*W)*W : Bsize
  repeat = N > B_normalized  
  N_temp = Core.ifelse(repeat, B_normalized, N)
  while true
    # println("Solve with N_temp = $N_temp and n = $n")
    rdiv_U!(spc, spa_rdiv, gesp(spu, (n,StaticInt{0}())), M, N_temp, StaticInt{X}(), Val(UNIT))
    repeat || break
    spa = gesp(spa, (StaticInt(0), B_normalized))
    spc = gesp(spc, (StaticInt(0), B_normalized))
    spu = gesp(spu, (StaticInt(0), B_normalized))
    nnext = n + B_normalized
    # N_temp = 
    n += B_normalized
    repeat = n + B_normalized < N
    N_temp = repeat ? N_temp : N - n
    # N_temp = min(n + B_normalized, N) - n
    # println("nmuladd with N_temp = $N_temp and n = $n")
    # mul!(
    #   copyto!(view(C, :, n+1:n+N_temp), view(A, :, n+1:n+N_temp)),
    #   view(C, :, 1:n),
    #   view(U, 1:n, n+1:n+N_temp),
    #   -1.0, 1.0
    # )
    nmuladd!(spc_base, spa, spu, M, n, N_temp)
    spa_rdiv = spc
  end
end
function rdiv_block_MandN!(
  spc::AbstractStridedPointer{T}, spa, spu, M, N, ::Val{UNIT}, ::StaticInt{X}
) where {T,UNIT,X}
  B = block_size(Val(T))
  W = VectorizationBase.pick_vector_width(T)
  B_normalized = VectorizationBase.vcld(N, VectorizationBase.vcld(N, B)*W)*W
  WUF = W*unroll_factor(W)
  B_m = VectorizationBase.vcld(M, VectorizationBase.vcld(M, B)*WUF)*WUF
  m = 0
  while m < M
    mu = m + B_m
    Mtemp = min(M, mu) - m
    rdiv_block_N!(
      spc, spa, spu, Mtemp, N, Val(UNIT), StaticInt{X}(),
      VectorizationBase.vcld(N, VectorizationBase.vcld(N, B)*W)*W
    )
    spa = gesp(spa, (B_m, StaticInt{0}()))
    spc = gesp(spc, (B_m, StaticInt{0}()))
    m = mu
  end
  nothing
end
function _nthreads()
  nc = VectorizationBase.num_cores()
  nt = VectorizationBase.num_threads()
  ifelse(Static.lt(nc,nt),nc,nt)
end
function m_thread_block_size(M, N, ::Val{T}) where {T}
  W = VectorizationBase.pick_vector_width(T)
  WUF = W * unroll_factor(W)
  nb = clamp(VectorizationBase.vdiv(M * N, StaticInt{256}() * W), 1, Int(_nthreads()))
  min(M, VectorizationBase.vcld(M, nb*W)*W)
end

function multithread_rdiv!(
  spc::AbstractStridedPointer{T}, spa, spu, M, N, mtb, ::Val{UNIT}, ::StaticInt{X}
) where {X,T,UNIT}
  mtb = 8
  (Md, Mr) = VectorizationBase.vdivrem(M, mtb)
  Nblock = Md + (Mr ≠ 0)
  Mrem = Core.ifelse(Mr ≠ 0, Mr, mtb)
  # @show mtb, Nblock, Mrem, Md, Mr
  # return
  let Md = Md, Mr = Mr, Nblock = Md + (Mr ≠ 0), Mrem = Core.ifelse(Mr ≠ 0, Mr, mtb), VUNIT = Val{UNIT}(), StaticX = StaticInt{X}()
    @batch for block in CloseOpen(Nblock)
      # for block in CloseOpen(Nblock)
      # let block = 0
      rdiv_block_MandN!(
        # rdiv_block_N!(
        gesp(spc, (mtb*block, StaticInt{0}())),
        gesp(spa, (mtb*block, StaticInt{0}())),
        spu, Core.ifelse(block == Nblock-1, Mrem, mtb), N, VUNIT, StaticX
        # spu, M, N, Val{UNIT}(), StaticInt{X}()
      )
    end
  end
  nothing
  # nlaunch = Md - (Mr == 0)
  # threads, torelease = Polyester.request_threads(Base.Threads.threadid(), nlaunch)
  # nthread = length(threads)
  # if (nthread % Int32) ≤ zero(Int32)
  #   return rdiv_block_MandN!(spc, spa, spu, M, N, Val(UNIT), StaticInt{X}())
  # end
  # nbatch = nthread + one(nthread)
  
end

# We're using `W x W` blocks, consuming `W` registers
# For each block we need to load 1 more value, plus another register is used for `B`. So:
# remaining_registers == register_count() - num_blocks * (W + 1) - 1
# 0 < register_count() - num_blocks * (W + 1) - 1
# num_blocks < (register_count() - 1) / (W + 1)
# num_blocks = (register_count() - 1) ÷ (W + 1)
function unroll_factor(::StaticInt{W}) where {W}
  num_blocks = (VectorizationBase.register_count() - StaticInt{1}()) ÷ (StaticInt{W}() + StaticInt{1}())
  ifelse(Static.lt(num_blocks, StaticInt{1}()), StaticInt{1}(), num_blocks)
end

function rdiv_U!(spc::AbstractStridedPointer{T}, spa::AbstractStridedPointer, spu::AbstractStridedPointer, M, N, ::StaticInt{var"#UNUSED#"}, ::Val{UNIT}) where {T,UNIT,var"#UNUSED#"}
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF*WS
  MU = UF > 1 ? M : 0
  Nd, Nr = VectorizationBase.vdivrem(N, WS)
  spb = lubuffer(Val(T), UF, N)
  m = 0
  while m < MU - WU + 1
    n = Nr
    if n > 0
      BdivU_small_kern_u!(spb, spc, spa, spu, n, UF, Val(UNIT))
    end
    for i ∈ 1:Nd
      rdiv_solve_W_u!(spb, spc, spa, spu, n, WS, UF, Val(UNIT))
      n += W
    end
    m += WU
    spa = gesp(spa, (WU,StaticInt(0)))
    spc = gesp(spc, (WU,StaticInt(0)))
  end
  finalmask = VectorizationBase.mask(WS, M)
  while m < M
    ubm = m+W
    nomaskiter = ubm < M
    mask = nomaskiter ? VectorizationBase.max_mask(WS) : finalmask
    n = Nr
    if n > 0
      BdivU_small_kern!(spb, spc, spa, spu, n, mask, Val(UNIT))
    end
    for i ∈ 1:Nd
      # @show C, n
      rdiv_solve_W!(spb, spc, spa, spu, n, i ≠ Nd, mask, Val(UNIT))
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

@static if VERSION >= v"1.8.0-beta1"
  let
    while true
      A = rand(1, 1)
      B = rand(1, 1)
      res = similar(A)
      rdiv!(res, A, UpperTriangular(B))
      rdiv!(res, A, UnitUpperTriangular(B))
      rdiv!(res, A, UpperTriangular(B), Val(false))
      rdiv!(res, A, UnitUpperTriangular(B), Val(false))

      __init__()
      ldiv!(res, LowerTriangular(B), A)
      ldiv!(res, UnitLowerTriangular(B), A)
      ldiv!(res, LowerTriangular(B), A, Val(false))
      ldiv!(res, UnitLowerTriangular(B), A, Val(false))
      break
    end
  end
end
end
