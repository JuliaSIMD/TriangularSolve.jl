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

@inline maybestore!(p, v, i) = vstore!(p, v, i)
@inline maybestore!(::Nothing, v, i) = nothing

@inline maybestore!(p, v, i, m) = vstore!(p, v, i, m)
@inline maybestore!(::Nothing, v, i, m) = nothing

@inline function store_small_kern!(spa, sp, v, _, i, n, mask, ::Val{true})
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
  L > length(buff) && resize!(buff, L%UInt)
  ptr = Base.unsafe_convert(Ptr{T}, buff)
  si = StrideIndex{2,(1,2),1}((VectorizationBase.static_sizeof(T), RSUF), (StaticInt(0),StaticInt(0)))
  stridedpointer(ptr, si, StaticInt{0}())
end
_canonicalize(x) = signed(x)
_canonicalize(::StaticInt{N}) where {N} = StaticInt{N}()
function div_dispatch!(C::AbstractMatrix{T}, A, U, nthread, ::Val{UNIT}) where {UNIT,T}
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
    mtb = m_thread_block_size(M, N, nthread, Val(T))
    if nthread > 1
      (M > mtb) && return multithread_rdiv!(spc, spa, spu, M, N, mtb, Val(UNIT), VectorizationBase.contiguous_axis(A))
    elseif N > block_size(Val(T))
      return rdiv_block_MandN!(spc, spa, spu, M, N, Val(UNIT), VectorizationBase.contiguous_axis(A))
    end
    return rdiv_U!(spc, spa, spu, M, N, VectorizationBase.contiguous_axis(A), Val(UNIT))
  end
end

_nthreads() = min(Int(VectorizationBase.num_cores())::Int, Threads.nthreads()::Int)
function rdiv!(A::AbstractMatrix{T}, U::UpperTriangular{T}, ::Val{true} = Val(true)) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, A, parent(U), _nthreads(), Val(false))
  return A
end
function rdiv!(A::AbstractMatrix{T}, U::UpperTriangular{T}, ::Val{false}) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, A, parent(U), static(0), Val(false))
  return A
end
function rdiv!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, U::UpperTriangular{T}, ::Val{true} = Val(true)) where {T<:Union{Float32,Float64}}
  div_dispatch!(C, A, parent(U), _nthreads(), Val(false))
  return C
end
function rdiv!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, U::UpperTriangular{T}, ::Val{false}) where {T<:Union{Float32,Float64}}
  div_dispatch!(C, A, parent(U), static(0), Val(false))
  return C
end
function rdiv!(A::AbstractMatrix{T}, U::UnitUpperTriangular{T}, ::Val{true} = Val(true)) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, A, parent(U), _nthreads(), Val(true))
  return A
end
function rdiv!(A::AbstractMatrix{T}, U::UnitUpperTriangular{T}, ::Val{false}) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, A, parent(U), static(0), Val(true))
  return A
end
function rdiv!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, U::UnitUpperTriangular{T}, ::Val{true} = Val(true)) where {T<:Union{Float32,Float64}}
  div_dispatch!(C, A, parent(U), _nthreads(), Val(true))
  return C
end
function rdiv!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, U::UnitUpperTriangular{T}, ::Val{false}) where {T<:Union{Float32,Float64}}
  div_dispatch!(C, A, parent(U), static(0), Val(true))
  return C
end
function ldiv!(U::LowerTriangular{T}, A::AbstractMatrix{T}, ::Val{true} = Val(true)) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(A), transpose(A), transpose(parent(U)), _nthreads(), Val(false))
  return A
end
function ldiv!(U::LowerTriangular{T}, A::AbstractMatrix{T}, ::Val{false}) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(A), transpose(A), transpose(parent(U)), static(0), Val(false))
  return A
end
function ldiv!(C::AbstractMatrix{T}, U::LowerTriangular{T}, A::AbstractMatrix{T}, ::Val{true} = Val(true)) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(C), transpose(A), transpose(parent(U)), _nthreads(), Val(false))
  return C
end
function ldiv!(C::AbstractMatrix{T}, U::LowerTriangular{T}, A::AbstractMatrix{T}, ::Val{false}) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(C), transpose(A), transpose(parent(U)), static(0), Val(false))
  return C
end
function ldiv!(U::UnitLowerTriangular{T}, A::AbstractMatrix{T}, ::Val{true} = Val(true)) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(A), transpose(A), transpose(parent(U)), _nthreads(), Val(true))
  return A
end
function ldiv!(U::UnitLowerTriangular{T}, A::AbstractMatrix{T}, ::Val{false}) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(A), transpose(A), transpose(parent(U)), static(0), Val(true))
  return A
end
function ldiv!(C::AbstractMatrix{T}, U::UnitLowerTriangular{T}, A::AbstractMatrix{T}, ::Val{true} = Val(true)) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(C), transpose(A), transpose(parent(U)), _nthreads(), Val(true))
  return C
end
function ldiv!(C::AbstractMatrix{T}, U::UnitLowerTriangular{T}, A::AbstractMatrix{T}, ::Val{false}) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(C), transpose(A), transpose(parent(U)), static(0), Val(true))
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
    rdiv_U!(spc, spa_rdiv, gesp(spu, (n,StaticInt{0}())), M, N_temp, StaticInt{X}(), Val{UNIT}())
    repeat || break
    spa = gesp(spa, (StaticInt(0), B_normalized))
    spc = gesp(spc, (StaticInt(0), B_normalized))
    spu = gesp(spu, (StaticInt(0), B_normalized))
    n += B_normalized
    repeat = n + B_normalized < N
    N_temp = repeat ? N_temp : N - n
    nmuladd!(spc_base, spa, spu, M, n, N_temp)
    spa_rdiv = spc
  end
end
function rdiv_block_MandN!(
  spc::AbstractStridedPointer{T}, spa, spu, M, N, ::Val{UNIT}, ::StaticInt{X}
) where {T,UNIT,X}
  B = block_size(Val(T))
  W = VectorizationBase.pick_vector_width(T)
  WUF = W*unroll_factor(W)
  B_m = VectorizationBase.vcld(M, VectorizationBase.vcld(M, B)*WUF)*WUF
  m = 0
  while m < M
    mu = m + B_m
    Mtemp = min(M, mu) - m
    rdiv_block_N!(
      spc, spa, spu, Mtemp, N, Val{UNIT}(), StaticInt{X}(),
      VectorizationBase.vcld(N, VectorizationBase.vcld(N, B)*W)*W
    )
    spa = gesp(spa, (B_m, StaticInt{0}()))
    spc = gesp(spc, (B_m, StaticInt{0}()))
    m = mu
  end
  nothing
end
function m_thread_block_size(M, N, nthreads, ::Val{T}) where {T}
  W = VectorizationBase.pick_vector_width(T)
  nb = clamp(VectorizationBase.vdiv(M * N, StaticInt{256}() * W), 1, nthreads)
  min(M, VectorizationBase.vcld(M, nb*W)*W)
end

struct RDivBlockMandNv2{UNIT,X} end
function (f::RDivBlockMandNv2{UNIT,X})(allargs, blockstart, blockstop) where {UNIT,X}
  spc, spa, spu, N, Mrem, Nblock, mtb = allargs
  for block = blockstart-1:blockstop-1
    rdiv_block_MandN!(
      gesp(spc, (mtb*block, StaticInt{0}())),
      gesp(spa, (mtb*block, StaticInt{0}())),
      spu, Core.ifelse(block == Nblock-1, Mrem, mtb), N, Val{UNIT}(), static(X)
    )
  end
end


function multithread_rdiv!(
  spc::AbstractStridedPointer{TC}, spa::AbstractStridedPointer{TA}, spu::AbstractStridedPointer{TU}, M::Int, N::Int, mtb::Int, ::Val{UNIT}, ::StaticInt{X}
) where {X,UNIT,TC,TA,TU}
  # Main._a[] = (spc, spa, spu, M, N, mtb, Val(UNIT), static(X));
  (Md, Mr) = VectorizationBase.vdivrem(M, mtb)
  Nblock = Md + (Mr ≠ 0)
  Mrem = Core.ifelse(Mr ≠ 0, Mr, mtb)
  f = RDivBlockMandNv2{UNIT,X}()
  batch(f, (Nblock,min(Nblock,Threads.nthreads())), spc, spa, spu, N, Mrem, Nblock, mtb)
  nothing
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
    for _ ∈ 1:Nd
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

using SnoopPrecompile
@static if VERSION >= v"1.8.0-beta1"
  @precompile_setup begin
    A = rand(1, 1)
    B = rand(1, 1)
    res = similar(A)
    @precompile_all_calls begin
      rdiv!(res, A, UpperTriangular(B))
      rdiv!(res, A, UnitUpperTriangular(B))
      rdiv!(res, A, UpperTriangular(B), Val(false))
      rdiv!(res, A, UnitUpperTriangular(B), Val(false))

      __init__()
      ldiv!(res, LowerTriangular(B), A)
      ldiv!(res, UnitLowerTriangular(B), A)
      ldiv!(res, LowerTriangular(B), A, Val(false))
      ldiv!(res, UnitLowerTriangular(B), A, Val(false))
    end
  end
end
