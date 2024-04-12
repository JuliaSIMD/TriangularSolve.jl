module TriangularSolve
using Base: @nexprs
if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
  @eval Base.Experimental.@max_methods 1
end

using LayoutPointers: stridedpointer_preserve, StrideIndex
using VectorizationBase, LinearAlgebra #LoopVectorization
using VectorizationBase:
  vfnmadd_fast,
  AbstractStridedPointer,
  AbstractMask,
  zero_offsets,
  gesp,
  StridedPointer
using CloseOpenIntervals: CloseOpen, SafeCloseOpen
using Static
using IfElse: ifelse
using LoopVectorization
using Polyester

@generated function solve_AU(
  A::VecUnroll{Nm1},
  spu::AbstractStridedPointer,
  noff,
  ::Val{UNIT}
) where {Nm1,UNIT}
  A_n_expr = UNIT ? :nothing : :(A_n = Base.FastMath.div_fast(A_n, U_n_n))
  N = Nm1 + 1
  quote
    $(Expr(:meta, :inline))
    Ad = VectorizationBase.data(A)
    Base.Cartesian.@nexprs $N n -> begin
      A_n = Ad[n]
      Base.Cartesian.@nexprs $(UNIT ? :(n - 1) : :n) m -> begin
        U_m_n = vload(spu, (noff + (m - 1), noff + (n - 1)))
      end
    end
    Base.Cartesian.@nexprs $N n -> begin
      Base.Cartesian.@nexprs n - 1 k -> begin
        A_n = vfnmadd_fast(A_k, U_k_n, A_n)
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
@inline store_small_kern!(spa, ::Nothing, v, spu, i, n, mask, ::Val{true}) =
  vstore!(spa, v, i, mask)

@inline function store_small_kern!(spa, sp, v, spu, i, n, mask, ::Val{false})
  x = v / vload(spu, (n, n))
  vstore!(spa, x, i, mask)
  vstore!(sp, x, i, mask)
end
@inline store_small_kern!(spa, ::Nothing, v, spu, i, n, mask, ::Val{false}) =
  vstore!(spa, v / vload(spu, (n, n)), i, mask)

@inline function store_small_kern!(spa, sp, v, spu, i, n, ::Val{true})
  vstore!(spa, v, i)
  vstore!(sp, v, i)
end
@inline store_small_kern!(spa, ::Nothing, v, spu, i, n, ::Val{true}) =
  vstore!(spa, v, i)

@inline function store_small_kern!(spa, sp, v, spu, i, n, ::Val{false})
  x = v / vload(spu, (n, n))
  vstore!(spa, x, i)
  vstore!(sp, x, i)
end
@inline store_small_kern!(spa, ::Nothing, v, spu, i, n, ::Val{false}) =
  vstore!(spa, v / vload(spu, (n, n)), i)

@inline function BdivU_small_kern!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spu::AbstractStridedPointer{T},
  N,
  mask::AbstractMask{W},
  ::Val{UNIT}
) where {T,UNIT,W}
  # W = VectorizationBase.pick_vector_width(T)
  for n ∈ CloseOpen(N)
    Amn = vload(spb, (MM{W}(StaticInt(0)), n), mask)
    for k ∈ SafeCloseOpen(n)
      Amn = vfnmadd_fast(
        vload(spa, (MM{W}(StaticInt(0)), k), mask),
        vload(spu, (k, n)),
        Amn
      )
    end
    store_small_kern!(
      spa,
      sp,
      Amn,
      spu,
      (MM{W}(StaticInt(0)), n),
      n,
      mask,
      Val{UNIT}()
    )
  end
end
@inline function BdivU_small_kern_u!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spu::AbstractStridedPointer{T},
  N,
  ::StaticInt{U},
  ::Val{UNIT}
) where {T,U,UNIT}
  W = Int(VectorizationBase.pick_vector_width(T))
  for n ∈ CloseOpen(N)
    Amn = vload(spb, Unroll{1,W,U,1,W,zero(UInt),1}((StaticInt(0), n)))
    for k ∈ SafeCloseOpen(n)
      Amk = vload(spa, Unroll{1,W,U,1,W,zero(UInt),1}((StaticInt(0), k)))
      Amn = vfnmadd_fast(Amk, vload(spu, (k, n)), Amn)
    end
    store_small_kern!(
      spa,
      sp,
      Amn,
      spu,
      Unroll{1,W,U,1,W,zero(UInt),1}((StaticInt(0), n)),
      n,
      Val{UNIT}()
    )
  end
end

@generated function rdiv_solve_W_u!(
  spc,
  spb,
  spa,
  spu,
  n,
  ::StaticInt{W},
  ::StaticInt{U},
  ::Val{UNIT}
) where {W,U,UNIT}
  z = static(0)
  quote
    $(Expr(:meta, :inline))
    # C = A / U; C * U = A
    # A_{i,j} = C_{i,j}U_{j,j} + \sum_{k=1}^{j-1} C_{i,k}U_{k,j}
    # C_{i,j} = (A_{i,j} - \sum_{k=1}^{j-1} C_{i,k}U_{k,j}) / U_{j,j}
    # Load A_{i,j}
    # Actually: (A_{i+[0,W*U), j+[0,W)}):
    # outer unroll are `W` columns
    # Inner unroll are `W*U` rows (U simd vecs)
    C11 = VectorizationBase.data(
      vload(
        spa,
        $(Unroll{2,1,W,1,W,zero(UInt),1})(
          $(Unroll{1,W,U,1,W,zero(UInt),1})(($z, n))
        )
      )
    )
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    for nk ∈ SafeCloseOpen(n) # nmuladd
      A11 = vload(spc, $(Unroll{1,W,U,1,W,zero(UInt),1})(($(StaticInt(0)), nk)))
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(spu, (nk, n + (c - 1))), C11_c)
    end
    C11vu =
      solve_AU(VecUnroll((Base.Cartesian.@ntuple $W C11)), spu, n, $(Val(UNIT)))
    i = $(Unroll{2,1,W,1,W,zero(UInt),1})(
      $(Unroll{1,W,U,1,W,zero(UInt),1})(($z, n))
    )
    vstore!(spc, C11vu, i)
    maybestore!(spb, C11vu, i)
  end
end
@generated function rdiv_solve_W!(
  spc,
  spb,
  spa,
  spu,
  n,
  storec::B,
  mask::AbstractMask{W},
  ::Val{UNIT}
) where {W,UNIT,B}
  storecexpr = if (B <: Bool)
    :(storec && vstore!(spc, C11, i, mask))
  else
    :(vstore!(spc, C11, i, mask))
  end
  z = static(0)
  quote
    $(Expr(:meta, :inline))
    # here, we just want to load the vectors
    C11 = VectorizationBase.data(
      vload(spa, $(Unroll{2,1,W,1,W,(-1 % UInt),1})(($z, n)), mask)
    )
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    for nk ∈ SafeCloseOpen(n) # nmuladd
      A11 = vload(spc, ($(MM{W}(z)), nk), mask)
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(spu, (nk, n + (c - 1))), C11_c)
    end
    C11 = VecUnroll((Base.Cartesian.@ntuple $W C11))
    C11 = solve_AU(C11, spu, n, $(Val(UNIT)))
    i = $(Unroll{2,1,W,1,W,(-1 % UInt),1})(($z, n))
    $storecexpr
    maybestore!(spb, C11, i, mask)
  end
end

@generated function ldiv_solve_W_u!(
  spc,
  spa,
  spu,
  n,
  ::StaticInt{W},
  ::StaticInt{U},
  ::Val{UNIT}
) where {W,U,UNIT}
  z = static(0)
  quote
    $(Expr(:meta, :inline))
    # C = L \ A; L * C = A
    # A_{i,j} = L_{i,i}*C_{i,j} + \sum_{k=1}^{i-1}L_{i,k}C_{k,j}
    # C_{i,j} = L_{i,i} \ (A_{i,j} - \sum_{k=1}^{i-1}L_{i,k}C_{k,j})
    # The inputs here are transposed, as the library was formulated in terms of `rdiv!`,
    # so we have
    # C_{j,i} = (A_{j,i} - \sum_{k=1}^{i-1}C_{j,k}U_{k,i}) / L_{i,i}
    # This solves for the block: C_{j+[0,W],i+[0,W*U)}
    # This can be viewed as `U` blocks that are each `W`x`W`
    # E.g. U=3, rough alg:
    # r=[0,W); c=[0,WU)
    # X = A_{j+r,i+c} - \sum_{k=1}^{i-1}C_{j+r,k}*U_{k,i+c}
    # C_{j+r,i+r} =  X[:, r] / U_{i+r,i+r}
    # C_{j+r,i+W+r} = (X[:, W+r] - C_{j+r,i+r}*U_{i+r,i+W+r}) / U_{i+W+r,i+W+r}
    # C_{j+r,i+2W+r} = (X[:, 2W+r] - C_{j+r,i+r}*U_{i+r,i+2W+r} - C_{j+r,i+W+r}*U_{i+W+r,i+2W+r}) / U_{i+2W+r,i+2W+r}
    #
    # outer unroll are `W` rows
    # Inner unroll are `W*U` columns (U simd vecs)
    # 
    A11 = getfield(
      vload(
        spa,
        $(Unroll{1,1,W,2,W,zero(UInt),1})(
          $(Unroll{2,W,U,2,W,zero(UInt),1})(($z, n))
        )
      ),
      :data
    )
    # The `W` rows
    Base.Cartesian.@nexprs $W c -> A11_c = getfield(A11, c)
    # compute
    # A_{j,i} - \sum_{k=1}^{i-1}U_{k,i}C_{j,k})
    # Each iter:
    # A_{j+[0,W), i+[0,W*U)} -= C_{j+[0,W),k}*U_{k,i+[0,W*U)}
    for nk ∈ SafeCloseOpen(n) # nmuladd
      U_ki = vload(spu, $(Unroll{2,W,U,2,W,zero(UInt),1})((nk, n)))
      Base.Cartesian.@nexprs $W c ->
        A11_c = vfnmadd_fast(U_ki, vload(spc, (static(c - 1), nk)), A11_c)
    end
    # solve AU wants:
    # outer unroll are `W` columns
    # Inner unroll are `W` rows (U simd vecs)
    # So, we'll use `U = 1`, and transpose blocks
    # We then have column-major multiplies
    Base.Cartesian.@nexprs $U u -> begin
      # take A[(u-1)*W,u*W), [0,W)]
      X_u = getfield(
        VectorizationBase.transpose_vecunroll(
          VecUnroll(
            Base.Cartesian.@ntuple $W w ->
              getfield(getfield(A11_w, :data), u)
          )
        ),
        :data
      )
      Base.Cartesian.@nexprs $W c -> X_u_c = getfield(X_u, c)
      Base.Cartesian.@nexprs (u - 1) j -> begin
        # subtract
        # r = W*(j-1)+[0,W)
        # A_{j+[0,W),i+r} -= C_{j+[0,W),r}*U_{r,i+r}
        # W x W matmul
        Base.Cartesian.@nexprs $W k -> begin # reduction
          Base.Cartesian.@nexprs $W c -> begin # cols
            U_u_j_k_c = vload(
              spu,
              (n + ((k - 1) + ((j - 1) * $W)), n + ((c - 1) + ((u - 1) * $W)))
            )
            X_u_c = vfnmadd_fast(C_j_k, U_u_j_k_c, X_u_c)
          end
        end
      end
      C_u = solve_AU(
        VecUnroll(Base.Cartesian.@ntuple $W X_u),
        spu,
        n + ((u - 1) * $W),
        $(Val(UNIT))
      )
      Cdata_u = getfield(C_u, :data)
      Base.Cartesian.@nexprs $W c -> C_u_c = getfield(Cdata_u, c)
    end
    # store at end (no aliasing)
    Base.Cartesian.@nexprs $U u -> begin
      vstore!(spc, C_u, $(Unroll{2,1,W,1,W,zero(UInt),1})(($z, n + (u - 1) * $W)))
    end
  end
end
@generated function ldiv_solve_W!(
  spc,
  spa,
  spu,
  n,
  ::StaticInt{W},
  ::Val{UNIT}
) where {W,UNIT}
  z = static(0)
  quote
    $(Expr(:meta, :inline))
    # Like `ldiv_solve_W_u!`, except no unrolling, just a `W`x`W` block
    #
    # C = L \ A; L * C = A
    # A_{i,j} = L_{i,i}*C_{i,j} + \sum_{k=1}^{i-1}L_{i,k}C_{k,j}
    # C_{i,j} = L_{i,i} \ (A_{i,j} - \sum_{k=1}^{i-1}L_{i,k}C_{k,j})
    # The inputs here are transposed, as the library was formulated in terms of `rdiv!`,
    # so we have
    # C_{j,i} = (A_{j,i} - \sum_{k=1}^{i-1}C_{j,k}U_{k,i}) / L_{i,i}
    # This solves for the block: C_{j+[0,W],i+[0,W)}
    # Rough alg:
    # r=[0,W);
    # X = A_{j+r,i+r} - \sum_{k=1}^{i-1}C_{j+r,k}*U_{k,i+r}
    # C_{j+r,i+r} =  X / U_{i+r,i+r}
    #
    # Load the `W`x`W` block...
    # what about masking?
    A11 =
      getfield(vload(spa, $(Unroll{1,1,W,2,W,zero(UInt),1})(($z, n))), :data)
    # The `W` rows
    Base.Cartesian.@nexprs $W c -> A11_c = getfield(A11, c)
    # compute
    # A_{j,i} - \sum_{k=1}^{i-1}U_{k,i}C_{j,k})
    # Each iter:
    # A_{j+[0,W), i+[0,W*U)} -= C_{j+[0,W),k}*U_{k,i+[0,W*U)}
    for nk ∈ SafeCloseOpen(n) # nmuladd
      U_ki = vload(spu, (nk, $(MM{W})(n)))
      Base.Cartesian.@nexprs $W c ->
        A11_c = vfnmadd_fast(U_ki, vload(spc, (static(c - 1), nk)), A11_c)
    end
    # Base.Cartesian.@nexprs $W c -> @show A11_c
    # solve AU wants us to transpose
    # We then have column-major multiplies
    # take A[(u-1)*W,u*W), [0,W)]
    X = VectorizationBase.transpose_vecunroll(
      VecUnroll(Base.Cartesian.@ntuple $W A11)
    )
    # @show X
    C_u = solve_AU(X, spu, n, $(Val(UNIT)))
    vstore!(spc, C_u, $(Unroll{2,1,W,1,W,zero(UInt),1})(($z, n)))
  end
end
@generated function ldiv_solve_W!(
  spc,
  spa,
  spu,
  n,
  ::StaticInt{W},
  ::Val{UNIT},
  ::StaticInt{R}
) where {W,UNIT,R}
  R <= 1 && throw("Remainder of `<= 1` shouldn't be called, but had $R.")
  R >= W && throw("Reaminderof `>= $W` shouldn't be called, but had $R.")
  z = static(0)
  WS = static(W)
  q = quote
    $(Expr(:meta, :inline))
    # Like `ldiv_solve_W_u!`, except no unrolling, just a `W`x`W` block
    #
    # C = L \ A; L * C = A
    # A_{i,j} = L_{i,i}*C_{i,j} + \sum_{k=1}^{i-1}L_{i,k}C_{k,j}
    # C_{i,j} = L_{i,i} \ (A_{i,j} - \sum_{k=1}^{i-1}L_{i,k}C_{k,j})
    # The inputs here are transposed, as the library was formulated in terms of `rdiv!`,
    # so we have
    # C_{j,i} = (A_{j,i} - \sum_{k=1}^{i-1}C_{j,k}U_{k,i}) / L_{i,i}
    # This solves for the block: C_{j+[0,R],i+[0,W)}
    # Rough alg:
    # r=[0,R); w=[0,W);
    # X = A_{j+r,i+w} - \sum_{k=1}^{i-1}C_{j+r,k}*U_{k,i+w}
    # C_{j+r,i+w} =  X / U_{i+r,i+w}
    #
    # Load the `W`x`W` block...
    # what about masking?
    A11 =
      getfield(vload(spa, $(Unroll{1,1,R,2,W,zero(UInt),1})(($z, n))), :data)
    # The `W` rows
    Base.Cartesian.@nexprs $R r -> A11_r = getfield(A11, r)
    # compute
    # A_{j,i} - \sum_{k=1}^{i-1}U_{k,i}C_{j,k})
    # Each iter:
    # A_{j+[0,W), i+[0,W*U)} -= C_{j+[0,W),k}*U_{k,i+[0,W*U)}
    for nk ∈ SafeCloseOpen(n) # nmuladd
      U_ki = vload(spu, (nk, $(MM{W})(n)))
      Base.Cartesian.@nexprs $R r ->
        A11_r = vfnmadd_fast(U_ki, vload(spc, (static(r - 1), nk)), A11_r)
    end
  end
  # pad with zeros
  Wpad = VectorizationBase.nextpow2(R)
  t = Expr(:tuple)
  for r = 1:R
    push!(t.args, Symbol(:A11_, r))
  end
  for _ = R+1:Wpad
    push!(t.args, :(zero(A11_1)))
  end
  q2 = quote
    # solve AU wants us to transpose
    # We then have column-major multiplies
    # take A[(u-1)*W,u*W), [0,W)]
    X = VectorizationBase.transpose_vecunroll(VecUnroll($t))
    C_u = solve_AU(X, spu, n, $(Val(UNIT)))
  end
  push!(q.args, q2)
  q3 = if R == Wpad
    quote
      i = $(Unroll{2,1,W,1,Wpad,zero(UInt),1})(($z, n))
      vstore!(spc, C_u, i)
    end
  else
    quote
      mask = VectorizationBase.mask($(static(Wpad)), $(static(R)))
      i = $(Unroll{2,1,W,1,Wpad,(-1 % UInt),1})(($z, n))
      vstore!(spc, C_u, i, mask)
    end
  end
  push!(q.args, q3)
  return q
end
@inline function rdiv_U!(
  spc::AbstractStridedPointer{T},
  spa::AbstractStridedPointer,
  spu::AbstractStridedPointer,
  M,
  N,
  ::StaticInt{1},
  ::StaticInt{1},
  ::Val{UNIT}
) where {T,UNIT}
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF * WS
  Nd, Nr = VectorizationBase.vdivrem(N, WS)
  m = 0
  if UF > 1
    while m < M - WU + 1
      n = Nr
      if n > 0
        BdivU_small_kern_u!(spc, nothing, spa, spu, n, UF, Val(UNIT))
      end
      for _ ∈ 1:Nd
        rdiv_solve_W_u!(spc, nothing, spa, spu, n, WS, UF, Val(UNIT))
        n += W
      end
      m += WU
      spa = gesp(spa, (WU, StaticInt(0)))
      spc = gesp(spc, (WU, StaticInt(0)))
    end
  end
  finalmask = VectorizationBase.mask(WS, M)
  while m < M
    ubm = m + W
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
    spa = gesp(spa, (WS, StaticInt(0)))
    spc = gesp(spc, (WS, StaticInt(0)))
    m = ubm
  end
  nothing
end

const buffer = Ref{Ptr{Cvoid}}(C_NULL)

function __init__()
  bp_size = 2 * sizeof(Int) * Threads.nthreads()
  buffer[] = bp = Libc.malloc(bp_size % UInt)
  Libc.memset(bp, 0, bp_size)
end

function _get_buffer_pointer(::StaticInt{UF}, N) where {UF}
  RS = VectorizationBase.register_size()
  RSUF = StaticInt{UF}() * RS
  L = RSUF * N
  tid = Threads.threadid() - 1
  bp = Ptr{Pair{Ptr{Cvoid},Int}}(buffer[]) + 2sizeof(Int) * tid
  (p, buff_current) = unsafe_load(bp)
  if buff_current < L
    p == C_NULL || Libc.free(p)
    buff_size = max(RSUF * 128, L)
    p = Libc.malloc((buff_size + RS - 1) % UInt)
    unsafe_store!(bp, p => buff_size)
  end
  return VectorizationBase.align(p, RS)
end

@inline function lubuffer(::Val{T}, ::StaticInt{UF}, N) where {T,UF}
  RS = VectorizationBase.register_size()
  RSUF = StaticInt{UF}() * RS
  ptr = Ptr{T}(_get_buffer_pointer(StaticInt{UF}(), N))
  si = StrideIndex{2,(1, 2),1}(
    (VectorizationBase.static_sizeof(T), RSUF),
    (StaticInt(0), StaticInt(0))
  )
  stridedpointer(ptr, si, StaticInt{0}()), nothing
end
@inline function lubuffer(
  ::Val{T},
  ::StaticInt{UF},
  ::StaticInt{N}
) where {T,UF,N}
  RSUF = StaticInt{UF}() * VectorizationBase.pick_vector_width(T)
  L = RSUF * N
  buf = Ref{NTuple{L,T}}()
  ptr = Base.unsafe_convert(Ptr{T}, buf)
  si = StrideIndex{2,(1, 2),1}(
    (
      VectorizationBase.static_sizeof(T),
      RSUF * VectorizationBase.static_sizeof(T)
    ),
    (StaticInt(0), StaticInt(0))
  )
  stridedpointer(ptr, si, StaticInt{0}()), buf
end
@inline _free(p::Ptr) = Libc.free(p)
_canonicalize(x) = signed(x)
_canonicalize(::StaticInt{N}) where {N} = StaticInt{N}()
function div_dispatch!(
  C::AbstractMatrix{T},
  A,
  U,
  nthread,
  ::Val{UNIT}
) where {UNIT,T}
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
  XC = VectorizationBase.contiguous_axis(C)
  XA = VectorizationBase.contiguous_axis(A)
  GC.@preserve spap spcp spup begin
    mtb = m_thread_block_size(M, N, nthread, Val(T))
    if nthread > 1
      (M > mtb) &&
        return multithread_rdiv!(spc, spa, spu, M, N, mtb, Val(UNIT), XC, XA)
    elseif N > block_size(Val(T))
      return rdiv_block_MandN!(spc, spa, spu, M, N, Val(UNIT), XC, XA)
    end
    return rdiv_U!(spc, spa, spu, M, N, XC, XA, Val(UNIT))
  end
end

_nthreads() =
  min(Int(VectorizationBase.num_cores())::Int, Threads.nthreads()::Int)
function rdiv!(
  A::AbstractMatrix{T},
  U::UpperTriangular{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, A, parent(U), _nthreads(), Val(false))
  return A
end
function rdiv!(
  A::AbstractMatrix{T},
  U::UpperTriangular{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, A, parent(U), static(1), Val(false))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UpperTriangular{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(C, A, parent(U), _nthreads(), Val(false))
  return C
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UpperTriangular{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(C, A, parent(U), static(1), Val(false))
  return C
end
function rdiv!(
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, A, parent(U), _nthreads(), Val(true))
  return A
end
function rdiv!(
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, A, parent(U), static(1), Val(true))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(C, A, parent(U), _nthreads(), Val(true))
  return C
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(C, A, parent(U), static(1), Val(true))
  return C
end
function ldiv!(
  U::LowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(A),
    transpose(A),
    transpose(parent(U)),
    _nthreads(),
    Val(false)
  )
  return A
end
function ldiv!(
  U::LowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(A),
    transpose(A),
    transpose(parent(U)),
    static(1),
    Val(false)
  )
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::LowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(C),
    transpose(A),
    transpose(parent(U)),
    _nthreads(),
    Val(false)
  )
  return C
end
function ldiv!(
  C::AbstractMatrix{T},
  U::LowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(C),
    transpose(A),
    transpose(parent(U)),
    static(1),
    Val(false)
  )
  return C
end
function ldiv!(
  U::UnitLowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(A),
    transpose(A),
    transpose(parent(U)),
    _nthreads(),
    Val(true)
  )
  return A
end
function ldiv!(
  U::UnitLowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(A),
    transpose(A),
    transpose(parent(U)),
    static(1),
    Val(true)
  )
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::UnitLowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(C),
    transpose(A),
    transpose(parent(U)),
    _nthreads(),
    Val(true)
  )
  return C
end
function ldiv!(
  C::AbstractMatrix{T},
  U::UnitLowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(C),
    transpose(A),
    transpose(parent(U)),
    static(1),
    Val(true)
  )
  return C
end

ldiv!(A, B, ::Val = Val(true)) = LinearAlgebra.ldiv!(A, B)
ldiv!(Y, A, B, ::Val = Val(true)) = LinearAlgebra.ldiv!(Y, A, B)
rdiv!(A, B, ::Val = Val(true)) = LinearAlgebra.rdiv!(A, B)

function block_size(::Val{T}) where {T}
  elements_l2 =
    (VectorizationBase.cache_size(StaticInt(2)) * StaticInt(19)) ÷
    (VectorizationBase.static_sizeof(T) * StaticInt(60))
  Static.floortostaticint(sqrt(elements_l2))
end

nmuladd!(C, A, U, M, K, N) = @turbo for n ∈ CloseOpen(N), m ∈ CloseOpen(M)
  Cmn = A[m, n]
  for k ∈ CloseOpen(K)
    Cmn -= C[m, k] * U[k, n]
  end
  C[m, K+n] = Cmn
end

function rdiv_block_N!(
  spc::AbstractStridedPointer{T},
  spa,
  spu,
  M,
  N,
  ::Val{UNIT},
  ::StaticInt{XC},
  ::StaticInt{XA},
  Bsize = nothing
) where {T,UNIT,XC,XA}
  spa_rdiv = spa
  spc_base = spc
  n = 0
  W = VectorizationBase.pick_vector_width(T)
  B_normalized =
    Bsize === nothing ?
    VectorizationBase.vcld(
      N,
      VectorizationBase.vcld(N, block_size(Val(T))) * W
    ) * W : Bsize
  repeat = N > B_normalized
  N_temp = Core.ifelse(repeat, B_normalized, N)
  while true
    # println("Solve with N_temp = $N_temp and n = $n")
    rdiv_U!(
      spc,
      spa_rdiv,
      gesp(spu, (n, StaticInt{0}())),
      M,
      N_temp,
      StaticInt{XC}(),
      StaticInt{XA}(),
      Val{UNIT}()
    )
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
  spc::AbstractStridedPointer{T},
  spa,
  spu,
  M,
  N,
  ::Val{UNIT},
  ::StaticInt{XC},
  ::StaticInt{XA}
) where {T,UNIT,XC,XA}
  B = block_size(Val(T))
  W = VectorizationBase.pick_vector_width(T)
  WUF = XC == XA == 2 ? W : W * unroll_factor(W)
  B_m = VectorizationBase.vcld(M, VectorizationBase.vcld(M, B) * WUF) * WUF
  m = 0
  while m < M
    mu = m + B_m
    Mtemp = min(M, mu) - m
    rdiv_block_N!(
      spc,
      spa,
      spu,
      Mtemp,
      N,
      Val{UNIT}(),
      StaticInt{XC}(),
      StaticInt{XA}(),
      VectorizationBase.vcld(N, VectorizationBase.vcld(N, B) * W) * W
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
  min(M, VectorizationBase.vcld(M, nb * W) * W)
end

struct RDivBlockMandNv2{UNIT,XC,XA} end
function (f::RDivBlockMandNv2{UNIT,XC,XA})(
  allargs,
  blockstart,
  blockstop
) where {UNIT,XC,XA}
  spc, spa, spu, N, Mrem, Nblock, mtb = allargs
  for block = blockstart-1:blockstop-1
    rdiv_block_MandN!(
      gesp(spc, (mtb * block, StaticInt{0}())),
      gesp(spa, (mtb * block, StaticInt{0}())),
      spu,
      Core.ifelse(block == Nblock - 1, Mrem, mtb),
      N,
      Val{UNIT}(),
      static(XC),
      static(XA)
    )
  end
end

function multithread_rdiv!(
  spc::AbstractStridedPointer{TC},
  spa::AbstractStridedPointer{TA},
  spu::AbstractStridedPointer{TU},
  M::Int,
  N::Int,
  mtb::Int,
  ::Val{UNIT},
  ::StaticInt{XC},
  ::StaticInt{XA}
) where {XC,XA,UNIT,TC,TA,TU}
  # Main._a[] = (spc, spa, spu, M, N, mtb, Val(UNIT), static(X));
  (Md, Mr) = VectorizationBase.vdivrem(M, mtb)
  Nblock = Md + (Mr ≠ 0)
  Mrem = Core.ifelse(Mr ≠ 0, Mr, mtb)
  f = RDivBlockMandNv2{UNIT,XC,XA}()
  batch(
    f,
    (Nblock, min(Nblock, Threads.nthreads())),
    spc,
    spa,
    spu,
    N,
    Mrem,
    Nblock,
    mtb
  )
  nothing
end

# We're using `W x W` blocks, consuming `W` registers
# For each block we need to load 1 more value, plus another register is used for `B`. So:
# remaining_registers == register_count() - num_blocks * (W + 1) - 1
# 0 < register_count() - num_blocks * (W + 1) - 1
# num_blocks < (register_count() - 1) / (W + 1)
# num_blocks = (register_count() - 1) ÷ (W + 1)
function unroll_factor(::StaticInt{W}) where {W}
  num_blocks =
    (VectorizationBase.register_count() - StaticInt{1}()) ÷
    (StaticInt{W}() + StaticInt{1}())
  ifelse(Static.lt(num_blocks, StaticInt{1}()), StaticInt{1}(), num_blocks)
end

function rdiv_U!(
  spc::AbstractStridedPointer{T},
  spa::AbstractStridedPointer,
  spu::AbstractStridedPointer,
  M,
  N,
  ::StaticInt{var"#UNUSED1#"},
  ::StaticInt{var"#UNUSED2#"},
  ::Val{UNIT}
) where {T,UNIT,var"#UNUSED1#",var"#UNUSED2#"}
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF * WS
  Nd, Nr = VectorizationBase.vdivrem(N, WS)
  spb, preserve = lubuffer(Val(T), UF, N)
  m = 0
  GC.@preserve preserve begin
    if UF > 1
      while m < M - WU + 1
        n = Nr
        if n > 0
          BdivU_small_kern_u!(spb, spc, spa, spu, n, UF, Val(UNIT))
        end
        for _ ∈ 1:Nd
          rdiv_solve_W_u!(spb, spc, spa, spu, n, WS, UF, Val(UNIT))
          n += W
        end
        m += WU
        spa = gesp(spa, (WU, StaticInt(0)))
        spc = gesp(spc, (WU, StaticInt(0)))
      end
    end
    finalmask = VectorizationBase.mask(WS, M)
    while m < M
      ubm = m + W
      nomaskiter = ubm < M
      mask = nomaskiter ? VectorizationBase.max_mask(WS) : finalmask
      n = Nr
      if n > 0
        BdivU_small_kern!(spb, spc, spa, spu, n, mask, Val(UNIT))
      end
      for i ∈ 1:Nd
        rdiv_solve_W!(spb, spc, spa, spu, n, i ≠ Nd, mask, Val(UNIT))
        n += W
      end
      spa = gesp(spa, (WS, StaticInt(0)))
      spc = gesp(spc, (WS, StaticInt(0)))
      m = ubm
    end
  end
  nothing
end

@generated function ldiv_remainder!(
  spc,
  spa,
  spu,
  M,
  N,
  m,
  Nr,
  ::StaticInt{W},
  ::Val{UNIT},
  ::StaticInt{r}
) where {W,UNIT,r}
  r <= 0 && throw("Remainder of `<= 0` shouldn't be called, but had $r.")
  r >= W && throw("Reaminderof `>= $W` shouldn't be called, but had $r.")
  if r == 1
    z = static(0)
    vlxj = :(vload(spc, ($z, j)))
    if UNIT
      vlxj = :(xj = $vlxj)
    else
      vlxj = quote
        xj = $vlxj / vload(spu, (j, j))
        vstore!(spc, xj, ($z, j))
      end
    end
    quote
      if pointer(spc) != pointer(spa)
        for n = 0:N-1
          vstore!(spc, vload(spa, ($z, n)), ($z, n))
        end
      end
      for j = 0:N-1
        $vlxj
        for i = (j+1):N-1
          xi = vload(spc, ($z, i))
          Uji = vload(spu, (j, i))
          vstore!(spc, xi - xj * Uji, ($z, i))
        end
      end
    end
  else
    WS = static(W)
    quote
      $(Expr(:meta, :inline))
      n = Nr # non factor of W remainder
      if n > 0
        mask = $(VectorizationBase.mask(WS, r))
        BdivU_small_kern!(spc, nothing, spa, spu, n, mask, $(Val(UNIT)))
      end
      # while n < N - $(W * U - 1)
      #   ldiv_solve_W_u!(spc, spa, spu, n, $WS, $US, Val(UNIT), Val(r))
      #   n += $(W * U)
      # end
      while n != N
        ldiv_solve_W!(spc, spa, spu, n, $WS, $(Val(UNIT)), $(StaticInt(r)))
        n += $W
      end
    end
  end
end
@generated function ldiv_remainder!(
  spc,
  spa,
  spu,
  M,
  N,
  m,
  Nr,
  ::StaticInt{W},
  # ::Val{U},
  ::Val{UNIT}
) where {W,UNIT}
  WS = static(W)
  # US = static(U)
  quote
    $(Expr(:meta, :inline))
    Base.Cartesian.@nif $W w -> m == M - w w -> ldiv_remainder!(
      spc,
      spa,
      spu,
      M,
      N,
      m,
      Nr,
      $WS,
      $(Val(UNIT)),
      StaticInt(w)
    )
  end
end

# spc = spa / spu
# spc' = (spu' \ spa')'
# This is ldiv
function rdiv_U!(
  spc::AbstractStridedPointer{T},
  spa::AbstractStridedPointer,
  spu::AbstractStridedPointer,
  M,
  N,
  ::StaticInt{2},
  ::StaticInt{2},
  ::Val{UNIT}
) where {T,UNIT}
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF * WS
  MU = UF > 1 ? M : 0
  Nd, Nr = VectorizationBase.vdivrem(N, WS)
  m = 0
  # @show M,N
  # m, no remainder
  while m < M - WS + 1
    n = Nr # non factor of W remainder
    if n > 0
      BdivU_small_kern_u!(spc, nothing, spa, spu, n, StaticInt(1), Val(UNIT))
    end
    while n < N - (WU - 1)
      ldiv_solve_W_u!(spc, spa, spu, n, WS, UF, Val(UNIT))
      n += WU
    end
    while n != N
      ldiv_solve_W!(spc, spa, spu, n, WS, Val(UNIT))
      n += W
    end
    m += W
    spa = gesp(spa, (W, StaticInt(0)))
    spc = gesp(spc, (W, StaticInt(0)))
  end
  # remainder on `m`
  m < M && ldiv_remainder!(spc, spa, spu, M, N, m, Nr, WS, Val(UNIT))
  # m < M && ldiv_remainder!(spc, spa, spu, M, N, m, Nr, WS, UF, Val(UNIT))
  nothing
end

#=
using PrecompileTools
@static if VERSION >= v"1.8.0-beta1"
  @setup_workload begin
    A = rand(1, 1)
    B = rand(1, 1)
    res = similar(A)
    @compile_workload begin
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
=#
end
