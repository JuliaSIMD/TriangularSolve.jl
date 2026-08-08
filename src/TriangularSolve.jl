module TriangularSolve
using Base: @nexprs, @ntuple
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

@inline function store_small_kern!(spa, sp, v, i, mask)
  vstore!(spa, v, i, mask)
  vstore!(sp, v, i, mask)
end
@inline store_small_kern!(spa, ::Nothing, v, i, mask) = vstore!(spa, v, i, mask)

@inline function store_small_kern!(spa, sp, v, i)
  vstore!(spa, v, i)
  vstore!(sp, v, i)
end
@inline store_small_kern!(spa, ::Nothing, v, i) = vstore!(spa, v, i)

@generated function BdivU_small_kern!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spu::AbstractStridedPointer{T},
  ::StaticInt{N},
  mask::AbstractMask{W},
  ::Val{UNIT}
) where {T,UNIT,W,N}
  z = static(0)
  if N == 1
    i = (MM{W}(z), z)
    Amn = :(vload(spb, $i, mask))
    if !UNIT
      Amn = :($Amn / vload(spu, $((z, z))))
    end
    quote
      $(Expr(:meta, :inline))
      store_small_kern!(spa, sp, $Amn, $i, mask)
    end
  else
    unroll = Unroll{2,1,N,1,W,(-1 % UInt),1}((z, z))
    tostore = :(VecUnroll(Base.Cartesian.@ntuple $N Amn))
    scale = UNIT ? nothing : :(Amn_n /= vload(spu, (n - 1, n - 1)))
    quote
      $(Expr(:meta, :inline))
      Amn = getfield(vload(spb, $unroll, mask), :data)
      Base.Cartesian.@nexprs $N n -> begin
        Amn_n = getfield(Amn, n)
        Base.Cartesian.@nexprs (n - 1) k -> begin
          Amn_n = vfnmadd_fast(Amn_k, vload(spu, (k - 1, n - 1)), Amn_n)
        end
        $scale
      end
      store_small_kern!(spa, sp, $tostore, $unroll, mask)
    end
  end
end
@generated function BdivU_small_kern_u!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spu::AbstractStridedPointer{T},
  ::StaticInt{N},
  ::StaticInt{U},
  ::Val{UNIT},
  ::StaticInt{W}
) where {T,U,UNIT,N,W}
  z = static(0)
  if N == 1
    unroll = Unroll{1,W,U,1,W,zero(UInt),1}((z, z))
    Amn = :(vload(spb, $unroll))
    if !UNIT
      Amn = :($Amn / vload(spu, $((z, z))))
    end
    quote
      $(Expr(:meta, :inline))
      store_small_kern!(spa, sp, $Amn, $unroll)
    end
  else
    double_unroll =
      Unroll{2,1,N,1,W,zero(UInt),1}(Unroll{1,W,U,1,W,zero(UInt),1}((z, z)))
    tostore = :(VecUnroll(Base.Cartesian.@ntuple $N Amn))
    scale = UNIT ? nothing : :(Amn_n /= vload(spu, (n - 1, n - 1)))
    quote
      $(Expr(:meta, :inline))
      Amn = getfield(vload(spb, $double_unroll), :data)
      Base.Cartesian.@nexprs $N n -> begin
        Amn_n = getfield(Amn, n)
        Base.Cartesian.@nexprs (n - 1) k -> begin
          Amn_n = vfnmadd_fast(Amn_k, vload(spu, (k - 1, n - 1)), Amn_n)
        end
        $scale
      end
      store_small_kern!(spa, sp, $tostore, $double_unroll)
    end
  end
end
@generated function BdivU_small_kern!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spu::AbstractStridedPointer{T},
  Nr::Int,
  mask::AbstractMask{W},
  ::Val{UNIT}
) where {T,UNIT,W}
  quote
    # $(Expr(:meta, :inline))
    Base.Cartesian.@nif $(W - 1) n -> n == Nr n ->
      BdivU_small_kern!(spa, sp, spb, spu, static(n), mask, $(Val(UNIT)))
  end
end
@generated function BdivU_small_kern_u!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spu::AbstractStridedPointer{T},
  Nr::Int,
  ::StaticInt{U},
  ::Val{UNIT},
  ::StaticInt{W}
) where {T,U,UNIT,W}
  su = static(U)
  vu = Val(UNIT)
  sw = static(W)
  quote
    # $(Expr(:meta, :inline))
    Base.Cartesian.@nif $(W - 1) n -> n == Nr n ->
      BdivU_small_kern_u!(spa, sp, spb, spu, static(n), $su, $vu, $sw)
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
    # $(Expr(:meta, :inline))
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
    # $(Expr(:meta, :inline))
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
    # solve AU wants us to transpose
    # We then have column-major multiplies
    # take A[(u-1)*W,u*W), [0,W)]
    X = VectorizationBase.transpose_vecunroll(
      VecUnroll(Base.Cartesian.@ntuple $W A11)
    )
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
  q = quote
    # $(Expr(:meta, :inline))
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
        BdivU_small_kern_u!(spc, nothing, spa, spu, n, UF, Val(UNIT), WS)
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
  GC.@preserve spap spcp spup begin
    mtb = m_thread_block_size(M, N, nthread, Val(T))
    if nthread > 1
      (M > mtb) && return multithread_rdiv!(spc, spa, spu, M, N, mtb, Val(UNIT))
    elseif N > block_size(Val(T))
      return rdiv_block_MandN!(spc, spa, spu, M, N, Val(UNIT))
    end
    return rdiv_U!(spc, spa, spu, M, N, Val(UNIT))
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

function rdiv!(
  A::StridedMatrix{T},
  L::LowerTriangular{T,<:StridedMatrix{T}},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(A, A, parent(L), _nthreads(), Val(false))
  return A
end
function rdiv!(
  A::StridedMatrix{T},
  L::LowerTriangular{T,<:StridedMatrix{T}},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(A, A, parent(L), static(1), Val(false))
  return A
end
function rdiv!(
  C::StridedMatrix{T},
  A::StridedMatrix{T},
  L::LowerTriangular{T,<:StridedMatrix{T}},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(C, A, parent(L), _nthreads(), Val(false))
  return C
end
function rdiv!(
  C::StridedMatrix{T},
  A::StridedMatrix{T},
  L::LowerTriangular{T,<:StridedMatrix{T}},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(C, A, parent(L), static(1), Val(false))
  return C
end
function rdiv!(
  A::StridedMatrix{T},
  L::UnitLowerTriangular{T,<:StridedMatrix{T}},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(A, A, parent(L), _nthreads(), Val(true))
  return A
end
function rdiv!(
  A::StridedMatrix{T},
  L::UnitLowerTriangular{T,<:StridedMatrix{T}},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(A, A, parent(L), static(1), Val(true))
  return A
end
function rdiv!(
  C::StridedMatrix{T},
  A::StridedMatrix{T},
  L::UnitLowerTriangular{T,<:StridedMatrix{T}},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(C, A, parent(L), _nthreads(), Val(true))
  return C
end
function rdiv!(
  C::StridedMatrix{T},
  A::StridedMatrix{T},
  L::UnitLowerTriangular{T,<:StridedMatrix{T}},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(C, A, parent(L), static(1), Val(true))
  return C
end
function ldiv!(
  U::UpperTriangular{T,<:StridedMatrix{T}},
  A::StridedMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(
    transpose(A),
    transpose(A),
    transpose(parent(U)),
    _nthreads(),
    Val(false)
  )
  return A
end
function ldiv!(
  U::UpperTriangular{T,<:StridedMatrix{T}},
  A::StridedMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(
    transpose(A),
    transpose(A),
    transpose(parent(U)),
    static(1),
    Val(false)
  )
  return A
end
function ldiv!(
  C::StridedMatrix{T},
  U::UpperTriangular{T,<:StridedMatrix{T}},
  A::StridedMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(
    transpose(C),
    transpose(A),
    transpose(parent(U)),
    _nthreads(),
    Val(false)
  )
  return C
end
function ldiv!(
  C::StridedMatrix{T},
  U::UpperTriangular{T,<:StridedMatrix{T}},
  A::StridedMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(
    transpose(C),
    transpose(A),
    transpose(parent(U)),
    static(1),
    Val(false)
  )
  return C
end
function ldiv!(
  U::UnitUpperTriangular{T,<:StridedMatrix{T}},
  A::StridedMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(
    transpose(A),
    transpose(A),
    transpose(parent(U)),
    _nthreads(),
    Val(true)
  )
  return A
end
function ldiv!(
  U::UnitUpperTriangular{T,<:StridedMatrix{T}},
  A::StridedMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(
    transpose(A),
    transpose(A),
    transpose(parent(U)),
    static(1),
    Val(true)
  )
  return A
end
function ldiv!(
  C::StridedMatrix{T},
  U::UnitUpperTriangular{T,<:StridedMatrix{T}},
  A::StridedMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(
    transpose(C),
    transpose(A),
    transpose(parent(U)),
    _nthreads(),
    Val(true)
  )
  return C
end
function ldiv!(
  C::StridedMatrix{T},
  U::UnitUpperTriangular{T,<:StridedMatrix{T}},
  A::StridedMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch_L!(
    transpose(C),
    transpose(A),
    transpose(parent(U)),
    static(1),
    Val(true)
  )
  return C
end

# Vector right-hand sides never defer to LinearAlgebra/BLAS and never touch
# the SIMD matrix kernels: the pure-Julia sweeps below beat the matrix
# drivers' M == 1 scalar remainder and beat BLAS trsv at every measured size
# (up to 7x below n = 256, 1.1-1.9x at n = 256..2000, AVX2).

# Column-oriented sweeps for parents whose columns are the contiguous
# direction, outer-unrolled rank-4 so each pass over `x` retires four
# columns: quartering the `x` store traffic measured ~1.5x faster than the
# rank-1 sweep. The unit-diagonal `x[j]` store discipline (elided forward,
# unconditional backward) also follows measurement; flipping either
# direction was up to 1.7x slower.
@inline function _naive_vsolve_fwd!(x, A, ::Val{UNIT}) where {UNIT}
  N = length(x)
  j = 1
  @inbounds while j < N - 2
    x0 = UNIT ? x[j] : x[j] / A[j, j]
    UNIT || (x[j] = x0)
    x1 = muladd(-x0, A[j+1, j], x[j+1])
    x1 = UNIT ? x1 : x1 / A[j+1, j+1]
    x[j+1] = x1
    x2 = muladd(-x1, A[j+2, j+1], muladd(-x0, A[j+2, j], x[j+2]))
    x2 = UNIT ? x2 : x2 / A[j+2, j+2]
    x[j+2] = x2
    x3 = muladd(
      -x2,
      A[j+3, j+2],
      muladd(-x1, A[j+3, j+1], muladd(-x0, A[j+3, j], x[j+3]))
    )
    x3 = UNIT ? x3 : x3 / A[j+3, j+3]
    x[j+3] = x3
    n0 = -x0
    n1 = -x1
    n2 = -x2
    n3 = -x3
    @simd ivdep for i = (j+4):N
      x[i] = muladd(
        n0,
        A[i, j],
        muladd(
          n1,
          A[i, j+1],
          muladd(n2, A[i, j+2], muladd(n3, A[i, j+3], x[i]))
        )
      )
    end
    j += 4
  end
  @inbounds while j <= N
    xj = UNIT ? x[j] : x[j] / A[j, j]
    UNIT || (x[j] = xj)
    nxj = -xj
    @simd ivdep for i = (j+1):N
      x[i] = muladd(nxj, A[i, j], x[i])
    end
    j += 1
  end
  nothing
end
@inline function _naive_vsolve_bwd!(x, A, ::Val{UNIT}) where {UNIT}
  N = length(x)
  j = N
  @inbounds while j > 3
    x0 = UNIT ? x[j] : x[j] / A[j, j]
    x[j] = x0
    x1 = muladd(-x0, A[j-1, j], x[j-1])
    x1 = UNIT ? x1 : x1 / A[j-1, j-1]
    x[j-1] = x1
    x2 = muladd(-x1, A[j-2, j-1], muladd(-x0, A[j-2, j], x[j-2]))
    x2 = UNIT ? x2 : x2 / A[j-2, j-2]
    x[j-2] = x2
    x3 = muladd(
      -x2,
      A[j-3, j-2],
      muladd(-x1, A[j-3, j-1], muladd(-x0, A[j-3, j], x[j-3]))
    )
    x3 = UNIT ? x3 : x3 / A[j-3, j-3]
    x[j-3] = x3
    n0 = -x0
    n1 = -x1
    n2 = -x2
    n3 = -x3
    @simd ivdep for i = 1:(j-4)
      x[i] = muladd(
        n0,
        A[i, j],
        muladd(
          n1,
          A[i, j-1],
          muladd(n2, A[i, j-2], muladd(n3, A[i, j-3], x[i]))
        )
      )
    end
    j -= 4
  end
  @inbounds while j >= 1
    xj = UNIT ? x[j] : x[j] / A[j, j]
    x[j] = xj
    nxj = -xj
    @simd ivdep for i = 1:(j-1)
      x[i] = muladd(nxj, A[i, j], x[i])
    end
    j -= 1
  end
  nothing
end
# Inner-product forms for parents whose rows are the contiguous direction.
@inline function _naive_vsolve_fwd_dot!(x, A, ::Val{UNIT}) where {UNIT}
  N = length(x)
  @inbounds for i = 1:N
    s = zero(eltype(x))
    @simd for j = 1:(i-1)
      s = muladd(A[i, j], x[j], s)
    end
    xi = x[i] - s
    x[i] = UNIT ? xi : xi / A[i, i]
  end
  nothing
end
@inline function _naive_vsolve_bwd_dot!(x, A, ::Val{UNIT}) where {UNIT}
  N = length(x)
  @inbounds for i = N:-1:1
    s = zero(eltype(x))
    @simd for j = (i+1):N
      s = muladd(A[i, j], x[j], s)
    end
    xi = x[i] - s
    x[i] = UNIT ? xi : xi / A[i, i]
  end
  nothing
end
@inline function _naive_vsolve!(x, A, ::Val{UNIT}, ::Val{UP}) where {UNIT,UP}
  colmajor = abs(stride(A, 1)) <= abs(stride(A, 2))
  if UP
    colmajor ? _naive_vsolve_bwd!(x, A, Val(UNIT)) :
    _naive_vsolve_bwd_dot!(x, A, Val(UNIT))
  else
    colmajor ? _naive_vsolve_fwd!(x, A, Val(UNIT)) :
    _naive_vsolve_fwd_dot!(x, A, Val(UNIT))
  end
end

for (wrap, UNIT, UP) in (
  (:LowerTriangular, false, false),
  (:UnitLowerTriangular, true, false),
  (:UpperTriangular, false, true),
  (:UnitUpperTriangular, true, true)
)
  @eval begin
    function ldiv!(
      U::$wrap{T,<:StridedMatrix{T}},
      b::StridedVector{T},
      ::Val = Val(true)
    ) where {T<:Union{Float32,Float64}}
      P = parent(U)
      N = length(b)
      if size(P, 1) != N
        throw(
          DimensionMismatch(
            "triangular matrix is $(size(P,1))×$(size(P,2)), right-hand side has length $N"
          )
        )
      end
      _naive_vsolve!(b, P, Val($UNIT), Val($UP))
      return b
    end
    function ldiv!(
      c::StridedVector{T},
      U::$wrap{T,<:StridedMatrix{T}},
      b::StridedVector{T},
      ::Val = Val(true)
    ) where {T<:Union{Float32,Float64}}
      P = parent(U)
      N = length(b)
      if size(P, 1) != N
        throw(
          DimensionMismatch(
            "triangular matrix is $(size(P,1))×$(size(P,2)), right-hand side has length $N"
          )
        )
      elseif length(c) != N
        throw(
          DimensionMismatch("destination has length $(length(c)), needs $N")
        )
      end
      c === b || copyto!(c, b)
      _naive_vsolve!(c, P, Val($UNIT), Val($UP))
      return c
    end
  end
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
  Bsize = nothing
) where {T,UNIT}
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
  spc::AbstractStridedPointer{T,<:Any,XC},
  spa::AbstractStridedPointer{T,<:Any,XA},
  spu::AbstractStridedPointer{T,<:Any,XU},
  M,
  N,
  ::Val{UNIT}
) where {T,UNIT,XC,XA,XU}
  B = block_size(Val(T))
  W = VectorizationBase.pick_vector_width(T)
  WUF = XC == XA == XA == 2 ? W : W * unroll_factor(W)
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

struct RDivBlockMandNv2{UNIT} end
function (f::RDivBlockMandNv2{UNIT})(
  allargs,
  blockstart,
  blockstop
) where {UNIT}
  spc, spa, spu, N, Mrem, Nblock, mtb = allargs
  for block = blockstart-1:blockstop-1
    rdiv_block_MandN!(
      gesp(spc, (mtb * block, StaticInt{0}())),
      gesp(spa, (mtb * block, StaticInt{0}())),
      spu,
      Core.ifelse(block == Nblock - 1, Mrem, mtb),
      N,
      Val{UNIT}()
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
  ::Val{UNIT}
) where {UNIT,TC,TA,TU}
  # Main._a[] = (spc, spa, spu, M, N, mtb, Val(UNIT), static(X));
  (Md, Mr) = VectorizationBase.vdivrem(M, mtb)
  Nblock = Md + (Mr ≠ 0)
  Mrem = Core.ifelse(Mr ≠ 0, Mr, mtb)
  batch(
    RDivBlockMandNv2{UNIT}(),
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

@generated function _ldiv_remainder!(
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
    sub = Base.FastMath.sub_fast
    mul = Base.FastMath.mul_fast
    div = Base.FastMath.div_fast
    vlxj = :(vload(spc, ($z, j)))
    if UNIT
      vlxj = :(xj = $vlxj)
    else
      vlxj = quote
        xj = $div($vlxj, vload(spu, (j, j)))
        vstore!(spc, xj, ($z, j))
      end
    end
    quote
      $(Expr(:meta, :inline))
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
          vstore!(spc, $sub(xi, $mul(xj, Uji)), ($z, i))
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
  if W == 2
    quote
      $(Expr(:meta, :inline))
      _ldiv_remainder!(
        spc,
        spa,
        spu,
        M,
        N,
        m,
        Nr,
        $WS,
        $(Val(UNIT)),
        $(static(1))
      )
    end
  else
    quote
      # $(Expr(:meta, :inline))
      Base.Cartesian.@nif $(W - 1) w -> m == M - w w -> _ldiv_remainder!(
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
end
@inline function rdiv_U!(
  spc::AbstractStridedPointer{T,2,2},
  spa::AbstractStridedPointer{T,2,2},
  spu::AbstractStridedPointer{T,2,2},
  M,
  N,
  ::Val{UNIT}
) where {T,UNIT}
  tup = (spc, spa, spu)
  _ldiv_L!(
    M,
    N,
    Val(UNIT),
    typeof(tup),
    LoopVectorization.flatten_to_tuple(tup)...
  )
end

# spc = spa / spu
# spc' = (spu' \ spa')'
# This is ldiv
function _ldiv_L!(
  M,
  N,
  ::Val{UNIT},
  ::Type{Args},
  args::Vararg{Any,K}
) where {UNIT,Args,K}
  spc, spa, spu = LoopVectorization.reassemble_tuple(Args, args)
  T = eltype(spc)
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF * WS
  Nd, Nr = VectorizationBase.vdivrem(N, WS)
  m = 0
  # m, no remainder
  while m < M - WS + 1
    n = Nr # non factor of W remainder
    if n > 0
      BdivU_small_kern_u!(
        spc,
        nothing,
        spa,
        spu,
        n,
        StaticInt(1),
        Val(UNIT),
        WS
      )
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

# Backward-substitution family: right-lower division `C = A / L` with `L` lower
# triangular, the dual of the right-upper `rdiv!` kernels above. Column blocks
# are solved last-to-first; within and across blocks every load/store uses
# ascending (positive-stride) indices, only block offsets descend.
# `ldiv!(::UpperTriangular, A)` is obtained from this family for free via
# `U \ B = (B' / U')'`, mirroring how `ldiv!(::LowerTriangular, A)` reuses the
# right-upper kernels.

# Solve `X * L = A` for the diagonal block `L[noff .+ (0:N-1), noff .+ (0:N-1)]`:
# X_n = (A_n - \sum_{k>n} X_k L_{k,n}) / L_{n,n}, iterating n = N, N-1, ..., 1.
@generated function solve_AL(
  A::VecUnroll{Nm1},
  spl::AbstractStridedPointer,
  noff,
  ::Val{UNIT}
) where {Nm1,UNIT}
  N = Nm1 + 1
  q = quote
    $(Expr(:meta, :inline))
    Ad = VectorizationBase.data(A)
  end
  for n = 1:N
    push!(q.args, :($(Symbol(:A_, n)) = Ad[$n]))
    for m = (UNIT ? n + 1 : n):N
      push!(
        q.args,
        :(
          $(Symbol(:L_, m, :_, n)) =
            vload(spl, (noff + $(m - 1), noff + $(n - 1)))
        )
      )
    end
  end
  for n = N:-1:1
    # descending k puts the freshest operand (A_{n+1}) last, so the older
    # updates overlap its latency and the critical path stays one fma per
    # column, matching `solve_AU`'s ascending-k order
    for k = N:-1:n+1
      push!(
        q.args,
        :(
          $(Symbol(:A_, n)) = vfnmadd_fast(
            $(Symbol(:A_, k)),
            $(Symbol(:L_, k, :_, n)),
            $(Symbol(:A_, n))
          )
        )
      )
    end
    UNIT || push!(
      q.args,
      :(
        $(Symbol(:A_, n)) =
          Base.FastMath.div_fast($(Symbol(:A_, n)), $(Symbol(:L_, n, :_, n)))
      )
    )
  end
  push!(
    q.args,
    Expr(:call, :VecUnroll, Expr(:tuple, [Symbol(:A_, n) for n = 1:N]...))
  )
  q
end

# Standalone solve of the trailing `N`-column block against the lower-triangular
# diagonal block at the pointers' origin (callers `gesp` to the block).
@generated function BdivL_small_kern!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spl::AbstractStridedPointer{T},
  ::StaticInt{N},
  mask::AbstractMask{W},
  ::Val{UNIT}
) where {T,UNIT,W,N}
  z = static(0)
  if N == 1
    i = (MM{W}(z), z)
    Amn = :(vload(spb, $i, mask))
    if !UNIT
      Amn = :($Amn / vload(spl, $((z, z))))
    end
    quote
      $(Expr(:meta, :inline))
      store_small_kern!(spa, sp, $Amn, $i, mask)
    end
  else
    unroll = Unroll{2,1,N,1,W,(-1 % UInt),1}((z, z))
    tostore = :(VecUnroll(Base.Cartesian.@ntuple $N Amn))
    q = quote
      $(Expr(:meta, :inline))
      Amn = getfield(vload(spb, $unroll, mask), :data)
      Base.Cartesian.@nexprs $N n -> Amn_n = getfield(Amn, n)
    end
    for n = N:-1:1
      for k = N:-1:n+1
        push!(
          q.args,
          :(
            $(Symbol(:Amn_, n)) = vfnmadd_fast(
              $(Symbol(:Amn_, k)),
              vload(spl, ($(k - 1), $(n - 1))),
              $(Symbol(:Amn_, n))
            )
          )
        )
      end
      UNIT || push!(
        q.args,
        :($(Symbol(:Amn_, n)) /= vload(spl, ($(n - 1), $(n - 1))))
      )
    end
    push!(q.args, :(store_small_kern!(spa, sp, $tostore, $unroll, mask)))
    q
  end
end
@generated function BdivL_small_kern_u!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spl::AbstractStridedPointer{T},
  ::StaticInt{N},
  ::StaticInt{U},
  ::Val{UNIT},
  ::StaticInt{W}
) where {T,U,UNIT,N,W}
  z = static(0)
  if N == 1
    unroll = Unroll{1,W,U,1,W,zero(UInt),1}((z, z))
    Amn = :(vload(spb, $unroll))
    if !UNIT
      Amn = :($Amn / vload(spl, $((z, z))))
    end
    quote
      $(Expr(:meta, :inline))
      store_small_kern!(spa, sp, $Amn, $unroll)
    end
  else
    double_unroll =
      Unroll{2,1,N,1,W,zero(UInt),1}(Unroll{1,W,U,1,W,zero(UInt),1}((z, z)))
    tostore = :(VecUnroll(Base.Cartesian.@ntuple $N Amn))
    q = quote
      $(Expr(:meta, :inline))
      Amn = getfield(vload(spb, $double_unroll), :data)
      Base.Cartesian.@nexprs $N n -> Amn_n = getfield(Amn, n)
    end
    for n = N:-1:1
      for k = N:-1:n+1
        push!(
          q.args,
          :(
            $(Symbol(:Amn_, n)) = vfnmadd_fast(
              $(Symbol(:Amn_, k)),
              vload(spl, ($(k - 1), $(n - 1))),
              $(Symbol(:Amn_, n))
            )
          )
        )
      end
      UNIT || push!(
        q.args,
        :($(Symbol(:Amn_, n)) /= vload(spl, ($(n - 1), $(n - 1))))
      )
    end
    push!(q.args, :(store_small_kern!(spa, sp, $tostore, $double_unroll)))
    q
  end
end
@generated function BdivL_small_kern!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spl::AbstractStridedPointer{T},
  Nr::Int,
  mask::AbstractMask{W},
  ::Val{UNIT}
) where {T,UNIT,W}
  quote
    # $(Expr(:meta, :inline))
    Base.Cartesian.@nif $(W - 1) n -> n == Nr n ->
      BdivL_small_kern!(spa, sp, spb, spl, static(n), mask, $(Val(UNIT)))
  end
end
@generated function BdivL_small_kern_u!(
  spa::AbstractStridedPointer{T},
  sp,
  spb::AbstractStridedPointer{T},
  spl::AbstractStridedPointer{T},
  Nr::Int,
  ::StaticInt{U},
  ::Val{UNIT},
  ::StaticInt{W}
) where {T,U,UNIT,W}
  su = static(U)
  vu = Val(UNIT)
  sw = static(W)
  quote
    # $(Expr(:meta, :inline))
    Base.Cartesian.@nif $(W - 1) n -> n == Nr n ->
      BdivL_small_kern_u!(spa, sp, spb, spl, static(n), $su, $vu, $sw)
  end
end

@generated function rdiv_solve_W_u_L!(
  spc,
  spb,
  spa,
  spl,
  n,
  nend,
  ::StaticInt{W},
  ::StaticInt{U},
  ::Val{UNIT}
) where {W,U,UNIT}
  z = static(0)
  quote
    $(Expr(:meta, :inline))
    # C = A / L; C * L = A
    # C_{i,j} = (A_{i,j} - \sum_{k=j+1}^{N} C_{i,k}L_{k,j}) / L_{j,j}
    # Solves the block C_{i+[0,W*U), n+[0,W)}; columns [n+W, nend) were already
    # solved, so the reduction runs over them before the diagonal-block solve.
    C11 = VectorizationBase.data(
      vload(
        spa,
        $(Unroll{2,1,W,1,W,zero(UInt),1})(
          $(Unroll{1,W,U,1,W,zero(UInt),1})(($z, n))
        )
      )
    )
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    # zero-based induction (pointers pre-offset to the solved columns) so LLVM
    # shares one scaled index across the broadcast streams, as in the forward
    # kernels' `SafeCloseOpen(n)` loops
    spck = gesp(spc, ($z, n + $W))
    splk = gesp(spl, (n + $W, n))
    # nk descends so the most recently solved (just stored) columns are read
    # last, giving their stores time to drain — the mirror of the forward
    # kernels' ascending reduction, which also ends on the freshest block
    nk = nend - n - $W
    while nk > 0
      nk -= 1
      A11 = vload(spck, $(Unroll{1,W,U,1,W,zero(UInt),1})(($(StaticInt(0)), nk)))
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(splk, (nk, $z + (c - 1))), C11_c)
    end
    C11vu =
      solve_AL(VecUnroll((Base.Cartesian.@ntuple $W C11)), spl, n, $(Val(UNIT)))
    i = $(Unroll{2,1,W,1,W,zero(UInt),1})(
      $(Unroll{1,W,U,1,W,zero(UInt),1})(($z, n))
    )
    vstore!(spc, C11vu, i)
    maybestore!(spb, C11vu, i)
  end
end
@generated function rdiv_solve_W_L!(
  spc,
  spb,
  spa,
  spl,
  n,
  nend,
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
    C11 = VectorizationBase.data(
      vload(spa, $(Unroll{2,1,W,1,W,(-1 % UInt),1})(($z, n)), mask)
    )
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    spck = gesp(spc, ($z, n + $W))
    splk = gesp(spl, (n + $W, n))
    nk = nend - n - $W
    while nk > 0
      nk -= 1
      A11 = vload(spck, ($(MM{W}(z)), nk), mask)
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(splk, (nk, $z + (c - 1))), C11_c)
    end
    C11 = VecUnroll((Base.Cartesian.@ntuple $W C11))
    C11 = solve_AL(C11, spl, n, $(Val(UNIT)))
    i = $(Unroll{2,1,W,1,W,(-1 % UInt),1})(($z, n))
    $storecexpr
    maybestore!(spb, C11, i, mask)
  end
end

@inline function rdiv_L!(
  spc::AbstractStridedPointer{T},
  spa::AbstractStridedPointer,
  spl::AbstractStridedPointer,
  M,
  N,
  ::Val{UNIT}
) where {T,UNIT}
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF * WS
  Nd, Nr = VectorizationBase.vdivrem(N, WS)
  Nrstart = N - Nr # trailing remainder columns [Nrstart, N) are solved first
  spl_diag = gesp(spl, (Nrstart, Nrstart))
  m = 0
  if UF > 1
    while m < M - WU + 1
      if Nr > 0
        BdivL_small_kern_u!(
          gesp(spc, (StaticInt(0), Nrstart)),
          nothing,
          gesp(spa, (StaticInt(0), Nrstart)),
          spl_diag,
          Nr,
          UF,
          Val(UNIT),
          WS
        )
      end
      n = Nrstart
      for _ ∈ 1:Nd
        n -= W
        rdiv_solve_W_u_L!(spc, nothing, spa, spl, n, N, WS, UF, Val(UNIT))
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
    if Nr > 0
      BdivL_small_kern!(
        gesp(spc, (StaticInt(0), Nrstart)),
        nothing,
        gesp(spa, (StaticInt(0), Nrstart)),
        spl_diag,
        Nr,
        mask,
        Val(UNIT)
      )
    end
    n = Nrstart
    for _ ∈ 1:Nd
      n -= W
      rdiv_solve_W_L!(spc, nothing, spa, spl, n, N, nothing, mask, Val(UNIT))
    end
    spa = gesp(spa, (WS, StaticInt(0)))
    spc = gesp(spc, (WS, StaticInt(0)))
    m = ubm
  end
  nothing
end

function div_dispatch_L!(
  C::AbstractMatrix{T},
  A,
  L,
  nthread,
  ::Val{UNIT}
) where {UNIT,T}
  _M, _N = size(A)
  M = _canonicalize(_M)
  N = _canonicalize(_N)
  if (size(L, 1) != N) | (size(L, 2) != N)
    throw(DimensionMismatch("triangular matrix is $(size(L,1))×$(size(L,2)), needs $N×$N"))
  elseif (size(C, 1) != M) | (size(C, 2) != N)
    throw(DimensionMismatch("destination is $(size(C,1))×$(size(C,2)), needs $M×$N"))
  end
  ((N == 0) | (M == 0)) && return nothing
  _spa, spap = stridedpointer_preserve(A)
  _spc, spcp = stridedpointer_preserve(C)
  _spl, splp = stridedpointer_preserve(L)
  spa = zero_offsets(_spa)
  spc = zero_offsets(_spc)
  spl = zero_offsets(_spl)
  GC.@preserve spap spcp splp begin
    mtb = m_thread_block_size(M, N, nthread, Val(T))
    if nthread > 1
      (M > mtb) &&
        return multithread_rdiv_L!(spc, spa, spl, M, N, mtb, Val(UNIT))
    elseif N > block_size(Val(T))
      return rdiv_block_MandN_L!(spc, spa, spl, M, N, Val(UNIT))
    end
    return rdiv_L!(spc, spa, spl, M, N, Val(UNIT))
  end
end

# Writes `spd[m, n] = spa[m, n] - sum_k spc[m, k] * spl[k, n]`, where `spc`
# points at the already-solved trailing columns and `spd` at the block being
# prepared for its triangular solve (`spc`/`spd` alias disjoint column ranges).
nmuladd_L!(spc, spa, spl, spd, M, K, N) =
  @turbo for n ∈ CloseOpen(N), m ∈ CloseOpen(M)
    Cmn = spa[m, n]
    for k ∈ CloseOpen(K)
      Cmn -= spc[m, k] * spl[k, n]
    end
    spd[m, n] = Cmn
  end

function rdiv_block_N_L!(
  spc::AbstractStridedPointer{T},
  spa,
  spl,
  M,
  N,
  ::Val{UNIT},
  Bsize = nothing
) where {T,UNIT}
  W = VectorizationBase.pick_vector_width(T)
  B_normalized =
    Bsize === nothing ?
    VectorizationBase.vcld(
      N,
      VectorizationBase.vcld(N, block_size(Val(T))) * W
    ) * W : Bsize
  repeat = N > B_normalized
  N_temp = Core.ifelse(repeat, B_normalized, N)
  nstart = N - N_temp
  rdiv_L!(
    gesp(spc, (StaticInt(0), nstart)),
    gesp(spa, (StaticInt(0), nstart)),
    gesp(spl, (nstart, nstart)),
    M,
    N_temp,
    Val{UNIT}()
  )
  repeat || return nothing
  while nstart > 0
    ns = nstart
    N_temp = min(B_normalized, ns)
    nstart = ns - N_temp
    nmuladd_L!(
      gesp(spc, (StaticInt(0), ns)),
      gesp(spa, (StaticInt(0), nstart)),
      gesp(spl, (ns, nstart)),
      gesp(spc, (StaticInt(0), nstart)),
      M,
      N - ns,
      N_temp
    )
    rdiv_L!(
      gesp(spc, (StaticInt(0), nstart)),
      gesp(spc, (StaticInt(0), nstart)),
      gesp(spl, (nstart, nstart)),
      M,
      N_temp,
      Val{UNIT}()
    )
  end
  nothing
end
function rdiv_block_MandN_L!(
  spc::AbstractStridedPointer{T,<:Any,XC},
  spa::AbstractStridedPointer{T,<:Any,XA},
  spl::AbstractStridedPointer{T,<:Any,XL},
  M,
  N,
  ::Val{UNIT}
) where {T,UNIT,XC,XA,XL}
  B = block_size(Val(T))
  W = VectorizationBase.pick_vector_width(T)
  WUF = XC == XA == XL == 2 ? W : W * unroll_factor(W)
  B_m = VectorizationBase.vcld(M, VectorizationBase.vcld(M, B) * WUF) * WUF
  m = 0
  while m < M
    mu = m + B_m
    Mtemp = min(M, mu) - m
    rdiv_block_N_L!(
      spc,
      spa,
      spl,
      Mtemp,
      N,
      Val{UNIT}(),
      VectorizationBase.vcld(N, VectorizationBase.vcld(N, B) * W) * W
    )
    spa = gesp(spa, (B_m, StaticInt{0}()))
    spc = gesp(spc, (B_m, StaticInt{0}()))
    m = mu
  end
  nothing
end

struct RDivBlockMandNv2L{UNIT} end
function (f::RDivBlockMandNv2L{UNIT})(
  allargs,
  blockstart,
  blockstop
) where {UNIT}
  spc, spa, spl, N, Mrem, Nblock, mtb = allargs
  for block = blockstart-1:blockstop-1
    rdiv_block_MandN_L!(
      gesp(spc, (mtb * block, StaticInt{0}())),
      gesp(spa, (mtb * block, StaticInt{0}())),
      spl,
      Core.ifelse(block == Nblock - 1, Mrem, mtb),
      N,
      Val{UNIT}()
    )
  end
end

function multithread_rdiv_L!(
  spc::AbstractStridedPointer{TC},
  spa::AbstractStridedPointer{TA},
  spl::AbstractStridedPointer{TL},
  M::Int,
  N::Int,
  mtb::Int,
  ::Val{UNIT}
) where {UNIT,TC,TA,TL}
  (Md, Mr) = VectorizationBase.vdivrem(M, mtb)
  Nblock = Md + (Mr ≠ 0)
  Mrem = Core.ifelse(Mr ≠ 0, Mr, mtb)
  batch(
    RDivBlockMandNv2L{UNIT}(),
    (Nblock, min(Nblock, Threads.nthreads())),
    spc,
    spa,
    spl,
    N,
    Mrem,
    Nblock,
    mtb
  )
  nothing
end

@generated function uldiv_solve_W_u!(
  spc,
  spa,
  spl,
  n,
  nend,
  ::StaticInt{W},
  ::StaticInt{U},
  ::Val{UNIT}
) where {W,U,UNIT}
  z = static(0)
  q = quote
    $(Expr(:meta, :inline))
    # C = A / L in the transposed (row-major) layout; equivalently, C' = L' \ A'
    # with L' upper triangular — backward substitution for the left-upper solve.
    # Solves the block C_{j+[0,W), n+[0,W*U)}; columns [n+W*U, nend) are already
    # solved, so the reduction runs over them, then the W*U block columns are
    # solved right-to-left in W-sized sub-blocks against L's diagonal block.
    # outer unroll are `W` rows
    # Inner unroll are `W*U` columns (U simd vecs)
    A11 = getfield(
      vload(
        spa,
        $(Unroll{1,1,W,2,W,zero(UInt),1})(
          $(Unroll{2,W,U,2,W,zero(UInt),1})(($z, n))
        )
      ),
      :data
    )
    Base.Cartesian.@nexprs $W c -> A11_c = getfield(A11, c)
    # zero-based induction (pointers pre-offset to the solved columns) so LLVM
    # shares one scaled index across the broadcast streams, as in the forward
    # kernels' `SafeCloseOpen(n)` loops
    spck = gesp(spc, ($z, n + $(W * U)))
    splk = gesp(spl, (n + $(W * U), n))
    nk = nend - n - $(W * U)
    while nk > 0
      nk -= 1
      L_ki = vload(splk, $(Unroll{2,W,U,2,W,zero(UInt),1})((nk, $z)))
      Base.Cartesian.@nexprs $W c ->
        A11_c = vfnmadd_fast(L_ki, vload(spck, (static(c - 1), nk)), A11_c)
    end
  end
  for u = U:-1:1
    Xtup = Expr(:tuple)
    for w = 1:W
      push!(Xtup.args, :(getfield(getfield($(Symbol(:A11_, w)), :data), $u)))
    end
    push!(
      q.args,
      :(
        $(Symbol(:X_, u)) = getfield(
          VectorizationBase.transpose_vecunroll(VecUnroll($Xtup)),
          :data
        )
      )
    )
    for c = 1:W
      push!(q.args, :($(Symbol(:X_, u, :_, c)) = getfield($(Symbol(:X_, u)), $c)))
    end
    # subtract the contributions of the already-solved sub-blocks to the right
    for j = U:-1:u+1
      for k = 1:W
        for c = 1:W
          push!(
            q.args,
            :(
              $(Symbol(:X_, u, :_, c)) = vfnmadd_fast(
                $(Symbol(:C_, j, :_, k)),
                vload(
                  spl,
                  (
                    n + $((k - 1) + (j - 1) * W),
                    n + $((c - 1) + (u - 1) * W)
                  )
                ),
                $(Symbol(:X_, u, :_, c))
              )
            )
          )
        end
      end
    end
    Xtup2 = Expr(:tuple)
    for c = 1:W
      push!(Xtup2.args, Symbol(:X_, u, :_, c))
    end
    push!(
      q.args,
      :(
        $(Symbol(:C_, u)) = solve_AL(
          VecUnroll($Xtup2),
          spl,
          n + $((u - 1) * W),
          $(Val(UNIT))
        )
      )
    )
    push!(q.args, :($(Symbol(:Cdata_, u)) = getfield($(Symbol(:C_, u)), :data)))
    for c = 1:W
      push!(
        q.args,
        :($(Symbol(:C_, u, :_, c)) = getfield($(Symbol(:Cdata_, u)), $c))
      )
    end
  end
  # store at end (no aliasing)
  for u = 1:U
    push!(
      q.args,
      :(vstore!(
        spc,
        $(Symbol(:C_, u)),
        $(Unroll{2,1,W,1,W,zero(UInt),1})(($z, n + $((u - 1) * W)))
      ))
    )
  end
  q
end
@generated function uldiv_solve_W!(
  spc,
  spa,
  spl,
  n,
  nend,
  ::StaticInt{W},
  ::Val{UNIT}
) where {W,UNIT}
  z = static(0)
  quote
    $(Expr(:meta, :inline))
    # Like `uldiv_solve_W_u!`, except no unrolling, just a `W`x`W` block
    A11 =
      getfield(vload(spa, $(Unroll{1,1,W,2,W,zero(UInt),1})(($z, n))), :data)
    Base.Cartesian.@nexprs $W c -> A11_c = getfield(A11, c)
    spck = gesp(spc, ($z, n + $W))
    splk = gesp(spl, (n + $W, n))
    nk = nend - n - $W
    while nk > 0
      nk -= 1
      L_ki = vload(splk, (nk, $(MM{W}(z))))
      Base.Cartesian.@nexprs $W c ->
        A11_c = vfnmadd_fast(L_ki, vload(spck, (static(c - 1), nk)), A11_c)
    end
    X = VectorizationBase.transpose_vecunroll(
      VecUnroll(Base.Cartesian.@ntuple $W A11)
    )
    C_u = solve_AL(X, spl, n, $(Val(UNIT)))
    vstore!(spc, C_u, $(Unroll{2,1,W,1,W,zero(UInt),1})(($z, n)))
  end
end
@generated function uldiv_solve_W!(
  spc,
  spa,
  spl,
  n,
  nend,
  ::StaticInt{W},
  ::Val{UNIT},
  ::StaticInt{R}
) where {W,UNIT,R}
  R <= 1 && throw("Remainder of `<= 1` shouldn't be called, but had $R.")
  R >= W && throw("Remainder of `>= $W` shouldn't be called, but had $R.")
  z = static(0)
  q = quote
    # Like the `R`-remainder `ldiv_solve_W!`, but backward: solves the block
    # C_{j+[0,R), n+[0,W)} with columns [n+W, nend) already solved.
    A11 =
      getfield(vload(spa, $(Unroll{1,1,R,2,W,zero(UInt),1})(($z, n))), :data)
    Base.Cartesian.@nexprs $R r -> A11_r = getfield(A11, r)
    spck = gesp(spc, ($z, n + $W))
    splk = gesp(spl, (n + $W, n))
    nk = nend - n - $W
    while nk > 0
      nk -= 1
      L_ki = vload(splk, (nk, $(MM{W}(z))))
      Base.Cartesian.@nexprs $R r ->
        A11_r = vfnmadd_fast(L_ki, vload(spck, (static(r - 1), nk)), A11_r)
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
    X = VectorizationBase.transpose_vecunroll(VecUnroll($t))
    C_u = solve_AL(X, spl, n, $(Val(UNIT)))
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

@generated function _uldiv_remainder!(
  spc,
  spa,
  spl,
  M,
  N,
  m,
  Nr,
  ::StaticInt{W},
  ::Val{UNIT},
  ::StaticInt{r}
) where {W,UNIT,r}
  r <= 0 && throw("Remainder of `<= 0` shouldn't be called, but had $r.")
  r >= W && throw("Remainder of `>= $W` shouldn't be called, but had $r.")
  if r == 1
    z = static(0)
    sub = Base.FastMath.sub_fast
    mul = Base.FastMath.mul_fast
    div = Base.FastMath.div_fast
    vlxj = :(vload(spc, ($z, j)))
    if UNIT
      vlxj = :(xj = $vlxj)
    else
      vlxj = quote
        xj = $div($vlxj, vload(spl, (j, j)))
        vstore!(spc, xj, ($z, j))
      end
    end
    quote
      $(Expr(:meta, :inline))
      if pointer(spc) != pointer(spa)
        for n = 0:N-1
          vstore!(spc, vload(spa, ($z, n)), ($z, n))
        end
      end
      for j = N-1:-1:0
        $vlxj
        for i = 0:j-1
          xi = vload(spc, ($z, i))
          Lji = vload(spl, (j, i))
          vstore!(spc, $sub(xi, $mul(xj, Lji)), ($z, i))
        end
      end
    end
  else
    WS = static(W)
    quote
      $(Expr(:meta, :inline))
      n = N - Nr
      if Nr > 0
        mask = $(VectorizationBase.mask(WS, r))
        BdivL_small_kern!(
          gesp(spc, (StaticInt(0), n)),
          nothing,
          gesp(spa, (StaticInt(0), n)),
          gesp(spl, (n, n)),
          Nr,
          mask,
          $(Val(UNIT))
        )
      end
      while n != 0
        n -= $W
        uldiv_solve_W!(spc, spa, spl, n, N, $WS, $(Val(UNIT)), $(StaticInt(r)))
      end
    end
  end
end
@generated function uldiv_remainder!(
  spc,
  spa,
  spl,
  M,
  N,
  m,
  Nr,
  ::StaticInt{W},
  ::Val{UNIT}
) where {W,UNIT}
  WS = static(W)
  if W == 2
    quote
      $(Expr(:meta, :inline))
      _uldiv_remainder!(
        spc,
        spa,
        spl,
        M,
        N,
        m,
        Nr,
        $WS,
        $(Val(UNIT)),
        $(static(1))
      )
    end
  else
    quote
      # $(Expr(:meta, :inline))
      Base.Cartesian.@nif $(W - 1) w -> m == M - w w -> _uldiv_remainder!(
        spc,
        spa,
        spl,
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
end
@inline function rdiv_L!(
  spc::AbstractStridedPointer{T,2,2},
  spa::AbstractStridedPointer{T,2,2},
  spl::AbstractStridedPointer{T,2,2},
  M,
  N,
  ::Val{UNIT}
) where {T,UNIT}
  tup = (spc, spa, spl)
  _ldiv_U!(
    M,
    N,
    Val(UNIT),
    typeof(tup),
    LoopVectorization.flatten_to_tuple(tup)...
  )
end

# spc = spa / spl with the pointers transposed (row-major)
# spc' = (spl' \ spa')'
# This is the left-upper ldiv
function _ldiv_U!(
  M,
  N,
  ::Val{UNIT},
  ::Type{Args},
  args::Vararg{Any,K}
) where {UNIT,Args,K}
  spc, spa, spl = LoopVectorization.reassemble_tuple(Args, args)
  T = eltype(spc)
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF * WS
  Nd, Nr = VectorizationBase.vdivrem(N, WS)
  Nrstart = N - Nr # trailing remainder columns [Nrstart, N) are solved first
  m = 0
  # m, no remainder
  while m < M - WS + 1
    if Nr > 0
      BdivL_small_kern_u!(
        gesp(spc, (StaticInt(0), Nrstart)),
        nothing,
        gesp(spa, (StaticInt(0), Nrstart)),
        gesp(spl, (Nrstart, Nrstart)),
        Nr,
        StaticInt(1),
        Val(UNIT),
        WS
      )
    end
    n = Nrstart
    while n >= WU
      n -= WU
      uldiv_solve_W_u!(spc, spa, spl, n, N, WS, UF, Val(UNIT))
    end
    while n != 0
      n -= W
      uldiv_solve_W!(spc, spa, spl, n, N, WS, Val(UNIT))
    end
    m += W
    spa = gesp(spa, (W, StaticInt(0)))
    spc = gesp(spc, (W, StaticInt(0)))
  end
  # remainder on `m`
  m < M && uldiv_remainder!(spc, spa, spl, M, N, m, Nr, WS, Val(UNIT))
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
