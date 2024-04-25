module TriangularSolve
using Base: @nexprs, @ntuple
if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
  @eval Base.Experimental.@max_methods 1
end

using LayoutPointers: stridedpointer_preserve
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
using StaticArrayInterface

const LPtr{T} = Core.LLVMPtr{T,0}
_lptr(x::Ptr{T}) where {T} = reinterpret(LPtr{T}, x)
_lptr(x) = x
_ptr(x::LPtr{T}) where {T} = reinterpret(Ptr{T}, x)
_ptr(x) = x
@inline reassemble_tup(::Type{T}, t) where {T} =
  LoopVectorization.reassemble_tuple(T, map(_ptr, t))
@inline flatten_to_tup(t) = map(_lptr, LoopVectorization.flatten_to_tuple(t))

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

@generated function BdivU_small_kern!(
  ::StaticInt{N},
  _mask::UInt32,
  ::StaticInt{W},
  ::Val{UNIT},
  ::Type{Args},
  args::Vararg{Any,K}
) where {UNIT,W,N,Args,K}
  z = static(0)
  if N == 1
    i = (MM{W}(z), z)
    Amn = :(vload(spa, $i, mask))
    if !UNIT
      Amn = :($Amn / vload(spu, $((z, z))))
    end
    quote
      $(Expr(:meta, :inline))
      mask = $(VectorizationBase.Mask{W})(_mask)
      spa, spu = reassemble_tup($Args, args)
      vstore!(spa, $Amn, $i, mask)
    end
  else
    unroll = Unroll{2,1,N,1,W,(-1 % UInt),1}((z, z))
    tostore = :(VecUnroll(Base.Cartesian.@ntuple $N Amn))
    scale = UNIT ? nothing : :(Amn_n /= vload(spu, (n - 1, n - 1)))
    quote
      $(Expr(:meta, :inline))
      spa, spu = reassemble_tup($Args, args)
      mask = $(VectorizationBase.Mask{W})(_mask)
      Amn = getfield(vload(spa, $unroll, mask), :data)
      Base.Cartesian.@nexprs $N n -> begin
        Amn_n = getfield(Amn, n)
        Base.Cartesian.@nexprs (n - 1) k -> begin
          Amn_n = vfnmadd_fast(Amn_k, vload(spu, (k - 1, n - 1)), Amn_n)
        end
        $scale
      end
      vstore!(spa, $tostore, $unroll, mask)
    end
  end
end
@generated function BdivU_small_kern!(
  Nr::Int,
  mask::UInt32,
  ::StaticInt{W},
  ::Val{UNIT},
  ::Type{Args},
  args::Vararg{Any,K}
) where {UNIT,W,Args,K}
  WS = static(W)
  quote
    # $(Expr(:meta, :inline))
    Base.Cartesian.@nif $(W - 1) n -> n == Nr n ->
      BdivU_small_kern!(static(n), mask, $WS, $(Val(UNIT)), $Args, args...)
  end
end
@generated function BdivU_small_kern_u!(
  ::StaticInt{N},
  ::StaticInt{U},
  ::Val{UNIT},
  ::StaticInt{W},
  ::Type{Args},
  args::Vararg{Any,K}
) where {U,UNIT,N,W,Args,K}
  z = static(0)
  if N == 1
    unroll = Unroll{1,W,U,1,W,zero(UInt),1}((z, z))
    Amn = :(vload(spa, $unroll))
    if !UNIT
      Amn = :($Amn / vload(spu, $((z, z))))
    end
    quote
      $(Expr(:meta, :inline))
      spa, spu = reassemble_tup($Args, args)
      vstore!(spa, $Amn, $unroll)
    end
  else
    double_unroll =
      Unroll{2,1,N,1,W,zero(UInt),1}(Unroll{1,W,U,1,W,zero(UInt),1}((z, z)))
    tostore = :(VecUnroll(Base.Cartesian.@ntuple $N Amn))
    scale = UNIT ? nothing : :(Amn_n /= vload(spu, (n - 1, n - 1)))
    quote
      $(Expr(:meta, :inline))
      spa, spu = reassemble_tup($Args, args)
      Amn = getfield(vload(spa, $double_unroll), :data)
      Base.Cartesian.@nexprs $N n -> begin
        Amn_n = getfield(Amn, n)
        Base.Cartesian.@nexprs (n - 1) k -> begin
          Amn_n = vfnmadd_fast(Amn_k, vload(spu, (k - 1, n - 1)), Amn_n)
        end
        $scale
      end
      vstore!(spa, $tostore, $double_unroll)
    end
  end
end

@generated function BdivU_small_kern_u!(
  Nr::Int,
  ::StaticInt{U},
  ::Val{UNIT},
  ::StaticInt{W},
  ::Type{Args},
  args::Vararg{Any,K}
) where {U,UNIT,W,Args,K}
  su = static(U)
  vu = Val(UNIT)
  sw = static(W)
  quote
    # $(Expr(:meta, :inline))
    Base.Cartesian.@nif $(W - 1) n -> n == Nr n ->
      BdivU_small_kern_u!(static(n), $su, $vu, $sw, $Args, args...)
  end
end

@generated function rdiv_solve_W_u!(
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
      A11 = vload(spa, $(Unroll{1,W,U,1,W,zero(UInt),1})(($(StaticInt(0)), nk)))
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(spu, (nk, n + (c - 1))), C11_c)
    end
    C11vu =
      solve_AU(VecUnroll((Base.Cartesian.@ntuple $W C11)), spu, n, $(Val(UNIT)))
    i = $(Unroll{2,1,W,1,W,zero(UInt),1})(
      $(Unroll{1,W,U,1,W,zero(UInt),1})(($z, n))
    )
    vstore!(spa, C11vu, i)
  end
end
@generated function rdiv_solve_W!(
  spa,
  spu,
  n,
  mask::AbstractMask{W},
  ::Val{UNIT}
) where {W,UNIT}
  z = static(0)
  quote
    $(Expr(:meta, :inline))
    # here, we just want to load the vectors
    C11 = VectorizationBase.data(
      vload(spa, $(Unroll{2,1,W,1,W,(-1 % UInt),1})(($z, n)), mask)
    )
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    for nk ∈ SafeCloseOpen(n) # nmuladd
      A11 = vload(spa, ($(MM{W}(z)), nk), mask)
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(spu, (nk, n + (c - 1))), C11_c)
    end
    C11 = VecUnroll((Base.Cartesian.@ntuple $W C11))
    C11 = solve_AU(C11, spu, n, $(Val(UNIT)))
    i = $(Unroll{2,1,W,1,W,(-1 % UInt),1})(($z, n))
    vstore!(spa, C11, i, mask)
  end
end

@generated function ldiv_solve_W_u!(
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
        A11_c = vfnmadd_fast(U_ki, vload(spa, (static(c - 1), nk)), A11_c)
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
      vstore!(spa, C_u, $(Unroll{2,1,W,1,W,zero(UInt),1})(($z, n + (u - 1) * $W)))
    end
  end
end
@generated function ldiv_solve_W!(
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
        A11_c = vfnmadd_fast(U_ki, vload(spa, (static(c - 1), nk)), A11_c)
    end
    # solve AU wants us to transpose
    # We then have column-major multiplies
    # take A[(u-1)*W,u*W), [0,W)]
    X = VectorizationBase.transpose_vecunroll(
      VecUnroll(Base.Cartesian.@ntuple $W A11)
    )
    C_u = solve_AU(X, spu, n, $(Val(UNIT)))
    vstore!(spa, C_u, $(Unroll{2,1,W,1,W,zero(UInt),1})(($z, n)))
  end
end
@inline _mask(x, y) = VectorizationBase.Mask(VectorizationBase.mask(x, y))
@generated function ldiv_solve_W!(
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
        A11_r = vfnmadd_fast(U_ki, vload(spa, (static(r - 1), nk)), A11_r)
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
      vstore!(spa, C_u, i)
    end
  else
    quote
      mask = _mask($(static(Wpad)), $(static(R)))
      i = $(Unroll{2,1,W,1,Wpad,(-1 % UInt),1})(($z, n))
      vstore!(spa, C_u, i, mask)
    end
  end
  push!(q.args, q3)
  return q
end

@inline function rdiv_U!(
  spa::AbstractStridedPointer{T},
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
        let t = (spa, spu), ft = flatten_to_tup(t)
          BdivU_small_kern_u!(n, UF, Val(UNIT), WS, typeof(t), ft...)
        end
      end
      for _ ∈ 1:Nd
        rdiv_solve_W_u!(spa, spu, n, WS, UF, Val(UNIT))
        n += W
      end
      m += WU
      spa = gesp(spa, (WU, StaticInt(0)))
    end
  end
  finalmask = _mask(WS, M)
  while m < M
    ubm = m + W
    nomaskiter = ubm < M
    mask = nomaskiter ? VectorizationBase.max_mask(WS) : finalmask
    n = Nr
    if n > 0
      let t = (spa, spu),
        ft = flatten_to_tup(t),
        mask = getfield(mask, :u) % UInt32

        BdivU_small_kern!(n, mask, WS, Val(UNIT), typeof(t), ft...)
      end
    end
    for _ ∈ 1:Nd
      rdiv_solve_W!(spa, spu, n, mask, Val(UNIT))
      n += W
    end
    spa = gesp(spa, (WS, StaticInt(0)))
    m = ubm
  end
  nothing
end

_canonicalize(x) = signed(x)
_canonicalize(::StaticInt{N}) where {N} = StaticInt{N}()
function div_dispatch!(
  A::AbstractMatrix{T},
  U,
  nthread,
  ::Val{UNIT}
) where {UNIT,T}
  _M, _N = size(A)
  M = _canonicalize(_M)
  N = _canonicalize(_N)
  ((N == 0) | (M == 0)) && return nothing
  _spa, spap = stridedpointer_preserve(A)
  _spu, spup = stridedpointer_preserve(U)
  spa = zero_offsets(_spa)
  spu = zero_offsets(_spu)
  GC.@preserve spap spup begin
    mtb = m_thread_block_size(M, N, nthread, Val(T))
    if nthread > 1
      (M > mtb) && return multithread_rdiv!(spa, spu, M, N, mtb, Val(UNIT))
    elseif N > block_size(Val(T))
      return rdiv_block_MandN!(spa, spu, M, N, Val(UNIT))
    end
    return rdiv_U!(spa, spu, M, N, Val(UNIT))
  end
end

_nthreads() =
  min(Int(VectorizationBase.num_cores())::Int, Threads.nthreads()::Int)
function rdiv!(
  A::AbstractMatrix{T},
  U::UpperTriangular{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, parent(U), _nthreads(), Val(false))
  return A
end
function rdiv!(
  A::AbstractMatrix{T},
  U::UpperTriangular{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, parent(U), static(1), Val(false))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UpperTriangular{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(copyto!(C, A), parent(U), _nthreads(), Val(false))
  return C
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UpperTriangular{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(copyto!(C, A), parent(U), static(1), Val(false))
  return C
end
function rdiv!(
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, parent(U), _nthreads(), Val(true))
  return A
end
function rdiv!(
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(A, parent(U), static(1), Val(true))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(copyto!(C, A), parent(U), _nthreads(), Val(true))
  return C
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(copyto!(C, A), parent(U), static(1), Val(true))
  return C
end
function ldiv!(
  U::LowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(A), transpose(parent(U)), _nthreads(), Val(false))
  return A
end
function ldiv!(
  U::LowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(A), transpose(parent(U)), static(1), Val(false))
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::LowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(copyto!(C, A)),
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
    transpose(copyto!(C, A)),
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
  div_dispatch!(transpose(A), transpose(parent(U)), _nthreads(), Val(true))
  return A
end
function ldiv!(
  U::UnitLowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  div_dispatch!(transpose(A), transpose(parent(U)), static(1), Val(true))
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::UnitLowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{true} = Val(true)
) where {T<:Union{Float32,Float64}}
  div_dispatch!(
    transpose(copyto!(C, A)),
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
    transpose(copyto!(C, A)),
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

struct Mat{T,ColMajor} <: AbstractMatrix{T}
  p::Ptr{T}
  x::Int
  M::Int
  N::Int
end
Base.size(A::Mat) = (A.M, A.N)
Base.axes(A::Mat) = (CloseOpen(A.M), CloseOpen(A.N))
Base.strides(A::Mat{T,true}) where {T} = (1, getfield(A, :x))
Base.strides(A::Mat{T,false}) where {T} = (getfield(A, :x), 1)
Base.transpose(A::Mat{T,true}) where {T} = Mat{T,false}(A.p, A.x, A.N, A.M)
Base.transpose(A::Mat{T,false}) where {T} = Mat{T,true}(A.p, A.x, A.N, A.M)
Base.pointer(A::Mat) = getfield(A, :p)
StaticArrayInterface.device(::Mat) = StaticArrayInterface.CPUPointer()
StaticArrayInterface.static_strides(A::Mat{T,true}) where {T} =
  (static(1), getfield(A, :x))
StaticArrayInterface.static_strides(A::Mat{T,false}) where {T} =
  (getfield(A, :x), static(1))
StaticArrayInterface.offsets(::Mat) = (static(0), static(0))
StaticArrayInterface.stride_rank(::Type{<:Mat{<:Any,true}}) = (static(1), static(2))
StaticArrayInterface.stride_rank(::Type{<:Mat{<:Any,false}}) = (static(2), static(1))
StaticArrayInterface.contiguous_batch_size(::Type{<:Mat}) = static(0)
StaticArrayInterface.dense_dims(::Type{<:Mat{<:Any,true}}) = (static(true),static(false))
StaticArrayInterface.dense_dims(::Type{<:Mat{<:Any,false}}) = (static(false),static(true))
StaticArrayInterface.contiguous_axis(::Type{<:Mat{<:Any,true}}) = static(1)
StaticArrayInterface.contiguous_axis(::Type{<:Mat{<:Any,false}}) = static(2)
@inline function Base.getindex(
  A::Mat{T,ColMajor},
  i::Int,
  j::Int
) where {T,ColMajor}
  (; p, x) = A
  offset = ColMajor ? i + j * x : i * x + j
  unsafe_load(p, offset + 1)
end
@inline function Base.setindex!(
  A::Mat{T,ColMajor},
  v::T,
  i::Int,
  j::Int
) where {T,ColMajor}
  (; p, x) = A
  offset = ColMajor ? i + j * x : i * x + j
  unsafe_store!(p, v, offset + 1)
  v
end
@inline function Mat(A::AbstractMatrix{T}) where {T}
  r, c = LoopVectorization.ArrayInterface.stride_rank(A)
  M, N = size(A)
  if r === static(1)
    Mat{T,true}(pointer(A), stride(A, 2), M, N)
  else
    @assert c === static(1)
    Mat{T,false}(pointer(A), stride(A, 1), M, N)
  end
end

# C -= A * B
@inline function _schur_complement!(C::Mat, A::Mat, B::Mat, ::Val{false})
  # _turbo_! will not be inlined
  @turbo warn_check_args = false for n in indices((C, B), 2),
    m in indices((C, A), 1)

    Cmn = C[m, n]
    for k in indices((A, B), (2, 1))
      Cmn -= A[m, k] * B[k, n]
    end
    C[m, n] = Cmn
  end
end
@inline function _schur_complement!(C::Mat, A::Mat, B::Mat, ::Val{true})
  # _turbo_! will not be inlined
  @tturbo warn_check_args = false for n in indices((C, B), 2),
    m in indices((C, A), 1)

    Cmn = C[m, n]
    for k in indices((A, B), (2, 1))
      Cmn -= A[m, k] * B[k, n]
    end
    C[m, n] = Cmn
  end
end
@inline function schur_complement!(
  C::Mat,
  A::Mat{<:Any,false},
  B::Mat{<:Any,false},
  ::Val{THREAD}
) where {THREAD}
  # C - A * B == (C' - B' * A')'
  _schur_complement!(transpose(C), transpose(B), transpose(A), Val(THREAD))
end
@inline function schur_complement!(
  C::Mat,
  A::Mat,
  B::Mat,
  ::Val{THREAD}
) where {THREAD}
  _schur_complement!(C, A, B, Val(THREAD))
end
@inline function schur_complement!(C, A, B, ::Val{THREAD}) where {THREAD}
  schur_complement!(Mat(C), Mat(A), Mat(B), Val(THREAD))
end

@inline function Mat(sp::StridedPointer{T,2,1}, M, N) where {T}
  x, y = strides(stridedpointer(sp))
  st = sizeof(T)
  @assert x == st
  Mat{T,true}(pointer(sp), y >>> trailing_zeros(st), M, N)
end
@inline function Mat(sp::StridedPointer{T,2,2}, M, N) where {T}
  x, y = strides(stridedpointer(sp))
  st = sizeof(T)
  @assert y == st
  Mat{T,false}(pointer(sp), x >>> trailing_zeros(st), M, N)
end

function rdiv_block_N!(
  spa::AbstractStridedPointer{T},
  spu,
  M,
  N,
  ::Val{UNIT},
  Bsize = nothing
) where {T,UNIT}
  spa_base = spa
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
    rdiv_U!(spa, gesp(spu, (n, StaticInt{0}())), M, N_temp, Val{UNIT}())
    repeat || break
    spa = gesp(spa, (StaticInt(0), B_normalized))
    spu = gesp(spu, (StaticInt(0), B_normalized))
    n += B_normalized
    repeat = n + B_normalized < N
    N_temp = repeat ? N_temp : N - n
    schur_complement!(
      Mat(spa, M, N_temp),
      Mat(spa_base, M, n),
      Mat(spu, n, N_temp),
      Val(false)
    )
  end
end
function rdiv_block_MandN!(
  spa::AbstractStridedPointer{T,<:Any,XA},
  spu::AbstractStridedPointer{T,<:Any,XU},
  M,
  N,
  ::Val{UNIT}
) where {T,UNIT,XA,XU}
  B = block_size(Val(T))
  W = VectorizationBase.pick_vector_width(T)
  WUF = XA == XA == 2 ? W : W * unroll_factor(W)
  B_m = VectorizationBase.vcld(M, VectorizationBase.vcld(M, B) * WUF) * WUF
  m = 0
  while m < M
    mu = m + B_m
    Mtemp = min(M, mu) - m
    rdiv_block_N!(
      spa,
      spu,
      Mtemp,
      N,
      Val{UNIT}(),
      VectorizationBase.vcld(N, VectorizationBase.vcld(N, B) * W) * W
    )
    spa = gesp(spa, (B_m, StaticInt{0}()))
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
  spa, spu, N, Mrem, Nblock, mtb = allargs
  for block = blockstart-1:blockstop-1
    rdiv_block_MandN!(
      gesp(spa, (mtb * block, StaticInt{0}())),
      spu,
      Core.ifelse(block == Nblock - 1, Mrem, mtb),
      N,
      Val{UNIT}()
    )
  end
end

function multithread_rdiv!(
  spa::AbstractStridedPointer{TA},
  spu::AbstractStridedPointer{TU},
  M::Int,
  N::Int,
  mtb::Int,
  ::Val{UNIT}
) where {UNIT,TA,TU}
  (Md, Mr) = VectorizationBase.vdivrem(M, mtb)
  Nblock = Md + (Mr ≠ 0)
  Mrem = Core.ifelse(Mr ≠ 0, Mr, mtb)
  batch(
    RDivBlockMandNv2{UNIT}(),
    (Nblock, min(Nblock, Threads.nthreads())),
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
  spa,
  spu,
  N,
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
    vlxj = :(vload(spa, ($z, j)))
    if UNIT
      vlxj = :(xj = $vlxj)
    else
      vlxj = quote
        xj = $div($vlxj, vload(spu, (j, j)))
        vstore!(spa, xj, ($z, j))
      end
    end
    quote
      $(Expr(:meta, :inline))
      for j = 0:N-1
        $vlxj
        for i = (j+1):N-1
          xi = vload(spa, ($z, i))
          Uji = vload(spu, (j, i))
          vstore!(spa, $sub(xi, $mul(xj, Uji)), ($z, i))
        end
      end
    end
  else
    WS = static(W)
    quote
      $(Expr(:meta, :inline))
      n = Nr # non factor of W remainder
      if n > 0
        let t = (spa, spu),
          ft = flatten_to_tup(t),
          mask = $(getfield(_mask(WS, r), :u) % UInt32)

          BdivU_small_kern!(n, mask, $WS, $(Val(UNIT)), typeof(t), ft...)
        end
      end
      # while n < N - $(W * U - 1)
      #   ldiv_solve_W_u!(spa, spa, spu, n, $WS, $US, Val(UNIT), Val(r))
      #   n += $(W * U)
      # end
      while n != N
        ldiv_solve_W!(spa, spu, n, $WS, $(Val(UNIT)), $(StaticInt(r)))
        n += $W
      end
    end
  end
end
@generated function ldiv_remainder!(
  M,
  N,
  m,
  Nr,
  ::StaticInt{W},
  # ::Val{U},
  ::Val{UNIT},
  ::Type{Args},
  args...
) where {W,UNIT,Args}
  WS = static(W)
  # US = static(U)
  if W == 2
    quote
      $(Expr(:meta, :inline))
      spa, spu = reassemble_tup(Args, args)
      _ldiv_remainder!(spa, spu, N, Nr, $WS, $(Val(UNIT)), $(static(1)))
      nothing
    end
  elseif W == 8
    quote
      # $(Expr(:meta, :inline))
      spa, spu = reassemble_tup(Args, args)
      if m == M - 1
        _ldiv_remainder!(spa, spu, N, Nr, static(8), $(Val(UNIT)), StaticInt(1))
      else
        if m == M - 2
          _ldiv_remainder!(
            spa,
            spu,
            N,
            Nr,
            static(8),
            $(Val(UNIT)),
            StaticInt(2)
          )
        else
          if m == M - 3
            _ldiv_remainder!(
              spa,
              spu,
              N,
              Nr,
              static(8),
              $(Val(UNIT)),
              StaticInt(3)
            )
          else
            if m == M - 4
              _ldiv_remainder!(
                spa,
                spu,
                N,
                Nr,
                static(8),
                $(Val(UNIT)),
                StaticInt(4)
              )
            else
              if m == M - 5
                _ldiv_remainder!(
                  spa,
                  spu,
                  N,
                  Nr,
                  static(8),
                  $(Val(UNIT)),
                  StaticInt(5)
                )
              else
                if m == M - 6
                  _ldiv_remainder!(
                    spa,
                    spu,
                    N,
                    Nr,
                    static(8),
                    $(Val(UNIT)),
                    StaticInt(6)
                  )
                else
                  _ldiv_remainder!(
                    spa,
                    spu,
                    N,
                    Nr,
                    static(8),
                    $(Val(UNIT)),
                    StaticInt(7)
                  )
                end
              end
            end
          end
        end
      end
      nothing
    end
  else
    quote
      # $(Expr(:meta, :inline))
      spa, spu = reassemble_tup(Args, args)
      Base.Cartesian.@nif $(W - 1) w -> m == M - w w ->
        _ldiv_remainder!(spa, spu, N, Nr, $WS, $(Val(UNIT)), static(w))
      nothing
    end
  end
end
@inline function rdiv_U!(
  spa::AbstractStridedPointer{T,2,2},
  spu::AbstractStridedPointer{T,2,2},
  M,
  N,
  ::Val{UNIT}
) where {T,UNIT}
  tup = (spa, spu)
  ftup = flatten_to_tup(tup)
  _ldiv_L!(M, N, Val(UNIT), typeof(tup), ftup...)
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
  spa, spu = reassemble_tup(Args, args)
  T = eltype(spa)
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF * WS
  Nr = VectorizationBase.vrem(N, WS)
  m = 0
  # m, no remainder
  while m < M - WS + 1
    n = Nr # non factor of W remainder
    if n > 0
      let t = (spa, spu), ft = flatten_to_tup(t)
        BdivU_small_kern_u!(n, StaticInt(1), Val(UNIT), WS, typeof(t), ft...)
      end
    end
    while n < N - (WU - 1)
      ldiv_solve_W_u!(spa, spu, n, WS, UF, Val(UNIT))
      n += WU
    end
    while n != N
      ldiv_solve_W!(spa, spu, n, WS, Val(UNIT))
      n += W
    end
    m += W
    spa = gesp(spa, (W, StaticInt(0)))
  end
  # remainder on `m`
  if m < M
    let tup = (spa, spu), ftup = flatten_to_tup(tup)
      ldiv_remainder!(M, N, m, Nr, WS, Val(UNIT), typeof(tup), ftup...)
    end
  end
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
