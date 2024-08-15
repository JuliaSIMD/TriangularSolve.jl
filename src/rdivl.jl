# rdivl and ldivu are basically equivalent, and names used somewhat interchangeably.

# Solve A / L
# [ A11 A12 A13    = [ B11 B12 B13    * [ L11
#   A21 A22 A23 ]      B21 B22 B23 ]      L21 L22
#                                         L31 L32 L33 ]
# 
# A and B are M x N, L is N x N
# A_{m,n} = \sum_{i=n}^N B_{m,i}L_{i,n}
# A_{m,n} = B_{m,n}L_{n,n} + \sum_{i=n+1}^N B_{m,i}L_{i,n}
# 
# B_{m,n} = (A_{m,n} - \sum_{i=n+1}^N B_{m,i}L_{i,n})/L_{n,n}
# 
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
  Asym = Vector{Symbol}(undef, N)
  Lsym = Matrix{Symbol}(undef, N, N)
  rett = Expr(:tuple)
  for n = 1:N
    A_n = Asym[n] = Symbol(:A_, n)
    push!(rett.args, A_n)
    push!(q.args, Expr(:(=), A_n, Expr(:call, getfield, :Ad, n)))
    for m = (UNIT ? n + 1 : n):N
      L_m_n = Lsym[m, n] = Symbol(:L_, n * N + m)
      push!(
        q.args,
        Expr(:(=), L_m_n, :(vload(spl, (noff + $(m - 1), noff + $(n - 1)))))
      )
    end
  end
  for n = N:-1:1
    A_n = Asym[n]
    for k = n+1:N
      push!(
        q.args,
        Expr(:(=), A_n, Expr(:call, vfnmadd_fast, Asym[k], Lsym[k, n], A_n))
      )
    end
    if !UNIT
      push!(
        q.args,
        Expr(:(=), A_n, Expr(:call, Base.FastMath.div_fast, A_n, Lsym[n, n]))
      )
    end
  end
  push!(q.args, Expr(:call, VecUnroll, rett))
  q
end
@generated function BdivL_small_kern!(
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
      Amn = :($Amn / vload(spl, $((z, z))))
    end
    quote
      $(Expr(:meta, :inline))
      mask = $(VectorizationBase.Mask{W})(_mask)
      spa, spl = reassemble_tup($Args, args)
      vstore!(spa, $Amn, $i, mask)
    end
  else
    unroll = Unroll{2,1,N,1,W,(-1 % UInt),1}((z, z))
    quote
      $(Expr(:meta, :inline))
      spa, spl = reassemble_tup($Args, args)
      mask = $(VectorizationBase.Mask{W})(_mask)
      Amn = vload(spa, $unroll, mask)
      vstore!(spa, solve_AL(Amn, spl, 0, $(Val(UNIT))), $unroll, mask)
    end
  end
end
@generated function BdivL_small_kern!(
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
      BdivL_small_kern!(static(n), mask, $WS, $(Val(UNIT)), $Args, args...)
  end
end
@generated function BdivL_small_kern_u!(
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
      Amn = :($Amn / vload(spl, $((z, z))))
    end
    quote
      $(Expr(:meta, :inline))
      spa, spl = reassemble_tup($Args, args)
      vstore!(spa, $Amn, $unroll)
    end
  else
    double_unroll =
      Unroll{2,1,N,1,W,zero(UInt),1}(Unroll{1,W,U,1,W,zero(UInt),1}((z, z)))
    quote
      $(Expr(:meta, :inline))
      spa, spl = reassemble_tup($Args, args)
      Amn = vload(spa, $double_unroll)
      vstore!(spa, solve_AL(Amn, spl, 0, $(Val(UNIT))), $double_unroll)
    end
  end
end
@generated function BdivL_small_kern_u!(
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
      BdivL_small_kern_u!(static(n), $su, $vu, $sw, $Args, args...)
  end
end

@generated function rdivl_solve_W_u!(
  spa,
  spl,
  n,
  ::StaticInt{W},
  ::StaticInt{U},
  ::Val{UNIT}
) where {W,U,UNIT}
  # n is num cols of `spa` to reduce
  z = static(0)
  error("not updated")
  quote
    $(Expr(:meta, :inline))
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
        C11_c = vfnmadd_fast(A11, vload(spl, (nk, n + (c - 1))), C11_c)
    end
    C11vu =
      solve_AU(VecUnroll((Base.Cartesian.@ntuple $W C11)), spl, n, $(Val(UNIT)))
    i = $(Unroll{2,1,W,1,W,zero(UInt),1})(
      $(Unroll{1,W,U,1,W,zero(UInt),1})(($z, n))
    )
    vstore!(spa, C11vu, i)
  end
end
@generated function rdivl_solve_W!(
  spa,
  spl,
  n,
  mask::AbstractMask{W},
  ::Val{UNIT}
) where {W,UNIT}
  z = static(0)
  error("not updated")
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
        C11_c = vfnmadd_fast(A11, vload(spl, (nk, n + (c - 1))), C11_c)
    end
    C11 = VecUnroll((Base.Cartesian.@ntuple $W C11))
    C11 = solve_AU(C11, spl, n, $(Val(UNIT)))
    i = $(Unroll{2,1,W,1,W,(-1 % UInt),1})(($z, n))
    vstore!(spa, C11, i, mask)
  end
end

@inline function rdivl_U!(
  spa::AbstractStridedPointer{T},
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
  m = 0
  z = StaticInt(0)
  if UF > 1
    while m < M - WU + 1
      n = Int(Nd * W)::Int
      if Nr > 0
        let t = (gesp(spa, (z, n)), gesp(spl, (n, n))), ft = flatten_to_tup(t)
          BdivL_small_kern_u!(Nr, UF, Val(UNIT), WS, typeof(t), ft...)
        end
      end
      for _ ∈ 1:Nd
        k = N - n
        n -= W
        rdivl_solve_W_u!(
          gesp(spa, (z, n)),
          gesp(spl, (n, n)),
          k,
          WS,
          UF,
          Val(UNIT)
        )
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
    n = Int(Nd * W)::Int
    if Nr > 0
      let t = (gesp(spa, (z, n)), gesp(spl, (n, n))),
        ft = flatten_to_tup(t),
        mask = getfield(mask, :u) % UInt32

        BdivL_small_kern!(n, mask, WS, Val(UNIT), typeof(t), ft...)
      end
    end
    for _ ∈ 1:Nd
      k = N - n
      n -= W
      rdivl_solve_W!(gesp(spa, (z, n)), gesp(spl, (n, n)), k, mask, Val(UNIT))
    end
    spa = gesp(spa, (WS, StaticInt(0)))
    m = ubm
  end
  nothing
end
function _ldivu_L!(
  M,
  N,
  ::Val{UNIT},
  ::Type{Args},
  args::Vararg{Any,K}
) where {UNIT,Args,K}
  error("not updated")
  spa, spl = reassemble_tup(Args, args)
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
      let t = (spa, spl), ft = flatten_to_tup(t)
        BdivU_small_kern_u!(n, StaticInt(1), Val(UNIT), WS, typeof(t), ft...)
      end
    end
    while n < N - (WU - 1)
      ldivl_solve_W_u!(spa, spl, n, WS, UF, Val(UNIT))
      n += WU
    end
    while n != N
      ldivl_solve_W!(spa, spl, n, WS, Val(UNIT))
      n += W
    end
    m += W
    spa = gesp(spa, (W, StaticInt(0)))
  end
  # remainder on `m`
  if m < M
    let tup = (spa, spl), ftup = flatten_to_tup(tup)
      ldivl_remainder!(M, N, m, Nr, WS, Val(UNIT), typeof(tup), ftup...)
    end
  end
  nothing
end
@inline function rdivl_U!(
  spa::AbstractStridedPointer{T,2,2},
  spl::AbstractStridedPointer{T,2,2},
  M,
  N,
  ::Val{UNIT}
) where {T,UNIT}
  tup = (spa, spl)
  ftup = flatten_to_tup(tup)
  _ldivu_L!(M, N, Val(UNIT), typeof(tup), ftup...)
end

# like rdivu, except we iterate backwards
function rdivl_block_N!(
  M,
  N,
  ::Val{UNIT},
  Bsize,
  ::Type{Args},
  args::Vararg{Any,K}
) where {K,Args,UNIT}
  spa, spl = reassemble_tup(Args, args)
  spa_base = spa
  T = eltype(spa)
  W = VectorizationBase.pick_vector_width(T)
  B_normalized =
    Bsize === nothing ?
    VectorizationBase.vcld(
      N,
      VectorizationBase.vcld(N, block_size(Val(T))) * W
    ) * W : Bsize

  Niter = VectorizationBase.vdiv(N, B_normalized)
  Nrem = N - B_normalized * Niter

  N_temp = Nrem != 0 ? Nrem : B_normalized
  Niter -= Nrem == 0
  n = Niter * B_normalized

  spa = gesp(spa, (StaticInt(0), n))
  spl = gesp(spl, (n, n))
  while true
    # println("Solve with N_temp = $N_temp and n = $n")
    rdivl_U!(spa, spl, M, N_temp, Val{UNIT}())
    n == 0 && break
    spa_prev = spa
    spa = gesp(spa, (StaticInt(0), -B_normalized))
    spl = gesp(spl, (StaticInt(0), -B_normalized))
    N_temp = B_normalized
    k = N - n
    n -= B_normalized
    schur_complement!(
      Mat(spa, M, B_normalized),
      Mat(spa_prev, M, k),
      Mat(spl, k, B_normalized),
      Val(false)
    )
    spl = gesp(spl, (-B_normalized, StaticInt(0)))
  end
end
function rdivl_block_MandN!(
  M,
  N,
  ::Val{UNIT},
  ::Type{Args},
  args::Vararg{Any,K}
) where {UNIT,Args,K}
  spa, spl = reassemble_tup(Args, args)
  T = eltype(spa)
  B = block_size(Val(T))
  W = VectorizationBase.pick_vector_width(T)
  XA = _contig_axis(spa)
  XA = _contig_axis(spl)
  WUF = XA == XA == 2 ? W : W * unroll_factor(W)
  B_m = VectorizationBase.vcld(M, VectorizationBase.vcld(M, B) * WUF) * WUF
  m = 0
  while m < M
    mu = m + B_m
    Mtemp = min(M, mu) - m
    let tup = (spa, spl), ftup = flatten_to_tup(tup)
      rdivl_block_N!(
        Mtemp,
        N,
        Val{UNIT}(),
        VectorizationBase.vcld(N, VectorizationBase.vcld(N, B) * W) * W,
        typeof(tup),
        ftup...
      )
    end
    spa = gesp(spa, (B_m, StaticInt{0}()))
    m = mu
  end
  nothing
end

function rdivl_dispatch!(A::AbstractMatrix{T}, U, ::Val{UNIT}) where {UNIT,T}
  _M, _N = size(A)
  M = _canonicalize(_M)
  N = _canonicalize(_N)
  ((N == 0) | (M == 0)) && return nothing
  _spa, spap = stridedpointer_preserve(A)
  _spl, splp = stridedpointer_preserve(U)
  spa = zero_offsets(_spa)
  spl = zero_offsets(_spl)
  GC.@preserve spap splp begin
    N <= block_size(Val(T)) && return rdivl_U!(spa, spl, M, N, Val(UNIT))
    let tup = (spa, spl), ftup = flatten_to_tup(tup)
      return rdivl_block_MandN!(M, N, Val(UNIT), typeof(tup), ftup...)
    end
  end
end

function rdiv!(
  A::AbstractMatrix{T},
  U::UpperTriangular{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(A, parent(U), Val(false))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UpperTriangular{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(copyto!(C, A), parent(U), Val(false))
  return C
end
function rdiv!(
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(A, parent(U), Val(true))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UnitUpperTriangular{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(copyto!(C, A), parent(U), Val(true))
  return C
end
function ldiv!(
  U::LowerTriangular{T},
  A::AbstractMatrix{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(A), transpose(parent(U)), Val(false))
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::LowerTriangular{T},
  A::AbstractMatrix{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(copyto!(C, A)), transpose(parent(U)), Val(false))
  return C
end
function ldiv!(
  U::UnitLowerTriangular{T},
  A::AbstractMatrix{T},
  ::Val{false}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(A), transpose(parent(U)), Val(true))
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::UnitLowerTriangular{T},
  A::AbstractMatrix{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(copyto!(C, A)), transpose(parent(U)), Val(true))
  return C
end
