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
  div = Base.FastMath.div_fast
  for n = 1:N
    A_n = Asym[n] = Symbol(:A_, n)
    push!(rett.args, A_n)
    push!(q.args, Expr(:(=), A_n, Expr(:call, getfield, :Ad, n)))
    for m = (UNIT ? n + 1 : n):N
      L_m_n = Lsym[m, n] = Symbol(:L_, n * N + m)
      push!(q.args, Expr(:(=), L_m_n, :(vload(spl, ($(m - 1), $(n - 1))))))
    end
    if !UNIT
      Lnn = Lsym[n, n]
      # latency is high, so start these early
      invLnn = Expr(:call, div, Expr(:call, one, Lnn), Lnn)
      push!(q.args, Expr(:(=), Lnn, invLnn))
    end
  end
  for n = N:-1:1
    A_n = Asym[n]
    if !UNIT
      push!(
        q.args,
        Expr(:(=), A_n, Expr(:call, Base.FastMath.mul_fast, A_n, Lsym[n, n]))
      )
    end
    for k = 1:n-1
      push!(
        q.args,
        Expr(
          :(=),
          Asym[k],
          Expr(:call, vfnmadd_fast, A_n, Lsym[n, k], Asym[k])
        )
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
      vstore!(spa, solve_AL(Amn, spl, $(Val(UNIT))), $unroll, mask)
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
      vstore!(spa, solve_AL(Amn, spl, $(Val(UNIT))), $double_unroll)
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
  quote
    $(Expr(:meta, :inline))
    i = $(Unroll{2,1,W,1,W,zero(UInt),1})(
      $(Unroll{1,W,U,1,W,zero(UInt),1})(($z, n))
    )
    C11 = VectorizationBase.data(vload(spa, i))
    Base.Cartesian.@nexprs $W c -> C11_c = $getfield(C11, c)
    for nk ∈ SafeCloseOpen(n) # nmuladd
      nkw = nk + $W
      A11 = vload(spa, $(Unroll{1,W,U,1,W,zero(UInt),1})(($z, nkw)))
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(spl, (nkw, c - 1)), C11_c)
    end
    # `solve_AL` solves A /= L, where
    # A is WU x W
    # L is W x W
    C11vu =
      solve_AL(VecUnroll((Base.Cartesian.@ntuple $W C11)), spl, $(Val(UNIT)))
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
  quote
    $(Expr(:meta, :inline))
    # here, we just want to load the vectors
    i = $(Unroll{2,1,W,1,W,(-1 % UInt),1})(($z, n))
    C11 = VectorizationBase.data(vload(spa, i, mask))
    Base.Cartesian.@nexprs $W c -> C11_c = $getfield(C11, c)
    for nk ∈ SafeCloseOpen(n) # nmuladd
      nkw = nk + $W
      A11 = vload(spa, ($(MM{W}(z)), nkw), mask)
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(spl, (nkw, n + (c - 1))), C11_c)
    end
    C11 = VecUnroll((Base.Cartesian.@ntuple $W C11))
    C11 = solve_AL(C11, spl, $(Val(UNIT)))
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
  finalmask = UInt32(getfield(_mask(WS, M), :u))
  maxmask = UInt32(getfield(VectorizationBase.max_mask(WS), :u))
  while m < M
    ubm = m + W
    nomaskiter = ubm < M
    mask = nomaskiter ? maxmask : finalmask
    n = Int(Nd * W)::Int
    if Nr > 0
      let t = (gesp(spa, (z, n)), gesp(spl, (n, n))), ft = flatten_to_tup(t)
        BdivL_small_kern!(Nr, mask, WS, Val(UNIT), typeof(t), ft...)
      end
    end
    for _ ∈ 1:Nd
      k = N - n
      n -= W
      rdivl_solve_W!(
        gesp(spa, (z, n)),
        gesp(spl, (n, n)),
        k,
        Mask{W}(mask),
        Val(UNIT)
      )
    end
    spa = gesp(spa, (WS, StaticInt(0)))
    m = ubm
  end
  nothing
end
@generated function ldivu_solve_W_u!(
  spa,
  spu,
  n,
  ::StaticInt{W},
  ::StaticInt{U},
  ::Val{UNIT}
) where {W,U,UNIT}
  z = static(0)
  error("not implemented")
  quote
    # $(Expr(:meta, :inline))
    # C = U \ A; U * C = A
    # A_{i,j} = U_{i,i}*C_{i,j} + \sum_{k=i+1}^{N}U_{i,k}C_{k,j}
    # C_{i,j} = U_{i,i} \ (A_{i,j} - \sum_{k=i+1}^{N}U_{i,k}C_{k,j})
    # The inputs here are transposed, as the library was formulated in terms of `rdiv!`,
    # so we have
    # C_{j,i} = (A_{j,i} - \sum_{k=i+1}^{N}C_{j,k}U_{k,i}) / L_{i,i}
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
@generated function ldivu_solve_W!(
  spa,
  spu,
  n,
  ::StaticInt{W},
  ::Val{UNIT}
) where {W,UNIT}
  z = static(0)
  error("not implemented")
  quote
    # $(Expr(:meta, :inline))
    # Like `ldivl_solve_W_u!`, except no unrolling, just a `W`x`W` block
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

@generated function _ldivu_remainder!(
  spa,
  spu,
  N,
  Nr,
  ::StaticInt{W},
  ::Val{UNIT},
  ::StaticInt{r}
) where {W,UNIT,r}
  error("not updated")
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

          BdivL_small_kern!(n, mask, $WS, $(Val(UNIT)), typeof(t), ft...)
        end
      end
      # while n < N - $(W * U - 1)
      #   ldivu_solve_W_u!(spa, spa, spu, n, $WS, $US, Val(UNIT), Val(r))
      #   n += $(W * U)
      # end
      while n != N
        ldivu_solve_W!(spa, spu, n, $WS, $(Val(UNIT)), $(StaticInt(r)))
        n += $W
      end
    end
  end
end
@generated function ldivu_remainder!(
  M,
  N,
  m,
  Nr,
  ::StaticInt{W},
  # ::Val{U},
  ::Val{UNIT},
  ::Type{Args},
  args::Vararg{Any,K}
) where {W,UNIT,Args,K}
  error("not updated")
  WS = static(W)
  # US = static(U)
  if W == 2
    quote
      $(Expr(:meta, :inline))
      spa, spu = reassemble_tup(Args, args)
      _ldivu_remainder!(spa, spu, N, Nr, $WS, $(Val(UNIT)), $(static(1)))
      nothing
    end
  elseif W == 8
    quote
      # $(Expr(:meta, :inline))
      spa, spu = reassemble_tup(Args, args)
      if m == M - 1
        _ldivu_remainder!(
          spa,
          spu,
          N,
          Nr,
          static(8),
          $(Val(UNIT)),
          StaticInt(1)
        )
      else
        if m == M - 2
          _ldivu_remainder!(
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
            _ldivu_remainder!(
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
              _ldivu_remainder!(
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
                _ldivu_remainder!(
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
                  _ldivu_remainder!(
                    spa,
                    spu,
                    N,
                    Nr,
                    static(8),
                    $(Val(UNIT)),
                    StaticInt(6)
                  )
                else
                  _ldivu_remainder!(
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
        _ldivu_remainder!(spa, spu, N, Nr, $WS, $(Val(UNIT)), static(w))
      nothing
    end
  end
end
# B_{n,m} = (A_{n,m} - \sum_{i=n+1}^N U_{n,i}B_{i,m})/U_{n,n}
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
# for ldiv, we unroll over `n`
  Nd, Nr = VectorizationBase.vdivrem(N, WS)
  z = StaticInt(0)
  m = 0
  # m, no remainder
  while m < M - WS + 1
    n = Int(Nd * W)::Int
    if Nr > 0
      let t = (gesp(spa, (n, z)), gesp(spl, (n, n))), ft = flatten_to_tup(t)
        BdivL_small_kern_u!(Nr, StaticInt(1), Val(UNIT), WS, typeof(t), ft...)
      end
    end
    for _ ∈ 1:Nd
      k = N - n
      n -= W
      ldivu_solve_W_u!(
        gesp(spa, (n, z)),
        gesp(spl, (n, n)),
        k,
        WS,
        UF,
        Val(UNIT)
      )
    end
    while n < N - (WU - 1)
      ldivu_solve_W_u!(spa, spl, n, WS, UF, Val(UNIT))
      n += WU
    end

    n = Nr # non factor of W remainder
    if n > 0
      let t = (spa, spl), ft = flatten_to_tup(t)
        BdivL_small_kern_u!(n, StaticInt(1), Val(UNIT), WS, typeof(t), ft...)
      end
    end
    while n < N - (WU - 1)
      ldivu_solve_W_u!(spa, spl, n, WS, UF, Val(UNIT))
      n += WU
    end
    while n != N
      ldivu_solve_W!(spa, spl, n, WS, Val(UNIT))
      n += W
    end
    m += W
    spa = gesp(spa, (W, z))
  end
  # remainder on `m`
  if m < M
    let tup = (spa, spl), ftup = flatten_to_tup(tup)
      ldivu_remainder!(M, N, m, Nr, WS, Val(UNIT), typeof(tup), ftup...)
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
  U::LowerTriangular{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(A, parent(U), Val(false))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::LowerTriangular{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(copyto!(C, A), parent(U), Val(false))
  return C
end
function rdiv!(
  A::AbstractMatrix{T},
  U::UnitLowerTriangular{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(A, parent(U), Val(true))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UnitLowerTriangular{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(copyto!(C, A), parent(U), Val(true))
  return C
end
function ldiv!(
  U::UpperTriangular{T},
  A::AbstractMatrix{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(A), transpose(parent(U)), Val(false))
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::UpperTriangular{T},
  A::AbstractMatrix{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(copyto!(C, A)), transpose(parent(U)), Val(false))
  return C
end
function ldiv!(
  U::UnitUpperTriangular{T},
  A::AbstractMatrix{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(A), transpose(parent(U)), Val(true))
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  A::AbstractMatrix{T}
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(copyto!(C, A)), transpose(parent(U)), Val(true))
  return C
end

function rdiv!(
  A::AbstractMatrix{T},
  U::LowerTriangular{T},
  ::Val
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(A, parent(U), Val(false))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::LowerTriangular{T},
  ::Val
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(copyto!(C, A), parent(U), Val(false))
  return C
end
function rdiv!(
  A::AbstractMatrix{T},
  U::UnitLowerTriangular{T},
  ::Val
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(A, parent(U), Val(true))
  return A
end
function rdiv!(
  C::AbstractMatrix{T},
  A::AbstractMatrix{T},
  U::UnitLowerTriangular{T},
  ::Val
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(copyto!(C, A), parent(U), Val(true))
  return C
end
function ldiv!(
  U::UpperTriangular{T},
  A::AbstractMatrix{T},
  ::Val
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(A), transpose(parent(U)), Val(false))
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::UpperTriangular{T},
  A::AbstractMatrix{T},
  ::Val
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(copyto!(C, A)), transpose(parent(U)), Val(false))
  return C
end
function ldiv!(
  U::UnitUpperTriangular{T},
  A::AbstractMatrix{T},
  ::Val
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(A), transpose(parent(U)), Val(true))
  return A
end
function ldiv!(
  C::AbstractMatrix{T},
  U::UnitUpperTriangular{T},
  A::AbstractMatrix{T},
  ::Val
) where {T<:Union{Float32,Float64}}
  rdivl_dispatch!(transpose(copyto!(C, A)), transpose(parent(U)), Val(true))
  return C
end
