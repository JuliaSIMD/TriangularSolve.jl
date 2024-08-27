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
      AdL = Expr(:call, Base.FastMath.mul_fast, A_n, Lsym[n, n])
      push!(q.args, Expr(:(=), A_n, AdL))
    end
    for k = 1:n-1
      AdL = Expr(:call, vfnmadd_fast, A_n, Lsym[n, k], Asym[k])
      push!(q.args, Expr(:(=), Asym[k], AdL))
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
  i = Unroll{2,1,W,1,W,zero(UInt),1}(Unroll{1,W,U,1,W,zero(UInt),1}((z, z)))
  quote
    $(Expr(:meta, :inline))
    C11 = VectorizationBase.data(vload(spa, $i))
    Base.Cartesian.@nexprs $W c -> C11_c = $getfield(C11, c)
    for nk ∈ SafeCloseOpen(n) # nmuladd
      nkw = nk + $W
      A11 = vload(spa, $(Unroll{1,W,U,1,W,zero(UInt),1})(($z, nkw)))
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(spl, (nkw, static(c - 1))), C11_c)
    end
    # `solve_AL` solves A /= L, where
    # A is WU x W
    # L is W x W
    C11vu =
      solve_AL(VecUnroll((Base.Cartesian.@ntuple $W C11)), spl, $(Val(UNIT)))
    vstore!(spa, C11vu, $i)
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
  i = Unroll{2,1,W,1,W,(-1 % UInt),1}((z, z))
  quote
    $(Expr(:meta, :inline))
    # here, we just want to load the vectors
    C11 = VectorizationBase.data(vload(spa, $i, mask))
    Base.Cartesian.@nexprs $W c -> C11_c = $getfield(C11, c)
    for nk ∈ SafeCloseOpen(n) # nmuladd
      nkw = nk + $W
      A11 = vload(spa, ($(MM{W}(z)), nkw), mask)
      Base.Cartesian.@nexprs $W c ->
        C11_c = vfnmadd_fast(A11, vload(spl, (nkw, static(c - 1))), C11_c) #=n +=#
    end
    C11 = VecUnroll((Base.Cartesian.@ntuple $W C11))
    C11 = solve_AL(C11, spl, $(Val(UNIT)))
    vstore!(spa, C11, $i, mask)
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
  spl,
  n,
  ::StaticInt{W},
  ::StaticInt{U},
  ::Val{UNIT}
) where {W,U,UNIT}
  z = static(0)
  # Actually a row-major rdivl
  # B_{m,n} = (A_{m,n} - \sum_{i=n+1}^N B_{m,i}L_{i,n})/L_{n,n}
  Aind = Unroll{1,1,W,2,W,zero(UInt),1}(Unroll{2,W,U,2,W,zero(UInt),1}((z, z)))
  q = quote
    # $(Expr(:meta, :inline))
    # 
    A11 = getfield(vload(spa, $Aind), :data)
    # The `W` rows
    Base.Cartesian.@nexprs $W c -> A11_c = getfield(A11, c)
    # compute
    # A_{j,i} - \sum_{k=1}^{i-1}U_{k,i}C_{j,k})
    # Each iter:
    # A_{j+[0,W), i+[0,W*U)} -= C_{j+[0,W),k}*U_{k,i+[0,W*U)}
    for nk ∈ SafeCloseOpen(n) # nmuladd
      nkw = nk + $(W * U)
      U_ki = vload(spl, $(Unroll{2,W,U,2,W,zero(UInt),1})((nkw, $z)))
      Base.Cartesian.@nexprs $W c ->
        A11_c = vfnmadd_fast(U_ki, vload(spa, (static(c - 1), nkw)), A11_c)
    end
  end
  # solve AL wants:
  # outer unroll are `W` columns
  # Inner unroll are `W` rows (U simd vecs)
  # So, we'll use `U = 1`, and transpose blocks
  # We then have column-major multiplies
  # This is needed, due to columns being independent, but rows dependent
  # Unroll across the `u` blocks, we move up
  Xu = Vector{Symbol}(undef, W)
  Csym = Vector{Symbol}(undef, U)
  for u = 1:U
    # X_u are future
    X_u = Symbol(:X_, u)
    push!(
      q.args,
      :(
        $X_u = getfield(
          VectorizationBase.transpose_vecunroll(
            VecUnroll(
              Base.Cartesian.@ntuple $W w ->
                getfield(getfield(A11_w, :data), $(U + 1 - u))
            )
          ),
          :data
        )
      )
    )
    # push!(q.args, :(println($X_u)))
    for c = 1:W
      X_u_c = Xu[c] = Symbol(:X_, u, :_, c)
      push!(q.args, Expr(:(=), X_u_c, Expr(:call, getfield, X_u, c)))
    end
    # take A[(U-u+1)*W,u*W), [0,W)]
    for j = 1:u-1 # iter over all blocks ordered after
      for k = 1:W # reduction dimension, reverse order
        for c = 1:W # columns of C
          urow = ((W - k) + ((U - j) * W))
          ucol = ((c - 1) + ((U - u) * W))
          # push!(q.args, Expr(:call, println, "Row = $urow; Col = $ucol"))
          Uexpr = :(vload(spl, ($urow, $ucol)))
          X_u_c = Xu[c]
          C_j_k = Symbol(:C_, j, :_, W + 1 - k)
          Xucexpr = Expr(:call, vfnmadd_fast, C_j_k, Uexpr, X_u_c)
          push!(q.args, Expr(:(=), X_u_c, Xucexpr))
        end
      end
    end
    o = (U - u) * W
    sp = Expr(:call, gesp, :spl, (o, o))
    Xut = Expr(:tuple)
    for c = 1:W
      push!(Xut.args, Xu[c])
    end
    Xuvu = Expr(:call, VecUnroll, Xut)
    C_u = Csym[u] = Symbol(:C_, u)
    push!(q.args, Expr(:(=), C_u, :(solve_AL($Xuvu, $sp, $(Val(UNIT))))))
    Cud = Symbol(C_u, :data)
    push!(q.args, Expr(:(=), Cud, Expr(:call, getfield, C_u, QuoteNode(:data))))
    for c = 1:W
      push!(
        q.args,
        Expr(:(=), Symbol(:C_, u, :_, c), Expr(:call, getfield, Cud, c))
      )
    end
  end
  for u = 1:U
    # u = 1 is last, first processed (reverse order)
    ui = Unroll{2,1,W,1,W,zero(UInt),1}((z, (U - u) * W))
    C_u = Csym[u]
    push!(q.args, :(vstore!(spa, $C_u, $ui)))
  end
  return q
end
@generated function ldivu_solve_W!(
  spa,
  spl,
  n,
  ::StaticInt{W},
  ::Val{UNIT},
  ::StaticInt{R}
) where {W,UNIT,R}
  # B_{n,m} = (A_{n,m} - \sum_{i=n+1}^N U_{n,i}B_{i,m})/U_{n,n}
  R <= 1 && throw("Remainder of `<= 1` shouldn't be called, but had $R.")
  R > W && throw("Reaminderof `> $W` shouldn't be called, but had $R.")
  z = static(0)
  Aind = Unroll{1,1,R,2,W,zero(UInt),1}((z, z))
  q = quote
    # $(Expr(:meta, :inline))
    # Like `ldivl_solve_W_u!`, except no unrolling, just a `W`x`W` block
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
    A11 = getfield(vload(spa, $Aind), :data)
    # The `W` rows
    Base.Cartesian.@nexprs $R r -> A11_r = getfield(A11, r)
    # compute
    # A_{j,i} - \sum_{k=1}^{i-1}U_{k,i}C_{j,k})
    # Each iter:
    # A_{j+[0,W), i+[0,W*U)} -= C_{j+[0,W),k}*U_{k,i+[0,W*U)}
    for nk ∈ SafeCloseOpen(n) # nmuladd
      nkw = nk + $W
      U_ki = vload(spl, (nkw, $(MM{W}(z))))
      Base.Cartesian.@nexprs $R r ->
        A11_r = vfnmadd_fast(U_ki, vload(spa, (static(r - 1), nkw)), A11_r)
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
    C_u = solve_AL(X, spl, $(Val(UNIT)))
  end
  push!(q.args, q2)
  q3 = if R == Wpad
    quote
      i = $(Unroll{2,1,W,1,Wpad,zero(UInt),1}((z, z)))
      vstore!(spa, C_u, i)
    end
  else
    quote
      mask = _mask($(static(Wpad)), $(static(R)))
      i = $(Unroll{2,1,W,1,Wpad,(-1 % UInt),1}((z, z)))
      vstore!(spa, C_u, i, mask)
    end
  end
  push!(q.args, q3)
  return q
end

@generated function _ldivu_remainder!(
  spa,
  spl,
  N,
  Nr,
  ::StaticInt{W},
  ::Val{UNIT},
  ::StaticInt{r}
) where {W,UNIT,r}
  # B_{n,m} = (A_{n,m} - \sum_{i=n+1}^N U_{n,i}B_{i,m})/U_{n,n}
  r <= 0 && throw("Remainder of `<= 0` shouldn't be called, but had $r.")
  r >= W && throw("Reaminderof `>= $W` shouldn't be called, but had $r.")
  z = static(0)
  if r == 1
    # one iter remaining
    sub = Base.FastMath.sub_fast
    mul = Base.FastMath.mul_fast
    div = Base.FastMath.div_fast
    vlxj = :(vload(spa, ($z, j)))
    if UNIT
      vlxj = :(xj = $vlxj)
    else
      vlxj = quote
        xj = $div($vlxj, vload(spl, (j, j)))
        vstore!(spa, xj, ($z, j))
      end
    end
    quote
      $(Expr(:meta, :inline))
      j = N
      while true
        j -= 1
        $vlxj
        i = j
        while i > 0
          i -= 1
          xi = vload(spa, ($z, i))
          Uji = vload(spl, (j, i))
          vstore!(spa, $sub(xi, $mul(xj, Uji)), ($z, i))
        end
        j == 0 && break
      end
    end
  else
    # `r` iters remaining, r > 1
    WS = static(W)
    quote
      $(Expr(:meta, :inline))
      mask = $(getfield(_mask(WS, r), :u) % UInt32)
      n = N - Nr
      if Nr > 0
        let t = (gesp(spa, ($z, n)), gesp(spl, (n, n))), ft = flatten_to_tup(t)
          BdivL_small_kern!(Nr, mask, $WS, $(Val(UNIT)), typeof(t), ft...)
        end
      end
      # non-U, order first as matmul kern is smaller than optimal
      while n != 0
        k = N - n
        n -= W
        ldivu_solve_W!(
          gesp(spa, ($z, n)),
          gesp(spl, (n, n)),
          k,
          $WS,
          Val(UNIT),
          $(StaticInt(r))
        )
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
  WS = static(W)
  if W == 2
    quote
      $(Expr(:meta, :inline))
      spa, spl = reassemble_tup(Args, args)
      _ldivu_remainder!(spa, spl, N, Nrr, Nru, $WS, $(Val(UNIT)), $(static(1)))
      nothing
    end
  elseif W == 8
    s8 = StaticInt(8)
    quote
      # $(Expr(:meta, :inline))
      spa, spl = reassemble_tup(Args, args)
      if m == M - 1
        _ldivu_remainder!(spa, spl, N, Nr, $s8, $(Val(UNIT)), $(StaticInt(1)))
      else
        if m == M - 2
          _ldivu_remainder!(spa, spl, N, Nr, $s8, $(Val(UNIT)), $(StaticInt(2)))
        else
          if m == M - 3
            _ldivu_remainder!(
              spa,
              spl,
              N,
              Nr,
              $s8,
              $(Val(UNIT)),
              $(StaticInt(3))
            )
          else
            if m == M - 4
              _ldivu_remainder!(
                spa,
                spl,
                N,
                Nr,
                $s8,
                $(Val(UNIT)),
                $(StaticInt(4))
              )
            else
              if m == M - 5
                _ldivu_remainder!(
                  spa,
                  spl,
                  N,
                  Nr,
                  $s8,
                  $(Val(UNIT)),
                  $(StaticInt(5))
                )
              else
                if m == M - 6
                  _ldivu_remainder!(
                    spa,
                    spl,
                    N,
                    Nr,
                    $s8,
                    $(Val(UNIT)),
                    $(StaticInt(6))
                  )
                else
                  _ldivu_remainder!(
                    spa,
                    spl,
                    N,
                    Nr,
                    $s8,
                    $(Val(UNIT)),
                    $(StaticInt(7))
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
      spa, spl = reassemble_tup(Args, args)
      Base.Cartesian.@nif $(W - 1) w -> m == M - w w ->
        _ldivu_remainder!(spa, spl, N, Nr, $WS, $(Val(UNIT)), static(w))
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
  # B_{n,m} = (A_{n,m} - \sum_{i=n+1}^N U_{n,i}B_{i,m})/U_{n,n}
  spa, spl = reassemble_tup(Args, args)
  T = eltype(spa)
  WS = pick_vector_width(T)
  W = Int(WS)
  UF = unroll_factor(WS)
  WU = UF * WS
  # for ldiv, we unroll over `n`
  Ndu, Nru = VectorizationBase.vdivrem(N, WU)
  Ndr, Nrr = VectorizationBase.vdivrem(Nru, WS)
  z = StaticInt(0)
  nstart = Int(Ndu * WU)::Int + Int(Ndr * W)::Int
  m = 0
  # compute = false
  compute = true
  # m, no remainder
  while m < M - WS + 1
    n::Int = nstart
    if Nrr > 0
      let t = (gesp(spa, (z, n)), gesp(spl, (n, n))), ft = flatten_to_tup(t)
        compute && BdivL_small_kern_u!(
          Nrr,
          StaticInt(1),
          Val(UNIT),
          WS,
          typeof(t),
          ft...
        )
      end
    end
    # non-U, order first as matmul kern is smaller than optimal
    for _ ∈ 1:Ndr
      k = N - n
      n -= W
      compute && ldivu_solve_W!(
        gesp(spa, (z, n)),
        gesp(spl, (n, n)),
        k,
        WS,
        Val(UNIT),
        WS
      )
    end
    # for _ ∈ 1:Ndu
    while n != 0
      k = N - n
      n -= Int(WU)
      compute && ldivu_solve_W_u!(
        gesp(spa, (z, n)),
        gesp(spl, (n, n)),
        k,
        WS,
        UF,
        Val(UNIT)
      )
    end
    m += W
    spa = gesp(spa, (W, StaticInt(0)))
  end
  # remainder on `m`
  if m < M
    let tup = (spa, spl), ftup = flatten_to_tup(tup)
      compute &&
        ldivu_remainder!(M, N, m, Nrr, WS, Val(UNIT), typeof(tup), ftup...)
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
