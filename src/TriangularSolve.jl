module TriangularSolve

using VectorizationBase, LinearAlgebra #LoopVectorization
using VectorizationBase: vfnmadd_fast, AbstractStridedPointer, AbstractMask, stridedpointer_preserve

@generated function solve_AU(A::VecUnroll{Nm1}, spu::AbstractStridedPointer, no) where {Nm1}
  N = Nm1 + 1
  quote
    $(Expr(:meta,:inline))
    Ad = VectorizationBase.data(A)
    Base.Cartesian.@nexprs $N n -> begin
      A_n = Ad[n]
      Base.Cartesian.@nexprs n m -> begin
        U_m_n = vload(spu, (no+m,no+n))
      end
    end
    Base.Cartesian.@nexprs $N n -> begin
      Base.Cartesian.@nexprs n-1 k -> begin
        A_n = Base.FastMath.sub_fast(A_n, Base.FastMath.mul_fast(A_k, U_k_n))
      end
      A_n = Base.FastMath.div_fast(A_n, U_n_n)
    end
    VecUnroll(Base.Cartesian.@ntuple $N A)
  end
end

# @generated function nmuladd(A::VecUnroll{Nm1},B::AbstractMatrix,C::VecUnroll{Nm1}) where {Nm1}
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
#         C_n = Base.FastMath.sub_fast(C_n, Base.FastMath.mul_fast(A_k, @inbounds(B[k,n])))
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


function BdivU_small_kern!(spa::AbstractStridedPointer{T}, spb::AbstractStridedPointer{T}, spu::AbstractStridedPointer{T}, m, M, N) where {T}
  W = VectorizationBase.pick_vector_width(T)
  while m < M
    ml = m+1
    mu = m+W
    maskiter = mu > M
    mask = maskiter ? VectorizationBase.mask(W, M) : VectorizationBase.max_mask(W)
    for n ∈ 1:N
      Amn = vload(spb, (MM(W, ml),n), mask)
      for k ∈ 1:n-1
        Amn = vfnmadd_fast(vload(spa, (MM(W, ml),k), mask), vload(spu, (k,n)), Amn)
      end
      vstore!(spa, Amn / vload(spu, (n,n)), (MM(W, ml),n), mask)
    end    
    m = mu
  end
end
function BdivU_small_kern_u3!(spa::AbstractStridedPointer{T}, spb::AbstractStridedPointer{T}, spu::AbstractStridedPointer{T}, m, N) where {T}
  WS = VectorizationBase.pick_vector_width(T)
  W = Int(WS)
  for n ∈ 1:N
    Amn = vload(spb, Unroll{1,W,3,1,W,0x0000000000000000,1}((m+1,n)))
    for k ∈ 1:n-1
      Amk = vload(spa, Unroll{1,W,3,1,W,0x0000000000000000,1}((m+1,k)))
      Amn = vfnmadd_fast(Amk, vload(spu, (k,n)), Amn)
    end
    vstore!(spa, Amn / vload(spu, (n,n)), Unroll{1,W,3,1,W,0x0000000000000000,1}((m+1,n)))
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
function nmuladd!(C,A,B,D)
  @turbo for n ∈ axes(C,2), m ∈ axes(C,1)
    Cmn = D[m,n]
    for k ∈ axes(B,1)
      Cmn -= A[m,k]*B[k,n]
    end
    C[m,n] = Cmn
  end
end

@generated function rdiv_solve_W_u3!(spc, spa, spu, lbm, n, ubn, ::StaticInt{W}) where {W}
  quote
    $(Expr(:meta,:inline))
    lbn = n + 1
    Nrange = lbn:ubn
    # here, we just want to load the vectors
    C11 = VectorizationBase.data(vload(spa, Unroll{2,1,$W,1,$W,0x0000000000000000,1}(Unroll{1,$W,3,1,$W,0x0000000000000000,1}((lbm,lbn)))))
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    for nk ∈ 1:n # nmuladd
      A11 = vload(spc, Unroll{1,$W,3,1,$W,0x0000000000000000,1}((lbm,nk)))
      Base.Cartesian.@nexprs $W c -> C11_c = vfnmadd_fast(A11, vload(spu, (nk,n+c)), C11_c)
    end
    C11vu = solve_AU(VecUnroll((Base.Cartesian.@ntuple $W C11)), spu, n)
    vstore!(spc, C11vu, Unroll{2,1,$W,1,$W,0x0000000000000000,1}(Unroll{1,$W,3,1,$W,0x0000000000000000,1}((lbm,lbn))))
  end
end
@generated function rdiv_solve_W_u1!(spc, spa, spu, lbm, n, ubn, mask::AbstractMask{W}) where {W}
  quote
    $(Expr(:meta,:inline))
    lbn = n + 1
    Nrange = lbn:ubn
    # here, we just want to load the vectors
    C11 = VectorizationBase.data(vload(spa, Unroll{2,1,$W,1,$W,0xffffffffffffffff,1}((lbm,lbn)), mask))
    Base.Cartesian.@nexprs $W c -> C11_c = C11[c]
    for nk ∈ 1:n # nmuladd
      A11 = vload(spc, (MM{$W}(lbm),nk), mask)
      Base.Cartesian.@nexprs $W c -> C11_c =  vfnmadd_fast(A11, vload(spu, (nk,n+c)), C11_c)
    end
    C11 = VecUnroll((Base.Cartesian.@ntuple $W C11))
    C11 = solve_AU(C11, spu, n)
    vstore!(spc, C11, Unroll{2,1,$W,1,$W,0xffffffffffffffff,1}((lbm,lbn)), mask)
  end
end

rdiv!(A, U::UpperTriangular) = rdiv_U!(A, A, parent(U))
rdiv!(C, A, U::UpperTriangular) = rdiv_U!(C, A, parent(U))
function rdiv_U!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, U::AbstractMatrix{T}) where {T}
  M, N = size(A)
  WS = pick_vector_width(T)
  W = Int(WS)
  W3 = StaticInt(3)*WS
  if VectorizationBase.register_count() > WS*StaticInt(3)
    # We have enough registers to use the 3W functions.
    Md, Mr = VectorizationBase.vdivrem(M, W3)
    Nd, Nr = VectorizationBase.vdivrem(N, W3)
  else
    # We do not have enough registers.
    Md = Nd = 0;
    Mr = M;
    Nr = N;
  end
  Mrd8, Mrr8 = VectorizationBase.vdivrem(Mr, WS)
  Nrd8, Nrr8 = VectorizationBase.vdivrem(Nr, WS)
  Nrd8main = Nrd8+Nd*3
  spa, spap = stridedpointer_preserve(A)
  spc, spcp = stridedpointer_preserve(C)
  spu, spup = stridedpointer_preserve(U)
  GC.@preserve spap spcp spup begin
    m = 0
    for _ ∈ 1:Md
      lbm = m+1
      ubm = m+W3
      Mrange = lbm:ubm
      
      n = 0
      if Nrr8 > 0
        n = Nrr8
        BdivU_small_kern_u3!(spc, spa, spu, m, Nrr8)
        # BdivU_small!(view(C,Mrange,1:n),view(A,Mrange,1:n),view(U,1:n,1:n))
      end
      for i ∈ 1:Nrd8main
        ubn = n+W
        rdiv_solve_W_u3!(spc, spa, spu, lbm, n, ubn, WS)
        n = ubn
      end
      #   for _ ∈ 1:Nd
      #     lbn = n+1
      #     ubn = n+W3
      #     Nrange = lbn:ubn
      #     if n == 0
      #       solve_3Wx3W!(view(C, Mrange, Nrange), view(A, Mrange, Nrange), view(U, Nrange, Nrange))
      #     else
      #       Cview = view(C, Mrange, Nrange)
      #       nmuladd!(Cview, view(C, Mrange, 1:n), view(U, 1:n, Nrange), view(A, Mrange, Nrange))
      #       solve_3Wx3W!(Cview, Cview, view(U, Nrange, Nrange))
      #     end
      #     n = ubn
      #   end
      m = ubm
    end
    for _ ∈ 1:Mrd8+(Mrr8>0)
      lbm = m+1
      ubm = m+W
      maskiter = ubm > M
      mask = maskiter ? VectorizationBase.mask(WS,Mrr8) : VectorizationBase.max_mask(WS)
      ubm = min(M,ubm)
      Mrange = lbm:ubm
      n = 0
      if Nrr8 > 0
        n = Nrr8
        BdivU_small_kern!(spc, spa, spu, m, ubm, n)
        # BdivU_small!(view(C,Mrange,1:n),view(A,Mrange,1:n),view(U,1:n,1:n))
      end
      for i ∈ 1:Nrd8main
        ubn = n+W
        rdiv_solve_W_u1!(spc, spa, spu, lbm, n, ubn, mask)
        n = ubn
      end
      # for _ ∈ 1:Nd
      #   lbn = n+1
      #   ubn = n+W3
      #   Nrange = lbn:ubn
      #   if n == 0
      #     solve_Wx3W!(spc, spa, view(U, Nrange, Nrange), lbm, lbn, mask)
      #   else
      #     nmuladd!(view(C, Mrange, Nrange), view(C, Mrange, 1:n), view(U, 1:n, Nrange), view(A, Mrange, Nrange))
      #     solve_Wx3W!(spc, spc, view(U, Nrange, Nrange), lbm, lbn, mask)
      #   end
      #   n = ubn
      # end
      m = ubm
    end
  end
  C
end




end
