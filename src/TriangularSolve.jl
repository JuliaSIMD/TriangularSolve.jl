module TriangularSolve
if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
  @eval Base.Experimental.@max_methods 1
end

using LayoutPointers:
  stridedpointer_preserve, AbstractStridedPointer, zero_offsets, StridedPointer
using VectorizationBase, LinearAlgebra #LoopVectorization
using VectorizationBase: vfnmadd_fast, AbstractMask, gesp
using CloseOpenIntervals: CloseOpen, SafeCloseOpen
using Static
using IfElse: ifelse
using LoopVectorization
using Polyester
using StaticArrayInterface

const LPtr{T} = Core.LLVMPtr{T,0}
# _lptr(x::Ptr{T}) where {T} = Base.bitcast(LPtr{T}, x)::LPtr{T}
# _ptr(x::LPtr{T}) where {T} = Base.bitcast(Ptr{T}, x)::Ptr{T}
# _lptr(x) = x
# _ptr(x) = x
const reassemble_tup = LoopVectorization.reassemble_tuple
const flatten_to_tup = LoopVectorization.flatten_to_tuple
# @inline reassemble_tup(::Type{T}, t) where {T} =
#   LoopVectorization.reassemble_tuple(T, map(_ptr, t))
# @inline flatten_to_tup(t) = map(_lptr, LoopVectorization.flatten_to_tuple(t))

#=
Av = ForwardDiff.value.(dA)
Ap = reinterpret(reshape,Float64,ForwardDiff.partials.(dA))
bv = ForwardDiff.value.(db)
bp = reinterpret(reshape, Float64, ForwardDiff.partials.(db))
F = lu!(Av);
c = F \ bv
cv0 = bp / F'
cv1 = reshape((c' * reshape((reshape(permutedims(Ap, (3,1,2)), (12,4))/F'), (4,12))), (3,4))
# cv1 = reshape((reshape(F\reshape(permutedims(Ap, (2,1,3)), (4,12)), (12,4)) * c), (4,3))
cv0 - cv1
=#

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

function m_thread_block_size(M, N, nthreads, ::Val{T}) where {T}
  W = VectorizationBase.pick_vector_width(T)
  nb = clamp(VectorizationBase.vdiv(M * N, StaticInt{256}() * W), 1, nthreads)
  min(M, VectorizationBase.vcld(M, nb * W) * W)
end

function block_size(::Val{T}) where {T}
  elements_l2 =
    (VectorizationBase.cache_size(StaticInt(2)) * StaticInt(19)) ÷
    (VectorizationBase.static_sizeof(T) * StaticInt(60))
  Static.floortostaticint(sqrt(elements_l2))
end


_nthreads() =
  min(Int(VectorizationBase.num_cores())::Int, Threads.nthreads()::Int)

_canonicalize(x) = signed(x)
_canonicalize(::StaticInt{N}) where {N} = StaticInt{N}()

include("schur_complement.jl")
include("rdivu.jl")
include("rdivl.jl")

ldiv!(A, B, ::Val = Val(true)) = LinearAlgebra.ldiv!(A, B)
ldiv!(Y, A, B, ::Val = Val(true)) = LinearAlgebra.ldiv!(Y, A, B)
rdiv!(A, B, ::Val = Val(true)) = LinearAlgebra.rdiv!(A, B)

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
