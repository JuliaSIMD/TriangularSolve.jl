
struct Mat{T,ColMajor} <: AbstractMatrix{T}
  p::Ptr{T}
  x::Int
  M::Int
  N::Int
end
Base.size(A::Mat)::Tuple{Int,Int} = (A.M, A.N)::Tuple{Int,Int}
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
StaticArrayInterface.stride_rank(::Type{<:Mat{<:Any,true}}) =
  (static(1), static(2))
StaticArrayInterface.stride_rank(::Type{<:Mat{<:Any,false}}) =
  (static(2), static(1))
StaticArrayInterface.contiguous_batch_size(::Type{<:Mat}) = static(0)
StaticArrayInterface.dense_dims(::Type{<:Mat{<:Any,true}}) =
  (static(true), static(false))
StaticArrayInterface.dense_dims(::Type{<:Mat{<:Any,false}}) =
  (static(false), static(true))
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
  @inbounds begin
    @turbo warn_check_args = false for n in indices((C, B), 2),
      m in indices((C, A), 1)

      Cmn = zero(eltype(C))
      for k in indices((A, B), (2, 1))
        Cmn -= A[m, k] * B[k, n]
      end
      C[m, n] += Cmn
    end
  end
end
@inline function _schur_complement!(C::Mat, A::Mat, B::Mat, ::Val{true})
  # _turbo_! will not be inlined
  @tturbo warn_check_args = false for n in indices((C, B), 2),
    m in indices((C, A), 1)

    Cmn = zero(eltype(C))
    for k in indices((A, B), (2, 1))
      Cmn -= A[m, k] * B[k, n]
    end
    C[m, n] += Cmn
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

