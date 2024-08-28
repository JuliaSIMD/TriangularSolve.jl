# TriangularSolve

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSIMD.github.io/TriangularSolve.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSIMD.github.io/TriangularSolve.jl/dev)
[![Build Status](https://github.com/JuliaSIMD/TriangularSolve.jl/workflows/CI/badge.svg)](https://github.com/JuliaSIMD/TriangularSolve.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSIMD/TriangularSolve.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSIMD/TriangularSolve.jl)


Performs some triangular solves. For example:
```julia
julia> using TriangularSolve, LinearAlgebra, MKL;

julia> BLAS.set_num_threads(1)

julia> BLAS.get_config().loaded_libs
1-element Vector{LinearAlgebra.BLAS.LBTLibraryInfo}:
 LBTLibraryInfo(libmkl_rt.so, ilp64)

julia> N = 100;

julia> A = rand(N,N); B = rand(N,N); C = similar(A);

julia> @benchmark TriangularSolve.rdiv!(copyto!($C, $A), UpperTriangular($B), Val(false)) # false means single threaded
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  12.106 μs … 108.141 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     12.841 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   12.862 μs ±   1.585 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

   █                                                            
  ▃█▆▂▁▁▁▁▂▃▇▄▂▂▇▄▂▂▂▃▃▃▂▂▂▃▃▂▂▂▂▁▁▁▁▂▂▁▁▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▂▂▂▂▂ ▂
  12.1 μs         Histogram: frequency by time         16.6 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark rdiv!(copyto!($C, $A), UpperTriangular($B))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  13.483 μs …  31.165 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     13.698 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   13.816 μs ± 791.937 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

   ▃▆██▅         ▃▄▃                                           ▂
  ▆██████▄▃▁▁▃▁▁▇████▃▃▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▄▁▄▆▇▆▆ █
  13.5 μs       Histogram: log(frequency) by time      16.8 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark TriangularSolve.rdiv!(copyto!($C, $A), LowerTriangular($B), Val(false)) # false means single threaded
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  12.990 μs …  35.762 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     13.094 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   13.170 μs ± 636.637 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▄██▅          ▂▂                                             ▂
  ████▆▄▃▁▁▁▁▁▁▄███▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▄▅▄▅▅▇ █
  13 μs         Histogram: log(frequency) by time        16 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark rdiv!(copyto!($C, $A), LowerTriangular($B))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  14.119 μs …  32.767 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     14.321 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   14.372 μs ± 461.912 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

    ▂█▆                                                         
  ▂▄███▇▃▂▂▂▁▁▂▁▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▂▂▂▂▂ ▂
  14.1 μs         Histogram: frequency by time         17.3 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark TriangularSolve.ldiv!(LowerTriangular($B), copyto!($C, $A), Val(false)) # false means single threaded
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  16.579 μs …  38.431 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     16.728 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   16.812 μs ± 640.766 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▁▅██▅              ▁▂▁                                       ▂
  █████▇▃▁▁▁▁▁▁▁▁▁▁▁▄███▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▃▅▃▅▇█ █
  16.6 μs       Histogram: log(frequency) by time      19.6 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark ldiv!(LowerTriangular($B), copyto!($C, $A))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  15.763 μs … 38.550 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     15.936 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   16.351 μs ±  1.103 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▃▇█▆         ▁▅▆▄   ▂▄▅▃    ▁▂▁                             ▂
  █████▄▄▁▁▁▃▁▄████▆▁▅████▆▁▃▆███▅▃▁▁▁▁▁▁▁▁▁▁▁▁▄▃▄▅▆▇█▆▅▆▅▅▆▄ █
  15.8 μs      Histogram: log(frequency) by time      19.5 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark TriangularSolve.ldiv!(UpperTriangular($B), copyto!($C, $A), Val(false)) # false means single threaded
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  15.921 μs …  34.610 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     16.180 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   16.290 μs ± 733.446 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

     ▅██▃            ▂▃▁                                       ▂
  ▇█▆████▅▄▁▁▁▁▁▁▁▃▁▅███▄▁▁▁▁▄▇█▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▆▅▄▃▇█ █
  15.9 μs       Histogram: log(frequency) by time      19.1 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark ldiv!(UpperTriangular($B), copyto!($C, $A))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  15.085 μs …  37.057 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     15.274 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   15.373 μs ± 733.025 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

   ▄▇█▇▄           ▃▂                                          ▂
  ██████▇▃▁▃▁▁▃▁▁▅████▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▆▃▅▅▆█ █
  15.1 μs       Histogram: log(frequency) by time      18.3 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.
```
Multithreaded benchmarks:
```julia
julia> BLAS.set_num_threads(min(Threads.nthreads(), TriangularSolve.VectorizationBase.num_cores()))

julia> @benchmark TriangularSolve.rdiv!($C, $A, UpperTriangular($B))
BenchmarkTools.Trial: 10000 samples with 3 evaluations.
 Range (min … max):  8.309 μs …  24.357 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     8.769 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   8.812 μs ± 382.702 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

               ▁▃▄▆▆██▇▆▅▃▁
  ▂▁▂▂▂▂▃▃▃▄▅▇██████████████▇▆▅▄▃▃▃▃▃▂▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂ ▄
  8.31 μs         Histogram: frequency by time         9.7 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark rdiv!(copyto!($C, $A), UpperTriangular($B))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  11.996 μs … 151.147 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     14.163 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   14.281 μs ±   2.372 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

          ▂▄▇███▇▆▅▃▂ ▁   ▂▄▄▅▅▅▆▃▃         ▁
  ▁▁▁▂▂▃▄▇██████████████████████████▇▆▅▄▅▆▇███▆▅▅▃▄▂▂▂▁▁▁▁▁▁▁▁ ▅
  12 μs           Histogram: frequency by time         17.3 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark TriangularSolve.ldiv!($C, LowerTriangular($B), $A)
BenchmarkTools.Trial: 10000 samples with 5 evaluations.
 Range (min … max):  7.903 μs …  22.442 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     9.871 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   9.789 μs ± 864.957 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▂▃  ▄▃  ▃▅   ▅▃   ▆▂   ▆▄   ▂▇▄   ▃█▅▂▂▁▁▄▆▃▁ ▁             ▂
  ██▅▂██▆▅██▆▆▆██▇▇███▇▇▇████▇█████▆██████████████▇███▇▇▆▇▆▅▆ █
  7.9 μs       Histogram: log(frequency) by time      11.8 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark ldiv!($C, LowerTriangular($B), $A)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  13.507 μs … 142.574 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     15.258 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   15.319 μs ±   2.045 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

     ▁▃   ▁▂   ▁▃▅▁  ▁▄▄▁  ▂▆█▆▃
  ▁▂▅███▆▇███▆▅████▆▅████▆▆█████▆▄▄▆▆▅▄▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁ ▄
  13.5 μs         Histogram: frequency by time         18.5 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> versioninfo()
Julia Version 1.8.0-DEV.438
Commit 88a6376e99* (2021-08-28 11:03 UTC)
Platform Info:
  OS: Linux (x86_64-redhat-linux)
  CPU: 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, tigerlake)
Environment:
  JULIA_NUM_THREADS = 8
```
Single-threaded benchmarks on an M1 mac:
```julia
julia> N = 100;

julia> A = rand(N,N); B = rand(N,N); C = similar(A);

julia> @benchmark TriangularSolve.rdiv!($C, $A, UpperTriangular($B), Val(false)) # false means single threaded
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  21.416 μs …  34.458 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     21.624 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   21.767 μs ± 491.788 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

    ▃ ▆██ ▆▄ ▁                 ▃▄ ▄▂          ▁            ▂▃▁ ▂
  ▃▇█▁███▁██▁█▆▁▁▁▁▁▁▁▁▁▁▁▁▁▃█▁██▁███▁▆▃▁▁▆▇▁██▁█▆▅▁▄▃▁▃▃▇▁███ █
  21.4 μs       Histogram: log(frequency) by time      23.2 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark rdiv!(copyto!($C, $A), UpperTriangular($B))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  39.124 μs … 57.749 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     46.166 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   46.274 μs ±  1.766 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

                              ▁▁▄▂▆▃█▅▇▄▇▅▃▃▁▃▁▂               
  ▂▁▁▂▂▂▂▂▁▂▂▂▂▂▂▃▃▃▃▃▄▄▅▅▆▅▇▇████████████████████▆▇▆▆▅▆▅▅▄▃▃ ▅
  39.1 μs         Histogram: frequency by time        50.2 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark ldiv!($C, LowerTriangular($B), $A)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  48.291 μs …  57.833 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     49.124 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   49.306 μs ± 802.143 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

    ▁▃▅▆▇██▇██▇▇▆▅▄▂▂▁▁▁▂▁▁▁▁▁▁▁ ▁▁▁                           ▃
  ▃████████████████████████████████████▇▆▄▂▄▃▂▃▃▄▄▃▆▅▇▇▇██▇█▇▇ █
  48.3 μs       Histogram: log(frequency) by time        53 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark TriangularSolve.ldiv!($C, LowerTriangular($B), $A, Val(false)) # false means single threaded
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  34.249 μs …  40.208 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     34.375 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   34.748 μs ± 774.675 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▆██▆▃▄▅▃                ▁▁▄▅▅▃▂▁                     ▂▃▂  ▁▂ ▂
  ████████▁▁▃▁▁▁▁▁▃▄▃▁▁▃██████████▇▅▄▅▅▆▄▄▄▄▄▅▄▄▃▅▃▄▃▅█████▇██ █
  34.2 μs       Histogram: log(frequency) by time      37.1 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.
```
Or
```julia
julia> @benchmark TriangularSolve.ldiv!($C, LowerTriangular($B), $A, Val(false)) # false means single threaded
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  23.750 μs …  30.541 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     23.875 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   23.948 μs ± 316.293 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

   ▃▁▆ █ ▇▆▆ ▄ ▁                               ▁ ▁         ▁ ▁ ▂
  ▅███▆█▁███▄█▁██▇▁▄▁▁▁▁▁▃▁▁▁▁▁▁▁▃▁▁▁▃▁▁▁▁▁▆▁▇▆█▁█▁▇▆▅▁▅▁▇▆█▁█ █
  23.8 μs       Histogram: log(frequency) by time        25 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.
```

For editing convenience (you can copy/paste the above into a REPL and it should automatically strip `julia> `s and outputs, but the above is less convenient to edit if you want to try changing the benchmarks):
```julia
using TriangularSolve, LinearAlgebra, MKL;
BLAS.set_num_threads(Threads.nthreads())
BLAS.get_config().loaded_libs
N = 100;

A = rand(N,N); B = rand(N,N); C = similar(A);

@benchmark TriangularSolve.rdiv!($C, $A, UpperTriangular($B), Val(false))
@benchmark rdiv!(copyto!($C, $A), UpperTriangular($B))

@benchmark TriangularSolve.ldiv!($C, LowerTriangular($B), $A, Val(false))
@benchmark ldiv!($C, LowerTriangular($B), $A)

BLAS.set_num_threads(TriangularSolve.VectorizationBase.num_cores())
@benchmark TriangularSolve.rdiv!($C, $A, UpperTriangular($B))
@benchmark rdiv!(copyto!($C, $A), UpperTriangular($B))

@benchmark TriangularSolve.ldiv!($C, LowerTriangular($B), $A)
@benchmark ldiv!($C, LowerTriangular($B), $A)

versioninfo()
```

Currently, `rdiv!` with `UpperTriangular` and `ldiv!` with `LowerTriangulra` matrices are the only supported configurations.


