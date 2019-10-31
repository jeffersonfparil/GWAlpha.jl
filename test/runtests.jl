using Pkg
Pkg.build() # Pkg.build(; verbose = true) for Julia 1.1 and up
Pkg.test(coverage=true)
using Test
using PoolGPAS

@test println("Testing!")
