using Pkg
Pkg.build(; verbose = true)
Pkg.test(coverage=true)
using Test
using PoolGPAS

@test println("Testing!")
