using Pkg
Pkg.build(; verbose = true)
using PoolGPAS

@test println("Testing!")
