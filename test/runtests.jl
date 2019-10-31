using Pkg
Pkg.build(; verbose = true)
Pkg.test(coverage=true)
using PoolGPAS

@test println("Testing!")
