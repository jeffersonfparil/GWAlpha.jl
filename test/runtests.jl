using Pkg
using Test
Pkg.build(; verbose = true)

@test println("Testing!")
