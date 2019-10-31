using Pkg
using Base.Test
using GWAlpha
Pkg.build(; verbose = true)

@test println("Testing!")
