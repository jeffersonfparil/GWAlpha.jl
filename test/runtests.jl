using Pkg
using Test
Pkg.build(; verbose = true)

Test.@test println("Testing!")
