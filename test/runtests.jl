using Pkg
using Test
Pkg.build() # Pkg.build(; verbose = true) for Julia 1.1 and up
Pkg.test(coverage=true)

@test println("Testing!")
