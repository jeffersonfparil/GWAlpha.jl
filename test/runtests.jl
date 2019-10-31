using Pkg
Pkg.build(; verbose = true)

@test println("Testing!")
