
using Test

testlist = [
    ("nodedynamics.jl", "Single Node Tests"),
    ("dynamicnodemacro.jl", "Dynamic Node Macro Tests"),
    ("complexview.jl", "Complex View Tests"),
    ("griddynamics.jl", "Grid Construction Tests"),
    ("outputanderrors.jl", "Output and Error Tests"),
    ("states.jl", "States Tests"),
]

@testset "$desc" for (file, desc) in testlist
    @time include(file)
end
