using ModiaMath
using PyCall
#https://github.com/ModiaSim/ModiaMath.jl/wiki/Installing-PyPlot-in-a-robust-way
pygui(:gtk)
using PyPlot
using Modia
  @model FirstOrder begin
     x = Variable(start=1)   # start means x(0)
     T = Parameter(0.5)      # Time constant
     u = 2.0                 # Same as Parameter(2.0)
  @equations begin
     T*der(x) + x = u        # der() means time derivative
     end
  end;

result = simulate(FirstOrder, 2);
@show result["x"][end];
ModiaMath.plot(result, "x")
