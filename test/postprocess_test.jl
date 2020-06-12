using eFEM, Test

mesh1 = squareMesh([-2,2,-1,1],12,1)
OperatorType = :Poisson2D

dNodes = Dirichlet(:left,:right,:bottom,:top)
dBCf = Dirichlet((x,y) -> x^4*y^4)
ff   = Forcing((x,y) -> -12*x^2*y^2*(x^2+y^2))

Nodes = [dNodes]; bcfun = [dBCf,ff]
prob = Problem(mesh1,Nodes,bcfun,OperatorType)
sol = solve(prob,mesh1)

L1   = DomainNorm(mesh1.xy,mesh1.cm,sol.u;normID="1")
L2   = DomainNorm(mesh1.xy,mesh1.cm,sol.u;normID="2")
Linf = DomainNorm(mesh1.xy,mesh1.cm,sol.u;normID="Inf")
hCalc(mesh1)

mesh2 = squareMesh([-2,2,-1,1],12,2)
OperatorType = :Poisson2D

dNodes = Dirichlet(:left,:right,:bottom,:top)
dBCf = Dirichlet((x,y) -> x^4*y^4)
ff   = Forcing((x,y) -> -12*x^2*y^2*(x^2+y^2))

Nodes = [dNodes]; bcfun = [dBCf,ff]
prob = Problem(mesh2,Nodes,bcfun,OperatorType)
sol = solve(prob,mesh2)

L1   = DomainNorm(mesh2.xy,mesh2.cm,sol.u;normID="1")
L2   = DomainNorm(mesh2.xy,mesh2.cm,sol.u;normID="2")
Linf = DomainNorm(mesh2.xy,mesh2.cm,sol.u;normID="Inf")
hCalc(mesh2)
