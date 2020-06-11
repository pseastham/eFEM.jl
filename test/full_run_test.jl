using eFEM, Test

@testset "Laplace" begin
  # define mesh through eFEMpart code
  mesh = squareMesh([-1.2,1.1,-1.2,1.1],12,1)

  # define variable name and operator type
  OperatorType = :Poisson2D

  # Define auxiliary information
  dNodes = Dirichlet(:left,:right)
  nNodes = Neumann(:top,:bottom)

  # Define dirichlet, neumann, and forcing functions
  dBCf = Dirichlet((x,y) -> (x==-1.2 ? -1.0 : y))
  nBCf = Neumann((x,y) -> 0.0)
  ff   = Forcing((x,y) -> 0.0)

  # required technicality for grouping auxiliary information
  Nodes = [dNodes,nNodes]; bcfun = [dBCf,nBCf,ff]

  # define problem variable
  prob = Problem(mesh,Nodes,bcfun,OperatorType)

  # solve problem
  sol = solve(prob,mesh)

  # save solution
  vtkname = Path("test-data","laplace_soln")
  sd = ScalarData(sol.u)
  sn = ScalarNames("testing_var")
  vtksave(mesh,sd,sn,vtkname)

  rm("test-data/laplace_soln.vtk")
end