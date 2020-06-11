using eFEM, Test

@testset "Laplace order 1" begin
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

@testset "Advection-Diffusion order 2" begin
  mesh = squareMesh([-1.2,1.1,-1.2,1.1],12,2)

  OperatorType = :AdvDiff2D
  
  dNodes = Dirichlet(:left,:top)
  nNodes = Neumann(:right,:bottom)

  dBCf   = Dirichlet((x,y) -> (x==-1.2 ? -0.3 : 0.0))
  ff     = Forcing((x,y) -> 0.0)
  
  Pe = 7.3
  windX(x,y) =  -y;  windY(x,y) = x
  wx = [windX.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
  wy = [windY.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
  param  = AdvDiffParam(wx,wy,Pe)

  Nodes = [dNodes]
  bcfun = [dBCf,ff]
  
  prob = Problem(mesh,Nodes,bcfun,OperatorType)
  sol = solve(prob,mesh,param)

  vtkname = Path("test-data","advdiff_soln")
  sd = ScalarData(sol.u)
  sn = ScalarNames("testing_var")
  vd = VectorData([wx,wy])
  vn = VectorNames("velocity")
  vtksave(mesh,sd,sn,vd,vn,vtkname)

  rm("test-data/advdiff_soln.vtk")
end