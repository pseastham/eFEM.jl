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

  FEMstats(mesh)
  FEMstats(mesh,prob)

  # save solution
  vtkname = Path("test-data","laplace_soln")
  sd = ScalarData(sol.u,sol.u)
  sn = ScalarNames("testing_var","testing_var_2")
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
  vd = VectorData([wx,wy],[wx,wy])
  vn = VectorNames("velocity","velocity_2")
  vtksave(mesh,sd,sn,vd,vn,vtkname)

  rm("test-data/advdiff_soln.vtk")
end

@testset "Stokes Example" begin
  mesh = squareMeshFluid([-1.2,1.1,-1.2,1.1],12)

  param = FluidParam(2.3)
  OperatorType = :Stokes2D

  dNodes = Dirichlet(:left,:bottom,:top)
  nNodes = Neumann(:right)
  Nodes = [dNodes,nNodes]

  dUBCf = Dirichlet((x,y) -> (x==-1.2 ? -(1.1-y)*(-1.2-y) : 0.0))
  dVBCf = Dirichlet((x,y) -> 0.0)
  bcfun = [dUBCf,dVBCf]

  prob = Problem(mesh,Nodes,bcfun,OperatorType)
  sol = solve(prob,mesh,param)

  # save solution
  vtkname = Path("test-data","stokes_soln")
  sd = ScalarData(sol.p)
  sn = ScalarNames("pressure")
  vd = VectorData([sol.u,sol.v])
  vn = VectorNames("velocity")
  vtksave(mesh,sd,sn,vd,vn,vtkname)

  rm("test-data/stokes_soln.vtk")
end

@testset "Stokes Axisymmetric Example" begin
  mesh = squareMeshFluid([-1.2,1.1,0.0,1.1],12)

  param = FluidParam(2.3)
  OperatorType = :StokesAS

  dUNodes = Dirichlet(:left,:top,:right)
  dVNodes = Dirichlet(:left,:top,:bottom,:right)
  Nodes = [dUNodes,dVNodes]

  dUBCf = Dirichlet((x,y) -> (x==-1.2 ? -(1.1-y)*(-1.2-y) : 0.0))
  dVBCf = Dirichlet((x,y) -> 0.0)
  bcfun = [dUBCf,dVBCf]

  prob = Problem(mesh,Nodes,bcfun,OperatorType)
  sol = solve(prob,mesh,param)

  # save solution
  vtkname = Path("test-data","stokes_as_soln")
  sd = ScalarData(sol.p)
  sn = ScalarNames("pressure")
  vd = VectorData([sol.u,sol.v])
  vn = VectorNames("velocity")
  vtksave(mesh,sd,sn,vd,vn,vtkname)

  rm("test-data/stokes_as_soln.vtk")
end

@testset "Multiphase Brinkman" begin
  x0 = 2.0
  mesh = squareMeshFluid([-x0,x0,-1,1],32)

  μ0 = 1.0
  a10 = 0.0
  a20 = 1.0

  μ  = μ0*ones(length(mesh.xy))
  a1  = a10*ones(length(mesh.xy))
  a2  = a20*ones(length(mesh.xy)) 
  G  = 1e0
  param = BrinkmanMPParam(μ,a1,a2)

  OperatorType = :BrinkmanMP2D

  dNodes = Dirichlet(:left,:top,:bottom)
  nNodes = Neumann(:right)
  Nodes = [dNodes,nNodes]

  dUBCf = Dirichlet((x,y) -> (x==-x0 ? 1.0-y^2 : 0.0))
  dVBCf = Dirichlet((x,y) -> 0.0)
  bcfun = [dUBCf,dVBCf]

  prob = Problem(mesh,Nodes,bcfun,OperatorType)
  sol = solve(prob,mesh,param)

  vtkname = Path("test-data","mpb_soln")
  vd = VectorData([sol.u,sol.v])
  vn = VectorNames("velocity")
  vtksave(mesh,vd,vn,vtkname)

  rm("test-data/mpb_soln.vtk")
end