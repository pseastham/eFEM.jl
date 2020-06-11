using eFEM, Test

# mesh input from gmsh
mesh1 = Mesh("test-data/step.msh")
mesh2 = FluidMesh("test-data/step.msh")
