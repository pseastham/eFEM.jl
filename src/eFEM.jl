module eFEM

# Custom types
include("MeshTypes.jl")
include("ParameterTypes.jl")
include("ProblemTypes.jl")
include("SolutionTypes.jl")

# Mesh input/generation
include("MeshGeneration.jl")
include("GMSHreader.jl")
include("MeshTransform.jl")

# Generate Matrices
include("Basis.jl")
include("BilinearForms.jl")
include("MatrixGeneration.jl")
include("BoundaryConditions.jl")

# Solve Problem
include("Assembly.jl")
include("TimeStepping.jl")
include("Solver.jl")

# PostProcessing
include("vtkExport.jl")
include("ErrorAnalysis.jl")
include("DomainQuadrature.jl")
include("PostProcessing.jl")
include("TracerGenerate.jl")

export AbstractMesh

export getBoundaries

# ParameterTypes.jl
export PoissonParam,
       HeatParam,
       AdvDiffParam,
       DarcyParam,
       FluidParam,
       BrinkmanParam,
       BrinkmanMPParam,
       FluidVarParam,
       AbstractVariableParameter,
       AbstractConstantParameter

# ProblemTypes.jl
export Problem,
       Dirichlet,
       Neumann,
       Robin,
       Forcing

# SolutionsTypes.jl
export ScalarSolution,
       FluidSolution

# MeshGeneration.jl
export squareMesh,
       squareMeshFluid,
       Mesh,
       FluidMesh,
       axisymMesh,
       axisymFluid

# MeshTransform.jl
export pointTransform,
       onSegment,
       doIntersect,
       isInsideDomain,
       meshTransform

# Basis.jl
export derivShape2D,
       shapeEval,
       GaussEdge,
       FEMstats

# MatrixGeneration.jl
export Mass2DMatrix

# TimeStepping.jl
export FEMForwardEuler,
       progressBar,
       sToHMS

# Solver.jl
export solve,
       solve!,
       GenerateSystem,
       ApplyBC!

# vtkExport.jl
export vtksave,
       Path,
       ScalarData,
       ScalarNames,
       VectorData,
       VectorNames,
       SURF_UNSTRUCTURED_TO_VTK

# ErrorAnalysis.jl
export DomainNorm,
       hCalc

# DomainQuadrature.jl
export DomainQuad,
       SurfaceQuad,
       SurfaceFlux,
       mySurface,
       VelocityFlux,
       SurfaceInterp

# PostProcessing.jl
export computeStress,
       computeDerivative,
       computeExtension,
       computeCompressibility

# TracerGenerate.jl
export GenerateTracers,
       Tracer,
       TracerInfo,
       TracerLineSource

end # module
