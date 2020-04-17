# File to be loaded into eFEMpart

abstract type AbstractSolution  end
abstract type AbstractLinearOperator end

############################
###### Solution Types ######
############################

mutable struct FluidSolution <: AbstractSolution
  u::Vector{Float64}
  v::Vector{Float64}
  p::Vector{Float64}
end

mutable struct ScalarSolution <: AbstractSolution
  u::Vector{Float64}
end

#####################################
###### LinearOperator Type ##########
#####################################

mutable struct LinearOperator{T1,T2} <: AbstractLinearOperator
  Op::SparseMatrixCSC{T1,T2}
  rhs::Vector{T1}
end

##### Constructor for LinearOperator ######
function LinearOperator(mesh::AbstractMesh,prob::AbstractProblem)
  Op  = spzeros(length(mesh.xy),length(mesh.xy))
  rhs = zeros(length(mesh.xy))

  return LinearOperator(Op,rhs)
end

function LinearOperator(mesh::FluidMesh,prob::AbstractProblem)
  nu = length(mesh.xy)
  np = length(mesh.xyp)
  nTotal = 2*nu+np
  Op  = spzeros(nTotal,nTotal)
  rhs = zeros(nTotal)

  return LinearOperator(Op,rhs)
end
