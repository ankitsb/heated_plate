using SpecialFunctions
using Pkg
## Pkg.add(name=”Oceananigans”, version=”0.82.0”) # use it once
Pkg.activate("/g/data/hh5/tmp/ab2462/julia_work/OOC_LES")
using Oceananigans
using Printf
using NCDatasets

Nx = Ny = 128 # 64
Nz = 256  # 128
Lx = Ly = 2
Lz = 1
L0 = Lz
grid = RectilinearGrid(size = (Nx, Nx, Nz), # GPU(), # only if gpu available
                                 x = (0, Lx),
                                 y = (0, Ly),
                                 z = (0, Lz))
ν = κ = 1.5e-5
B0 = 8.5e10 * ν * κ ^2 / L0^4
b0 = 9e8 * ν * κ / L0^3
δi = κ * b0 / B0
B_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(B0),
                                       top = GradientBoundaryCondition(0))
u_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(0), bottom = GradientBoundaryCondition(0))

bottom_sponge_layer = Relaxation(; rate = 1/1,
                                   target = LinearTarget{:z}(intercept=0, gradient=0),
                                   mask = GaussianMask{:z}(center=Lz, width=Lz/6))


model = NonhydrostaticModel(; grid,
                                   advection = WENO(order=5),
                                   tracers = (:b), buoyancy = BuoyancyTracer(),
                                   coriolis = nothing,
                                   closure = AnisotropicMinimumDissipation(),
                                   boundary_conditions = (u=u_bcs, b=B_bcs),
 forcing=(b=bottom_sponge_layer, w=bottom_sponge_layer, u= bottom_sponge_layer))

#Ξ(z) = randn() * z / model.grid.Lz * ( (1 - z) / model.grid.Lz) # noise
Ξ(z) = randn() * ( 1 - tanh(π * z / model.grid.Lz / 100) ) # noise

bᵢ(x, y, z) = (1 - erf(sqrt(π)* z / 2δi)) * b0 * (1 + 1e-3 * Ξ(z))
uᵢ(x, y, z) = (1e-3 * Ξ(z))
set!(model, u=uᵢ, v=uᵢ, w=uᵢ, b=bᵢ)
simulation = Simulation(model, Δt=0.0001, stop_time=50)
wizard = TimeStepWizard(cfl=0.1, max_change=1.1, max_Δt=0.01)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e ms⁻¹, wall time: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt),
 maximum(abs, sim.model.velocities.w), prettytime(sim.run_wall_time))

simulation.callbacks[:print] = Callback(progress_message, IterationInterval(20))

u, v, w = model.velocities
b = model.tracers.b
nu = model.diffusivity_fields.νₑ
outputs = (; u, v, w, b, nu)
simulation.output_writers[:fields] =  NetCDFOutputWriter(model, outputs;
                                      filename = joinpath(@__DIR__, "heated_bottom_flat_plate.nc"),
                                      schedule = TimeInterval(10),
                                      overwrite_existing = true)

run!(simulation)
