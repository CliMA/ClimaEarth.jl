#=
# Various things
julia> atmos.integrator.p.precomputed.sfc_conditions.
buoyancy_flux   obukhov_length
ts              ustar
ρ_flux_h_tot    ρ_flux_q_tot

# Thermodynamic state:
julia> atmos.integrator.p.precomputed.sfc_conditions.ts
Thermodynamics.PhaseNonEquil{Float32}-valued Field:
  e_int: Float32[-52994.0, -54332.2, -56341.9, -57425.2, -49488.4, -50800.1, -52759.9, -53877.0, -43064.2, -44231.1  …  -61447.5, -62595.5, -54915.9, -56354.6, -58530.2, -59695.0, -52994.0, -54332.2, -56341.9, -57425.2]
  ρ: Float32[1.12521, 1.11691, 1.1047, 1.09776, 1.14439, 1.13845, 1.1275, 1.12067, 1.17122, 1.16716  …  1.07314, 1.06576, 1.11283, 1.10381, 1.09103, 1.08377, 1.12521, 1.11691, 1.1047, 1.09776]
  q:
    tot: Float32[0.00775389, 0.00743542, 0.00696459, 0.00671588, 0.00861411, 0.00828583, 0.00780791, 0.00754089, 0.010279, 0.00996838  …  0.00581741, 0.00557031, 0.00729899, 0.00696366, 0.00646453, 0.00620355, 0.00775389, 0.00743542, 0.00696459, 0.00671588]
    liq: Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ice: Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Horizontal velocity
julia> atmos.integrator.u.c.uₕ
ClimaCore.Geometry.Covariant12Vector{Float32}-valued Field:
  components:
    data:
      1: Float32[64069.2, 61040.6, 59449.0, 56331.9, 51847.4, 46551.9, 41137.2, 36922.8, 34819.9, 35484.4  …  94650.0, 95722.5, 84585.7, 54860.3, 24166.5, 37060.7, 94621.0, 1.0586f5, -7395.97, -280682.0]
      2: Float32[-62789.3, -52728.5, -52538.0, -49500.9, -44785.2, -40040.0, -36200.9, -34041.4, -33711.5, -34926.7  …  37862.9, 56130.8, 81812.8, 87400.3, 32924.8, -83292.4, -161303.0, -57706.3, -51483.9, -118627.0]

# "int" for "interior"
# level 1: lowest level
julia> uₕ_int = CC.Geometry.UVVector.(CC.Spaces.level(sim.integrator.u.c.uₕ, 1))
# zonal component: uₕ_int.components.data.:1
# meridional component: uₕ_int.components.data.:2
=#
import ClimaAtmos as CA
import ClimaCore as CC
using ClimaParams
using Thermodynamics

# For regridding:
using ClimaUtilities
using Interpolations

AtmosInterpolations = Base.get_extension(ClimaUtilities, :ClimaUtilitiesClimaCoreInterpolationsExt)
atmos = CA.get_simulation(CA.AtmosConfig("simple_atmos_simulation.yml"))

#=
time_step!(atmos::CA.AtmosSimulation) = CA.SciMLBase.step!(atmos.integrator)

# time step
Δt = atmos.integrator.dt

# interpolation
interpolated_array = interpolate(field, hcoords, zcoords)

# λ = λnodes(grid, Center(), Center(), Center())
# φ = φnodes(grid, Center(), Center(), Center())
hcoords = @. Geometry.LatLongPoint(φ, λ)

uₕ_int = CC.Geometry.UVVector.(CC.Spaces.level(atmos.integrator.u.c.uₕ, 1))
u = CC.Remapping.interpolate(uₕ_int.components.data.:1, hcoords, nothing)
v = CC.Remapping.interpolate(uₕ_int.components.data.:2, hcoords, nothing)
# T, q, p from thermodynamic state ts?
ts_a = CC.Spaces.level(atmos.integrator.p.precomputed.ᶜts, 1)
ℂₐ = atmos.integrator.p.params.thermodynamics_params
# ℂₐ = CA.Parameters.thermodynamics_params(atmos.integrator.p.params)
T_a = Thermodynamics.air_temperature.(ℂₐ, ts_a)
p_a = Thermodynamics.air_pressure.(ℂₐ, ts)
q_a = Thermodynamics.total_specific_humidity.(ℂₐ, ts)

T = CC.Remapping.interpolate(ts_a, hcoords, nothing)

# Regridding back to atmos
space3 = axes(atmos.integrator.p.precomputed.sfc_conditions.ts)
space2 = ClimaCore.Spaces.SpectralElementSpace2D(space3.grid.full_grid.horizontal_grid)
regridder = ClimaUtilities.Regridders.InterpolationsRegridder(space2)

Q = Field{Center, Center, Nothing}(grid)

# Interpolate to intermediate Oceananigans grid
intermediate_grid = LatitudeLongitudeGrid(Float32,
                                          size=(360, 180),
                                          longitude=(0, 360),
                                          latitude=(-90, 90),
                                          topology=(Periodic, Bounded, Flat))

Qi = Field{Center, Center, Nothing}(intermediate_grid)
λi = λnodes(intermediate_grid, Center(), Center(), Center())
φi = φnodes(intermediate_grid, Center(), Center(), Center())
Oceananigans.Fields.interpolate!(Qi, Q)

# Regrid to ClimaCore grid
using ClimaCore.Utilities: half
Qa = ClimaUtilities.Regridders.regrid(regridder, interior(Qint, :, :, 1), (λi, φi))

# Project onto a vector!
# :eyes https://github.com/CliMA/ClimaEarth.jl/pull/5/files
c = atmos.integrator.p.scratch.ᶠtemp_scalar
𝒢 = ClimaCore.Fields.level(ClimaCore.Fields.local_geometry_field(c), half)
ρwh = atmos.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot
@. ρwh = CA.SurfaceConditions.vector_from_component(Qa, 𝒢)
@. ρτ = CA.SurfaceConditions.tensor_from_components(ρτxz, ρτyz, 𝒢)

#=
dict = z_max: 60000.0
z_elem: 31
dz_bottom: 50.0
dt: "400secs"
t_end: "1days"
topography: "Earth"
rayleigh_sponge: true
implicit_diffusion: true
approximate_linear_solve_iters: 2
moist: "equil"
surface_setup: "DefaultMoninObukhov"
vert_diff: "DecayWithHeightDiffusion"
precip_model: "0M"
cloud_model: "grid_scale"
toml: [sphere_aquaplanet.toml]
=#
=#
