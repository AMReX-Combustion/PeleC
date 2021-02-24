# EB-C4: STEADY STATE THERMAL DIFFUSION TEST

## Case description

|                    | description                                         |
|:-------------------|:----------------------------------------------------|
| Problem definition | Heat conduction from body into quiescent environ    |
| EB geometry        | Same as C2: Cylinder intersecting Sphere            |
| Mechanism          | LiDryer                                             |
| Multi-level        | no                                                  |
| Metric             | steady-state temp matches body surface temp         |

## Notes

BC on body surface : 301 K
BC on domain bundaries : outflow or insulated walls
IC in domain : 300 K

No species transport so LiDryer mechanism is not needed, but is used for 
compatibility with C2 and C3.

## Results


# EB-C5: Stokes 1st Problem / Thermal Diffusion into a semi-infinite medium

## Case description

|                    | description                                         |
|:-------------------|:----------------------------------------------------|
| Problem definition | diffusive processes into semi-infinite medium       |
| EB geometry        | planar surface at arbitrary angle in 3D             |
| Mechanism          | none                                                |
| Multi-level        | no                                                  |
| Metric             | temporally evolving profiles match the self-similar |
|                    |   analytical solution                               |

## Notes

Identical problem setup can be used to validate thermal diffusion and
viscous diffusion (impulsively started plate). Boundary conditions would
be more of an issue for the viscous diffusion problem.

Due to bounary effects, only validate profiles in the normal direction
near the center of the planar surface where boundary effect will be smallest
for either case.

BC on planar surface : 301 K stationary or 300K w/in-plane velocity
BC on domain bundaries : insulated walls or outflow
IC in domain : 300 K, quiescent

## Results
