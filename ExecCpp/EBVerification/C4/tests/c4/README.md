# STEADY STATE THERMAL DIFFUSION TEST

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

