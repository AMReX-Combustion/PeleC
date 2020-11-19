# EB Verification test suite

|    | Description                                                                                    | Type                                  | H | D | S | R |
|----|------------------------------------------------------------------------------------------------|---------------------------------------|---|---|---|---|
| 1  | [EB MMS](C1/README.md)                                                                         | convergence                           | x | x |   |   |
| 2  | [Quiescent flow on EB cylinder intersect with sphere](C2/README.md)                            | maintain zero flow                    | x | x | x |   |
| 3  | [Auto-ignition of Case 2 but diffuse to zero](C3/README.md)                                    | ignition delay                        | x |   |   | x |
| 4  | [Case 2 but isothermal EB at 310K, flow at 300K (outflow BC)](C4/README.md)                    | steady-state should be 310K           | x | x | x |   |
| 5  | [Quiescent flow on EB tilted plane but with case 4 BC](C5/README.md)                           | examine transient                     | x | x | x |   |
| 6  | [3 orthogonal planes, not centered on grid, species diffusion with different BC](C6/README.md) | testing homog. neu. at the species EB | x | x | x |   |
| 7  | [gamma-law shock tube at an angle with AMR](C7/README.md)                                      | analytic left/right states            | x |   |   |   |
| 8  | [multi-species shock tube at an angle with AMR](C8/README.md)                                  | literature                            | x |   | x |   |
| 9  | [acoustic pulse in circular channel](C9/README.md)                                             | analytical                            | x |   |   |   |
| 10 | [Hagen–Poiseuille flow](C10/README.md)                                                         | analytical                            | x | x |   |   |
| 11 | [Multi-species shock tube with end box](C11/README.md)                                         | literature                            | x | x | x |   |
| 12 | [Advection of smooth density](C12/README.md)                                                   | analytical                            | x | x |   |   |
