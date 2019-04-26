# Laminar channel flow with embedded boundaries

This case is meant to be used to study the pressure driven laminar
channel flow (i.e. laminar viscous flow between two parallel plates)
where the channel is rotated at various angles with respect to the
mesh.

The `gen_channel_input.py` script generates the input lines for the
embedded boundary definition and the `probin` file for various
cases. The typical way to invoke the python script is

```console
$ ./PeleC3d.gnu.MPI.ex inputs.3d "$(python gen_channel_input.py -a 45 -r 100 -n 8)"
```

This runs a Re=100 channel oriented at a 45 degree with a resolution
of 8 cells per channel half width. The default Mach number is 0.1 and
the default Prandtl number is 0.71. Special care is taken at the
inflow and outflow to apply the pressure boundary condition to
correctly represent the pressure driven channel. The initial flow
field is a 10% perturbation of the exact incompressible solution for
this flow.
