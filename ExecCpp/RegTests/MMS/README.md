# Method of Manufactured Solutions

This uses the [MASA](https://github.com/manufactured-solutions/MASA)
and auto-differention tools to implement the Method of Manufactured
Solutions into Pele.

## Building
Make sure you change the install paths to match what you want it to
be.

### MetaPhysicl
1. Get [Metaphysicl](https://github.com/roystgnr/MetaPhysicL)
2. `./bootstrap` (on Eagle, do  `spack load autoconf` first)
3. `./configure --prefix=${HOME}/combustion/install/MetaPhysicL`
4. `make`
5. `make install`

### MASA

1. Get [MASA](https://github.com/manufactured-solutions/MASA)
2. `./bootstrap` (on Eagle, do  `spack load autoconf` first)
3. `./configure --enable-fortran-interfaces METAPHYSICL_DIR=${HOME}/combustion/install/MetaPhysicL --prefix=${HOME}/combustion/install/MASA --enable-python-interfaces`
4. `make`
5. `make check`
6. `make install`

## Using

MASA must be linked to when building Pele. We will be using
`ad_cns_3d_les` to manufacture a source and initial condition for
Pele. The user must have defined the `MASA_HOME` variable to point to
the install location of MASA.
