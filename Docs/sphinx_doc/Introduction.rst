.. highlight:: rst

.. Warning:: This documentation is a placeholder, and contents should currently be considered a work in progress, out of context, or just plain wrong until this note is removed!


Introduction
============

Pele solves the reacting compressible Navier-Stokes on a structured grid with, optionally, embedded boundary geometry treatment and non-ideal gas equations of state. A variety of time advance schemes are implemented, notably an operator-split (Strang) approach and an SDC based approach. A variety of examples are included to provide a model setup for the various options. These are discussed futher in the :ref:`Getting Started<GettingStarted>` section.

Dependencies
------------

Pele is built on AMReX (available at `https://github.com/AMReX-Codes/amrex <https://github.com/AMReX-Codes/amrex>`_), an adaptive mesh refinement software framework, which provides the underlying software infrastructure for block structured AMR operations. The full AMReX documentation can be found `here <https://amrex-codes.github.io/AMReXUsersGuide.pdf>`_. 

Modules for describing the equation of state, diffusion transport model, and reaction kinetics are localized to the ``PelePhysics`` repository. For the purpose of this Users' Guide  ``PelePhysics`` is considered part of ``PeleC`` but it needs to be obtained through a separate checkout as described in the :ref:`readme <README>`.


Development
-----------

A separate developers guide does not yet exist; along with the algorithmic description in this Users' Guide, doxygen documentation exists in place and an input file exists in `PeleC/Docs` that can be build using:

::

	doxygen Doxyfile
