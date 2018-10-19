.. _amrex-geometry:



There are a few key AMReX facilities used to enable the EB capability.

.. function:: AMREX_USE_EB

    Preprocessor variable to enable AMReX-EB related operations. In PeleC used to control call of initialization routines for EB geometry (in main.cpp).


.. doxygenclass:: amrex::EBFArrayBox
   :members:
   :undoc-members:

.. doxygenfunction:: amrex::EBFArrayBox::getEBCellFlagFab

.. doxygenenum:: amrex::FabType

.. _geom:

.. doxygenclass:: amrex::GeometryShop


.. f:function:: amrex_ebcellflag_module/get_neighbor_cells
    :p flag: Flag encoding information about neighbor cells
    :p nbr: A 3x3x3 (3x3 in 2D) returning 0 if cell is covered

    This returns 1 if the neighbor cell is connected, and 0 if not. There is an alternative version `get_neighbor_cells_real` that returns floating point values.


.. f:function:: amrex_ebcellflag_module/is_regular_cell
    :p flag: Flag encoding information about neighbor cells

    This returns 1 if the cell is a regular cell, and 0 if it is a cut cell. 

