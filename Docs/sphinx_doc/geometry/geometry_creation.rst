.. _geometry-creation:



GeometryShop and Implicit Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


One of the greatest advantages of EB technology is that grid generation is robust and fast and can be done to any accuracy as described by `Schwartz et al. <http://dx.doi.org/10.2140/camcos.2015.10.83>`_ The foundation class that AMReX uses for
geometry generation is called :ref:`GeometryShop <geom>`. This class is used to initialize geometric information and associated connectivity graph stored in a distributed database class `EBIndexSpace`. Historically, the `EBIndexSpace` database was developed to be used throughout a calculation. Here, we use it only to populate datastructures that can be accessed efficiently in the patterns representative of the Pele motivating problem space. 

 Given an implicit function :math:`I`, ``GeometryShop`` interprets the surface upon which :math:`I(\mathbf{x}) = 0` as the surface with which to cut the grid cells. ``GeometryShop`` interprets the positive regions of the implicit function (:math:`\mathbf{x}: I(\mathbf{x}) > 0`) as covered by the geometry and negative regions (:math:`\mathbf{x}: I(\mathbf{x}) < 0`) as part of  the solution domain.  For example, if one defines her implicit function :math:`S` as

.. math:

   S(\mathbf{x}) = x^2 + y^2 + z^2 - R^2,

the solution domain would be the interior of a sphere of radius :math:`R`. Reverse the sign of :math:`S` and the solution domain would be the exterior of the sphere.   

To define a geometry shop, one needs to send it a BaseIF, which describes an implicit function. 

.. code-block:: c

    GeometryShop(const BaseIF& a_localGeom)

The implicit function's interface needs to be able to do two things: create a copy of itself and return the value of the function at any point in space.

.. code-block:: c

    /// Return the value of the function at a_point.  
    virtual Real value(const RealVect& a_point) const = 0;

    ///   Return a newly allocated derived class.  
    virtual BaseIF* newImplicitFunction() const = 0;


To continue with the example above, if one wants to define a geometry
as a domain with  a sphere cut out of it, one uses the ``SphereIF`` class, the functions of which are:

.. code-block:: c

    SphereIF::
    SphereIF(const Real&     a_radius,
             const RealVect& a_center,
             const bool&     a_inside)
    {
     m_radius  = a_radius;
     m_radius2 = m_radius*m_radius;
     m_inside  = a_inside;
     m_center  = a_center;
    }

    Real
    SphereIF::
    value(const RealVect& a_point) const
    {
      RealVect dist = a_point - m_center;
      Real distance2 = dist.radSquared();
      Real retval = distance2 - m_radius2;
      // Change the sign to change inside to outside
      if (!m_inside)
        {
          retval = -retval;
        }
  
      return retval;
    }
    BaseIF* 
    SphereIF::
    newImplicitFunction() const
    {
      SphereIF* spherePtr = new SphereIF(m_radius,
                                         m_center,
                                         m_inside);
  
      return static_cast<BaseIF*>(spherePtr);
  }

One can get by with very simple implicit functions such as sphere and plane because one can create more complicated geometries by composition of simple implicit functions. AMReX contains the following classes which are used compose implicit functions:

* ``TransformIF`` allows for translations and rotations of an implicit function.
* ``UnionIF`` produces the union of two implicit functions.  
* ``IntersectionIF`` produces the intersection of two implicit functions.
* ``LatheIF`` creates a 3D implicit function as the surface of revolution of a 2D implicit function.

Here is an example that uses many of these tools.  This example creates a geometry with multiple spheres cut out.

.. code-block:: c

   //say you have a bunch of radii and centers of spheres
   /* fill these in however you like */
   vector<Real>     radius(numSpheres);
   vector<RealVect> center(numSpheres);
   ...
   //create an implicit function for each sphere
   vector<BaseIF*>  spheres(numSpheres);
   
   for(int isphere = 0; isphere < numSpheres; isphere++)
   {
     //create each sphere at the origin and translate
     SphereIF sphereAtZero(radius[isphere], RealVect::Zero, false);
     TransformIF* movedSphere = new TransformIF(sphereAtZero);
     movedSphere->translate(center[isphere]);
     spheres[isphere] = static_cast<BaseIF*>(movedSphere);
   }
   //create an implicit function that is the intersection of all the spheres
   IntersectionIF impMultisphere(spheres);
   //we want the fluid to be the complement (the space outside the sphere
   ComplementIF sideImpMultisphere(impMultisphere, false);
   //create the geometryshop
   GeometryShop workshop(sideImpMultisphere)

Specifics of EBDemo3D Geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To make a geometry similar to the combustor geometry in EBDemo3D, we first read in from the input file the points of the polygons that form the two-dimensional (2D) description of the geometry (see :ref:`polygons`).   We create each polygon as a union of line implicit functions.   The union of the polygons forms the 2D geometric description.   For the 3D geometric description, we form the two dimensional description and then use the surface of revolution of the 2D about the $z$ axis as our 3D surface.  :ref:`surface` shows the 3D surface and :ref:`slice` shows slices of the 3D domain along coordinate planes. 

The geometry, included in the inputs file is:

.. code::

   geom_type = "polygon_revolution"
   num_poly = 3
   #far wall (all these polys are done as fractions of the domain)
   poly_0_num_pts  = 4
   poly_0_point_0  = 0.45  -1.0 0.0
   poly_0_point_1  = 2.0    -1.0 0.0
   poly_0_point_2  = 2.0     2.0 0.0
   poly_0_point_3  = 0.45   2.0 0.0
   
   #ramping bit
   poly_1_num_pts  = 5
   poly_1_point_0  = 0.1    -1.0 0.0
   poly_1_point_1  = 2.0    -1.0 0.0
   poly_1_point_2  = 2.0     0.6 0.0
   poly_1_point_3  = 0.25    0.6 0.0
   poly_1_point_4  = 0.1     0.1 0.0
   
   #pipe
   poly_2_num_pts  = 4
   poly_2_point_0  = 0.06 -1.0 0.0
   poly_2_point_1  = 0.08    -1.0 0.0
   poly_2_point_2  = 0.08    0.35 0.0
   poly_2_point_3  = 0.06   0.35 0.0

Here is how the code to generate this geometry in 3D looks.

.. code-block:: c


    // Data for polygons making up nozzle
    vector<vector<RealVect> > polygons;
    ///fill the polygons point information from the input file
    ...
    
    // For building each polygon
    
    int num_poly;
    RealVect translation;
    //this is the amount to translate to get the center of the domain
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      translation[idir] = 0.5*n_cell[idir]*fine_dx;
    }
    
    // Make the vector of (convex) polygons (vectors of points) into a union
    // of convex polygons, each made from the intersection of a set of half
    // planes/spaces - all represented by implicit functions.
    
    UnionIF* crossSection = makeCrossSection(polygons);
    
    // In 3D rotate about the z-axis and complement if necessary
    LatheIF lathe(*crossSection,false);
    //we are starting around the z axis so we need to translate
    //over to the center 
    
    translation[2] = 0;
    TransformIF implicit(lathe);
    implicit.translate(translation);
    impfunc = implicit.newImplicitFunction();
    
    //create the geometryshop object
    GeometryShop geom(*impfunc);

.. _polygons:

.. figure:: polygons.pdf
   :alt: Polygons
   :width: 500

   Polygons

   Polygons which form 2D combustion geometry.  The three polygons are cut out of the solution domain.   For a 3D geometry, we use the surface of revolution of these polygons to cut the geometry then translate that surface to the center of the domain.

.. _slice:

.. figure:: ./slice3d.pdf
   :alt: 3dslice
   :width: 500

   Slices

   Slices along coordinate planes of 3D combustor geometry.


.. _surface:

.. figure:: surface3d.pdf
   :alt: 3dsurface
   :width: 500

   Cutting surface



