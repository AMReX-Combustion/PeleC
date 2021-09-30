#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
#include "Factory.H"
#include "EB.H"
#include "Utilities.H"

namespace pele {
namespace pelec {

class Geometry : public physics::Factory<Geometry>
{
public:
  static const std::string base_identifier() { return "Geometry"; }
  virtual void
  build(const amrex::Geometry& geom, const int max_coarsening_level) = 0;
  virtual ~Geometry() = default;
};

class FlatPlate : public Geometry::Register<FlatPlate>
{
public:
  static std::string identifier() { return "flat_plate"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

class Ramp : public Geometry::Register<Ramp>
{
public:
  static std::string identifier() { return "ramp"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

class Combustor : public Geometry::Register<Combustor>
{
public:
  static std::string identifier() { return "combustor"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

class ICEPistonBowl : public Geometry::Register<ICEPistonBowl>
{
public:
  static std::string identifier() { return "ICE_PistonBowl"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

class ExtrudedTriangles : public Geometry::Register<ExtrudedTriangles>
{
public:
  static std::string identifier() { return "extruded_triangles"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

class PolygonRevolution : public Geometry::Register<PolygonRevolution>
{
public:
  static std::string identifier() { return "polygon_revolution"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

class MovingPlane : public Geometry::Register<MovingPlane>
{
public:
  static std::string identifier() { return "moving_plane"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

class QuarterCircle : public Geometry::Register<QuarterCircle>
{
public:
  static std::string identifier() { return "quarter-circle"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

} // namespace pelec
} // namespace pele
#endif