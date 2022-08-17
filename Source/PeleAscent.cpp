#include <AMReX_ParmParse.H>
#include "PeleAscent.H"

namespace pele {
PeleAscent::PeleAscent()
{
  {
    amrex::ParmParse pp("ascent");
    pp.query("plot_int", plot_int);
  }
}
} // namespace pele
