#include "prob.H"

namespace ProbParm {
} // namespace ProbParm

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{
}
}

#ifdef DO_PROBLEM_POST_TIMESTEP
void
PeleC::problem_post_timestep()
{
}
#endif

#ifdef DO_PROBLEM_POST_INIT
void
PeleC::problem_post_init()
{
}
#endif
