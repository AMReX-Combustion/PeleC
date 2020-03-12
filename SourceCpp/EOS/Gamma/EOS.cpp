#include "EOS.H"

namespace EOS {

AMREX_GPU_DEVICE_MANAGED int upassMap[NPASSIVE];
AMREX_GPU_DEVICE_MANAGED int qpassMap[NPASSIVE];
AMREX_GPU_DEVICE_MANAGED amrex::Real gamma = 1.4;

void
init()
{
  amrex::ParmParse pp("eos");
  pp.query("gamma", gamma);
  AMREX_D_PICK(upassMap[0] = UMY;
	       qpassMap[0] = QV;
	       upassMap[1] = UMZ;
	       qpassMap[1] = QW;
	       int curMapIndx = 2;,
	       upassMap[0] = UMZ;
	       qpassMap[0] = QW;
	       int curMapIndx = 1;,
	       int curMapIndx = 0;);
  for (int i = 0; i != NUM_ADV; ++i) {
    upassMap[curMapIndx] = i + UFA;
    qpassMap[curMapIndx] = i + QFA;
    curMapIndx++;
  }
  for (int i = 0; i != NUM_SPECIES; ++i) {
    upassMap[curMapIndx] = i + UFS;
    qpassMap[curMapIndx] = i + QFS;
    curMapIndx++;
  }
  AMREX_ASSERT(curMapIndx == NPASSIVE);
}

} // namespace EOS

void
pc_eos_init()
{
  EOS::init();
}


