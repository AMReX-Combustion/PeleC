#include "IndexDefines.H"

void
init_pass_map(std::unique_ptr<PassMap>& pmap)
{
  AMREX_D_PICK(
    pmap->upassMap[0] = UMY; pmap->qpassMap[0] = QV; pmap->upassMap[1] = UMZ;
    pmap->qpassMap[1] = QW; int curMapIndx = 2;, pmap->upassMap[0] = UMZ;
    pmap->qpassMap[0] = QW; int curMapIndx = 1;, int curMapIndx = 0;);
  for (int i = 0; i < NUM_ADV; ++i) {
    pmap->upassMap[curMapIndx] = i + UFA;
    pmap->qpassMap[curMapIndx] = i + QFA;
    curMapIndx++;
  }
  for (int i = 0; i < NUM_SPECIES; ++i) {
    pmap->upassMap[curMapIndx] = i + UFS;
    pmap->qpassMap[curMapIndx] = i + QFS;
    curMapIndx++;
  }
  for (int i = 0; i < NUM_AUX; ++i) {
    pmap->upassMap[curMapIndx] = i + UFX;
    pmap->qpassMap[curMapIndx] = i + QFX;
    curMapIndx++;
  }
}
