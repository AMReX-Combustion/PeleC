set(PELE_PHYSICS_EOS_MODEL Fuego)
set(PELE_PHYSICS_CHEMISTRY_MODEL SootReaction)
set(PELE_PHYSICS_TRANSPORT_MODEL Simple)
set(PELE_PHYSICS_ENABLE_SOOT ON)
set(PELE_PHYSICS_NUM_SOOT_MOMENTS 3)
set(PELE_PHYSICS_ENABLE_SPRAY OFF)
set(PELE_PHYSICS_SPRAY_FUEL_NUM 0)
include(BuildExeAndLib)
