/* problem-specific stuff goes here */
#include <AMReX_BLFort.H>

#ifndef _PeleC_Problem_F_H_
#define _PeleC_Problem_F_H_
#include <AMReX_BLFort.H>

BL_FORT_PROC_DECL(PROBLEM_CHECKPOINT,problem_checkpoint)
  (int * int_dir_name, int * len);

BL_FORT_PROC_DECL(PROBLEM_RESTART,problem_restart)
  (int * int_dir_name, int * len);

#endif


