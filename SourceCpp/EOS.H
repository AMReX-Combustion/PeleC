#ifndef _EOS_H_
#define _EOS_H_

#ifdef PELEC_EOS_FUEGO
#include <EOS_FUEGO.H>
#elif PELEC_EOS_GAMMA
#include <EOS_GAMMA.H>
#endif

void pc_eos_init();

#endif
