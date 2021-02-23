#ifndef _Derive_H_
#define _Derive_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#ifdef PELEC_USE_MASA
#include <masa.h>
using namespace MASA;
#endif

void pc_dervelx(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dervely(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dervelz(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dermagvel(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dermagmom(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derkineng(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dereint1(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dereint2(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derlogden(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derspec(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dermagvort(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derdivu(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derenstrophy(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dernull(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dermolefrac(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dersoundspeed(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derentropy(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dermachnumber(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derpres(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dertemp(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derspectrac(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level,
  const int idx);

void pc_derradialvel(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

#ifdef PELEC_USE_MASA
void pc_derrhommserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derummserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dervmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derwmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derpmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);
#endif
#endif
