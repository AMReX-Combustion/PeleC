#ifndef _FILTER_H_
#define _FILTER_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>

#include "Constants.H"
#include "Utilities.H"

#ifdef _OPENMP
#include <omp.h>
#endif

// Filter types
enum filter_types {
  no_filter = 0,
  box,                           // 1
  gaussian,                      // 2
  box_3pt_approx,                // 3
  box_5pt_approx,                // 4
  box_3pt_optimized_approx,      // 5
  box_5pt_optimized_approx,      // 6
  gaussian_3pt_approx,           // 7
  gaussian_5pt_approx,           // 8
  gaussian_3pt_optimized_approx, // 9
  gaussian_5pt_optimized_approx, // 10
  num_filter_types
};

class Filter
{

public:
  // Default constructor
  Filter(const int type = box, const int fgr = 2) : _type(type), _fgr(fgr)
  {

    switch (_type) {

    case box:
      set_box_weights();
      break;

    case gaussian:
      set_gaussian_weights();
      break;

    case box_3pt_approx:
    case gaussian_3pt_approx: // same as box_3pt_approx
      set_box_3pt_approx_weights();
      break;

    case box_5pt_approx:
      set_box_5pt_approx_weights();
      break;

    case box_3pt_optimized_approx:
      set_box_3pt_optimized_approx_weights();
      break;

    case box_5pt_optimized_approx:
      set_box_5pt_optimized_approx_weights();
      break;

    case gaussian_5pt_approx:
      set_gaussian_5pt_approx_weights();
      break;

    case gaussian_3pt_optimized_approx:
      set_gaussian_3pt_optimized_approx_weights();
      break;

    case gaussian_5pt_optimized_approx:
      set_gaussian_5pt_optimized_approx_weights();
      break;

    case no_filter:
    default:
      _fgr = 1;
      _ngrow = 0;
      _nweights = 2 * _ngrow + 1;
      _weights.resize(_nweights);
      _weights[0] = 1.;
      break;

    } // end switch
  }

  // Default destructor
  ~Filter() = default;

  int get_filter_ngrow() const { return _ngrow; }

  void apply_filter(const amrex::MultiFab& in, amrex::MultiFab& out);

  void apply_filter(
    const amrex::MultiFab& in,
    amrex::MultiFab& out,
    const int nstart,
    const int ncnt);

  void apply_filter(
    const amrex::Box& cbox, const amrex::FArrayBox& in, amrex::FArrayBox& out);

  void apply_filter(
    const amrex::Box& cbox,
    const amrex::FArrayBox& in,
    amrex::FArrayBox& out,
    const int nstart,
    const int ncnt);

  void apply_filter(
    const amrex::Box& box,
    const amrex::FArrayBox& in,
    amrex::FArrayBox& out,
    const int nstart,
    const int ncnt,
    const int ncomp);

private:
  int _type;
  int _fgr;
  int _ngrow;
  int _nweights;
  amrex::Vector<amrex::Real> _weights;

  void set_box_weights();

  void set_gaussian_weights();

  void set_box_3pt_approx_weights();

  void set_box_5pt_approx_weights();

  void set_box_3pt_optimized_approx_weights();

  void set_box_5pt_optimized_approx_weights();

  void set_gaussian_5pt_approx_weights();

  void set_gaussian_3pt_optimized_approx_weights();

  void set_gaussian_5pt_optimized_approx_weights();
};

#endif /*_FILTER_H_*/
