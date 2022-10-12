#ifndef SparseData_H
#define SparseData_H

#include <AMReX_Vector.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
getIndex(const int i, const int comp, const int size)
{
  return comp * size + i;
}

// SparseData is a templated data holder defined over a vector of Cell objects.
template <class T, class Cell>
class SparseData
{
public:
  SparseData() = default;

  ~SparseData();

  // Defining constructor.  Specifies the irregular domain
  // and the number of data components per index. The
  // contents are uninitialized.  Calls full define function.
  SparseData(const amrex::Gpu::DeviceVector<Cell>& region, int nComp);

  // Full define function.  Specifies the irregular domain and the
  // number of data components per index.  The contents
  // are uninitialized.  If it has previously been defined, the old definition
  // data is overwritten and lost.
  void define(const amrex::Gpu::DeviceVector<Cell>& region, int nComp);

  AMREX_FORCE_INLINE T* dataPtr(int comp = 0)
  {
    return &(m_data.data()[getIndex(0, comp, m_region_size)]);
  }

  const T* dataPtr(int comp = 0) const
  {
    return &(m_data.data()[getIndex(0, comp, m_region_size)]);
  }

  void setVal(const T& val);

  void setVal(const T& val, int comp, int ncomp = 1);

  void merge(
    const SparseData& thdlocal,
    int comp,
    int ncomp,
    const amrex::Gpu::DeviceVector<int> mask);

  int numPts() const { return m_region_size; }

  int nComp() const { return m_ncomp; }

private:
  int m_ncomp = 0;
  int m_region_size = 0;
  amrex::Gpu::DeviceVector<Cell> m_region;
  amrex::Gpu::DeviceVector<T> m_data;
};

template <class T, class Cell>
AMREX_FORCE_INLINE
SparseData<T, Cell>::SparseData(
  const amrex::Gpu::DeviceVector<Cell>& _region, int _nComp)
{
  define(_region, _nComp);
  amrex::Print() << "Init SparseData with ncomp = " << _nComp << std::endl;
}

template <class T, class Cell>
AMREX_FORCE_INLINE SparseData<T, Cell>::~SparseData()
{
  m_data.clear();
  m_ncomp = 0;
  m_region.clear();
  m_region_size = 0;
}

template <class T, class Cell>
AMREX_FORCE_INLINE void
SparseData<T, Cell>::define(
  const amrex::Gpu::DeviceVector<Cell>& _region, int _nComp)
{
  // Note, could use ref/ptr but dangerous, however this copy might be expensive
  m_region = _region;
  m_region_size = m_region.size();
  m_ncomp = _nComp;
  m_data.resize(numPts() * m_ncomp);
}

template <class T, class Cell>
AMREX_FORCE_INLINE void
SparseData<T, Cell>::setVal(const T& val)
{
  for (int i = 0; i < m_ncomp; ++i) {
    setVal(val, i);
  }
}

template <class T, class Cell>
AMREX_FORCE_INLINE void
SparseData<T, Cell>::setVal(const T& val, int comp, int ncomp)
{
  AMREX_ASSERT(comp + ncomp <= m_ncomp);
  auto* d_m_data = m_data.data();
  const int captured_m_region_size = m_region_size;
  amrex::ParallelFor(m_region_size, [=] AMREX_GPU_DEVICE(int i) {
    for (int n = 0; n < ncomp; ++n) {
      d_m_data[getIndex(i, comp + n, captured_m_region_size)] = val;
    }
  });
}

template <class T, class Cell>
AMREX_FORCE_INLINE void
SparseData<T, Cell>::merge(
  const SparseData& thdlocal,
  int comp,
  int ncomp,
  const amrex::Gpu::DeviceVector<int> mask)
{
  AMREX_ASSERT(comp + ncomp <= m_ncomp);
  const int captured_m_region_size = m_region_size;
  auto* d_data = m_data.data();
  auto* d_thdlocal_data = thdlocal.m_data.data();
  const auto* d_mask = mask.data();
  amrex::ParallelFor(captured_m_region_size, [=] AMREX_GPU_DEVICE(int i) {
    if (d_mask[i] != 0) {
      for (int n = 0; n < ncomp; ++n) {
#ifdef AMREX_USE_OMP
#pragma omp atomic write
#endif
        d_data[getIndex(i, comp + n, captured_m_region_size)] =
          d_thdlocal_data[getIndex(i, comp + n, captured_m_region_size)];
      }
    }
  });
}

#endif
