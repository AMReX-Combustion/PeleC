#ifndef _SparseData_H_
#define _SparseData_H_

#include <vector>
///
/**
   SparseData is a templated data holder defined over a vector of Cell objects.
*/
    template <class T, class Cell>
    class SparseData
    {
    public:
        typedef T value_type;
        ///
        /**
           Null constructor.  
        */
        SparseData() {}

        /**
           Destructor 
        */
        ~SparseData();

        ///
        /**
           Defining constructor.  Specifies the irregular domain
           and the number of data components per index. The
           contents are uninitialized.  Calls full define function.
        */
        SparseData(const std::vector<Cell>& region,
                   int                      nComp);

        ///
        /**
           Full define function.  Specifies the irregular domain and the
           number of data components per index.  The contents
           are uninitialized.  If it has previously been defined, the old definition
           data is overwritten and lost.
        */
        void define(const std::vector<Cell>& region,
                    int                      nComp);

        const std::vector<Cell>&  getRegion() const;

        inline T* dataPtr(int comp=0) {return &(m_data[getIndex(0,comp)]);}

        const T* dataPtr(int comp=0) const {return &(m_data[getIndex(0,comp)]);}
  
        void setVal(const T& val);

        void setVal(const T& val, int comp, int ncomp=1);
        void merge(const SparseData& thdlocal, int comp, int ncomp, const std::vector<int> mask);

        ///
        T& operator() (int i, int comp) {return m_data[getIndex(i,comp)];}

        ///
        const T& operator() (int i, int comp) const {return m_data[getIndex(i,comp)];}

        int numPts() const {return m_region.size();}

        int nComp() const {return m_ncomp;}

    private:

    protected:
        int getIndex(int i, int comp) const {return comp*m_region.size() + i;}

        int m_ncomp = 0;
        std::vector<Cell> m_region;
        std::vector<T> m_data;
  };

template <class T, class Cell> inline
SparseData<T,Cell>::SparseData(const std::vector<Cell>& _region,
                               int                      _nComp)
{
    define(_region,_nComp);
    amrex::Print() << "Init SparseData with ncomp = " << _nComp <<std::endl;
}

template <class T, class Cell> inline
SparseData<T,Cell>::~SparseData()
{
    m_data.clear();
    m_ncomp = 0;
    m_region.clear();
}

template <class T, class Cell> inline
void
SparseData<T,Cell>::define(const std::vector<Cell>& _region,
                           int                      _nComp)
{
    // Note, could use ref/ptr but dangerous, however this copy might be expensive
    m_region = _region;
    m_ncomp = _nComp;
    m_data.resize(numPts() * m_ncomp);
}

template <class T, class Cell> inline
void
SparseData<T,Cell>::setVal(const T& val)
{
    for (int i=0; i<m_ncomp; ++i)
    {
        setVal(val,i);
    }
}

template <class T, class Cell> inline
void
SparseData<T,Cell>::setVal(const T& val, int comp, int ncomp)
{
    BL_ASSERT(comp+ncomp <= m_ncomp);
    for (int n=0; n<ncomp; ++n)
    {
      for (int i=0; i<m_region.size(); ++i)
      {
        m_data[getIndex(i,comp+n)] = val;
      }
    }
}

template <class T, class Cell> inline
void
SparseData<T,Cell>::merge(const SparseData& thdlocal, int comp, int ncomp, const std::vector<int> mask)
{
    BL_ASSERT(comp+ncomp <= m_ncomp);

    for (int n=0; n<ncomp; ++n)
    {
        for (int i=0; i<m_region.size(); ++i)
        {
            if (mask[i]) {
#ifdef _OPENMP
#pragma omp atomic write
#endif                
                m_data[getIndex(i,comp+n)] = thdlocal.m_data[getIndex(i,comp+n)];
            }
        }
    }
}

#endif
