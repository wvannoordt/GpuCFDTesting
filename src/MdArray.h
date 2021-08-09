#ifndef CMF_MD_ARRAY_H
#define CMF_MD_ARRAY_H

#include <type_traits>
#include "CallError.h"
#include "CudaHeaders.h"

template <typename arType, const int arRank = 1> struct MdArray
{
	arType* data;
	int rank = arRank;
	int dims[arRank];
    int idxCoeff[arRank];
    size_t totalSize;
    
    MdArray(const MdArray& rhs)
    {
        data = rhs.data;
        rank = rhs.rank;
        totalSize = rhs.totalSize;
        for (int i = 0; i < arRank; i++)
        {
            dims[i] = rhs.dims[i];
            idxCoeff[i] = rhs.idxCoeff[i];
        }
    }
    
    MdArray& operator =(const MdArray& rhs)
    {
        data = rhs.data;
        rank = rhs.rank;
        totalSize = rhs.totalSize;
        for (int i = 0; i < arRank; i++)
        {
            dims[i] = rhs.dims[i];
            idxCoeff[i] = rhs.idxCoeff[i];
        }
        return *this;
    }
    
    template <typename arType2, const int arRank2 = 1> MdArray<arType2, arRank2> ReCast(int index)
    {
        //Expand this to reCast in rank as well?
        static_assert(arRank2==arRank, "Incorrect rank when performing ReCast");
        MdArray<arType2, arRank2> output;
        for (int i = 0; i < arRank2; i++) output.dims[i] = dims[i];
        output.dims[index] = (output.dims[index]*sizeof(arType))/sizeof(arType2);
        output.data = (arType2*)data;
        output.idxCoeff[0] = 1;
        for (int i = 1; i < arRank; i++)
        {
            output.idxCoeff[i] = output.idxCoeff[i-1]*output.dims[i-1];
        }
        return output;
    }

	void Ralloc(int lev)
	{
		totalSize = 0;
	}

	template <typename Ts> void Ralloc(int lev, Ts t)
	{
		dims[lev] = t;
        idxCoeff[0] = 1;
        for (int i = 1; i < arRank; i++)
        {
            idxCoeff[i] = idxCoeff[i-1]*dims[i-1];
        }
        totalSize = 1;
        for (int i = 0; i < arRank; i++)
        {
            totalSize *= dims[i];
        }
	}

	template <typename T, typename... Ts> void Ralloc(int lev, T t, Ts... ts)
	{
		static_assert(std::is_integral<T>::value, "Integral type required for dimension initialization.");
		dims[lev] = t;
		Ralloc(lev+1, ts...);
	}

	template <typename... Ts> MdArray(Ts... ts)
	{
		Ralloc(0, ts...);
	}

	template <typename... Ts> MdArray(arType* ptr, Ts... ts)
	{
        data = ptr;
		Ralloc(0, ts...);
	}
	
    template <typename T> _f_hybrid inline arType * idxC(int lev, T t)
    {
        return data+idxCoeff[lev]*t;
    }

    template <typename T, typename... Ts> _f_hybrid inline arType * idxC(int lev, T t, Ts... ts)
    {
        static_assert(std::is_integral<T>::value, "Integral type required for indexing");
        return idxC(lev+1, ts...) + t*idxCoeff[lev];
    }

    template <typename... Ts> _f_hybrid inline arType & operator () (Ts... ts)
    {
        static_assert(sizeof...(Ts)==arRank, "Incorrect rank in array index");
        return *(idxC(0, ts...));
    }
    
    template <typename T> _f_hybrid inline size_t offsetInternal(int lev, T t)
    {
        return idxCoeff[lev]*t;
    }

    template <typename T, typename... Ts> _f_hybrid inline size_t offsetInternal(int lev, T t, Ts... ts)
    {
        static_assert(std::is_integral<T>::value, "Integral type required for indexing");
        return offsetInternal(lev+1, ts...) + t*idxCoeff[lev];
    }

    template <typename... Ts> _f_hybrid inline size_t offset(Ts... ts)
    {
        static_assert(sizeof...(Ts)==arRank, "Incorrect rank in array index");
        return offsetInternal(0, ts...);
    }
};

#endif