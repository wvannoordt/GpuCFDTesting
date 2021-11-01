#ifndef STATIC_LOOP_H
#define STATIC_LOOP_H
#include <type_traits>
#include "CudaHeaders.h"

template <const int start, const int finish, typename loopCall> static _f_hybrid inline void static_for(const loopCall& loopObj)
{
	if constexpr (start < finish)
	{
		loopObj(std::integral_constant<int, start>{});
		static_for<start+1, finish>(loopObj);
	}
};

#endif