// #ifndef CONSTMATH_H
// #define CONSTMATH_H
// #include "CudaHeaders.h"
// #include <type_traits>
// namespace ConstMath
// {
//     const double pi = 3.14159265359;
//     const double epsilon = 1e-10;
//     const int maxIterations = 50;
//     template <typename numericType> constexpr _f_hybrid numericType mod(const numericType a, const numericType b)
//     {
//         if (a<0.0) return mod(a+b, b);
//         else if (a>=b) return mod(a-b, b);
//         else return a;
//     }
//     template <typename numericType, const int ITERLEVEL> constexpr _f_hybrid numericType powfac(const numericType x)
//     {
//         if constexpr (ITERLEVEL<=0) return 1.0;
//         else return (x/((numericType)ITERLEVEL))*powfac<numericType, ITERLEVEL-1>(x);
//     }
//     template <typename numericType> constexpr _f_hybrid numericType abs(const numericType x)
//     {
//         return (x < 0.0) ? (-x) : (x);
//     }
//     template <typename numericType, const int ITERLEVEL> constexpr _f_hybrid numericType cos_recurse(const numericType eps, const numericType sum, const numericType xmod)
//     {
//         static_assert(ITERLEVEL>=0, "Invalid iteration index (cos_recurse)!");
//         const int isOdd = ITERLEVEL&1;
//         const numericType sign = 1.0-2.0*isOdd;
//         const numericType newTerm = powfac<numericType, ITERLEVEL>(xmod);
//         const numericType incremented = sum+sign*newTerm;
//         if constexpr (abs<numericType>(newTerm)<eps) return incremented;
//         if constexpr (ITERLEVEL > maxIterations) return sum;
//         else return incremented+cos_recurse<numericType, ITERLEVEL+1>(eps, sum, xmod);
//     }
//     template <typename numericType> constexpr _f_hybrid numericType cos(const numericType x)
//     {
//         numericType xmod = mod(x, 2.0*pi);
//         numericType sum = 0.0;
//         return cos_recurse<numericType, 0>(epsilon, sum, xmod);
//     }
// }
// 
// #endif