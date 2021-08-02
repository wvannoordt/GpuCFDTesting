#ifndef PRINT_DEBUG_H
#define PRINT_DEBUG_H
#include <iostream>
#include <ostream>

template <typename T> void print_recursive(std::ostream& stm, T t)
{
    stm << t << std::endl;
}

template <typename T, typename... Ts> void print_recursive(std::ostream& stm, T t, Ts... ts)
{
    stm << t << " ";
    print_recursive(stm, ts...);
}

template <typename... Ts> void print(Ts... ts)
{
    print_recursive(std::cout, ts...);
}

template <typename... Ts> void PrintToStream(std::ostream& stm, Ts... ts)
{
    print_recursive(stm, ts...);
}

#endif