#pragma once

#include <complex>

using complex = std::complex<double>;

const complex I(0,1); 
constexpr double Pi=4*atan(1);

const int initial_random=time(NULL); 
extern int INT_RANDOM;
