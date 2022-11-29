#ifndef __INCLUDE_GUARD_Utility_hh__
#define __INCLUDE_GUARD_Utility_hh__
#include <iomanip>      // std::setprecision(), std::put:time
#include <iostream>     // cout
#include <fstream>      // file input/output
#include <cmath>        // signbit
#include <algorithm>    // clamp
#include <sstream>      // stringstream
#include <filesystem>   // folder/file management
#include <unistd.h>     // usleep
#include <chrono>       // time
#include "eigen/Eigen/Dense"


// Fast itneger exponentiation.
template<int N>
inline __attribute__((always_inline)) double MyPow(double a)
{ return a * MyPow<N-1>(a); }
template<>
inline __attribute__((always_inline)) double MyPow<0>(double a)
{ return 1; }

// Get sign of T
template <typename T>
inline __attribute__((always_inline)) int sgn(T val)
{ return (T(0) < val) - (val < T(0)); }



// Marker for debugging
inline __attribute__((always_inline)) void Marker(std::string name="", bool newline=true)
{
    static int i = 0;
    std::cout << name << ": " << i << std::endl;
    if(newline)
        std::cout << std::endl;
    i++;
}



// Number of digits of given number.
inline __attribute__((always_inline)) int numDigits(int number)
{
    int digits = 0;
    while (number)
	{
        number /= 10;
        digits++;
    }
    return digits;
}



// Converts frame number to correct format, e.g. 1 -> 0001.
inline __attribute__((always_inline)) std::string FrameNumber(unsigned int f)
{
	std::string frameNumber;
	int maxDigits = 4;
	int fDigits   = numDigits(f);

	for(int i=0; i<(maxDigits-fDigits); i++)
		frameNumber += "0";
    if(f!=0)
    	frameNumber += std::to_string(f);

	return frameNumber;
}



// Add leading space if positive.
inline __attribute__((always_inline)) std::string Format(const double n, const int precision=6)
{
    std::string output;

    // Leading space if positive:
    if(n == 0 and std::signbit(n) == 0){output = " ";}
    else if(n == 0 and std::signbit(n) == 1){output = "";}
	else if(n >= 0){output = " ";}
	else   {output = "";}

    // Leading space if 2 digit number:
    //if(abs(n)<100){output += " ";}
    // Leading space if 1 digit number:
    //if(abs(n)<10){output += " ";}

    // Number of decimal digits:
    std::ostringstream ss;
    ss.precision(precision);
    ss << std::fixed << n; 
    output += ss.str();

    return output;    
}



// Print formatted double.
inline __attribute__((always_inline)) void PrintDouble(const double d, std::string name, bool newline=false, const int precision=6)
{
    std::cout << name << " = "<< Format(d,precision) << "\n";
    if(newline) std::cout << "\n";
}



// Sleep for given amount of seconds.
inline __attribute__((always_inline)) void sleep(double seconds)
{
    int microsecondsInSecond = 1000000;
    usleep(seconds*microsecondsInSecond);
}



// Exit programm with Error Message.
inline __attribute__((always_inline)) void exit_on_error(const char* const msg="")
{
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(errno);
}



// Minimum value of given array
inline __attribute__((always_inline)) double MinValue(const double* const array, int length)
{
    double min = 1e16;
    for(int i=0; i<length; i++)
        if(array[i]<min)
            min = array[i];
    return min;
}
// Maximum value of given array
inline __attribute__((always_inline)) double MaxValue(const double* const array, int length)
{
    double max = -1e16;
    for(int i=0; i<length; i++)
        if(max<array[i])
            max = array[i];
    return max;
}


// Inverse lerp.
inline __attribute__((always_inline)) double Map01(double x, double min, double max)
{ return std::clamp((x - min) / (max - min), 0.0, 1.0); }



// Print formatted Matrix.
inline __attribute__((always_inline)) void PrintMat(const Eigen::MatrixXd& A, int precision=6)
{
    for(int i=0; i<A.cols(); i++)
    {
        std::cout << "(" << Format(A(i,0),precision);
        for(int j=1; j<A.rows(); j++)
            std::cout << ", " << Format(A(i,j),precision);
        std::cout << ")" << std::endl;
    }
}



// File management:
// Creates directory in current working directory.
// Nested directories possible, e.g: path="data/output/xValues"
inline __attribute__((always_inline)) bool CreateDirectory(std::string path)
{
	namespace fs = std::filesystem;
    fs::path currentPath = fs::current_path();
	std::string directoryPath = currentPath/path;
    return fs::create_directories(directoryPath);
}
#endif //__INCLUDE_GUARD_Utility_hh__