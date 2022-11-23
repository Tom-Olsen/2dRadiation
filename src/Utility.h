#ifndef __INCLUDE_GUARD_Utility_h__
#define __INCLUDE_GUARD_Utility_h__
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



// Get sign of int/float/double.
template <typename T>
int sgn(T val);

// marker to debug code.
void Marker(std::string name="", bool newline=false);

// return number of digits in given int.
int numDigits(int number);

// Converts frame number to correct format, e.g. 1 -> 0001.
std::string FrameNumber(unsigned int);

// Add leading space if positive.
std::string Format(const double n, int precision=6);

// Print double.
void PrintDouble(const double d, std::string name, bool newline=false, int precision=6);

// Sleep for given amount of seconds.
void sleep(double seconds);

// Exit Programm with Error Message.
void exit_on_error(const char* const msg = "");

// Minimum value of given array
double MinValue(const double* const array, int length);

// Maximum value of given array
double MaxValue(const double* const array, int length);

// Inverse Lerp: returns 0 if x=min, 1 if x=max, and linear interpolated everywhere else.
double Map01(double x, double min, double max);



// Debugging:
void PrintMat(const Eigen::MatrixXd& A);




// -------------------- File management --------------------
// Creates directory in current working directory.
// Nested directories possible, e.g: path="data/output/xValues"
bool CreateDirectory(std::string path);
// ---------------------------------------------------------


// Timer:
template<int ID>
struct Timer
{
    static double time;
    static std::string name;
    static std::chrono::steady_clock::time_point p;

    static void Start();
    static void End();
    static void Reset(std::string name_);
    static void Print();
};

#endif //__INCLUDE_GUARD_Utility_h__