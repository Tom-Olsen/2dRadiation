#ifndef __INCLUDE_GUARD_ControlFlow_hh__
#define __INCLUDE_GUARD_ControlFlow_hh__

// Coordinate Systems:
class xy {};
class rph {};

// Frames:
class IF {};
class LF {};
class Tetrad {}; // the tetrad is actually a mixed tensor having one IF and one LF index, e.g. e^mu'_mu

// Profiling macros:
#define PROFILING 1
#if PROFILING
    #define PROFILE_SCOPE(name) Profiler::Timer timer##__LINE__(name)
    #define PROFILE_FUNCTION() PROFILE_SCOPE(__PRETTY_FUNCTION__)
#else
    #define PROFILE_FUNCTION()
#endif

#endif //__INCLUDE_GUARD_ControlFlow_hh__