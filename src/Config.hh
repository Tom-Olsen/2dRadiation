#ifndef __INCLUDE_GUARD_Config_hh__
#define __INCLUDE_GUARD_Config_hh__
#include "Utility.hh"

// Input system:
enum StreamingType
{
    FlatFixed,
    FlatAdaptive,
    CurvedFixed,
    CurvedAdaptive,
    GeodesicFixed
};
inline std::string StreamingName(int n)
{
    std::string name("unknown");
    switch (n)
    {
    case 0:
        name = "FlatFixed";
        break;
    case 1:
        name = "FlatAdaptive";
        break;
    case 2:
        name = "CurvedFixed";
        break;
    case 3:
        name = "CurvedAdaptive";
        break;
    case 4:
        name = "GeodesicFixed";
        break;
    default:
        ExitOnError("Invalid StreamingType");
    }
    return name;
}
enum InitialDataType
{
    Moments,
    Intensities
};

struct Config
{
    std::string name;
    double t0;
    double simTime;
    double writePeriod;
    bool updateFourierHarmonics;
    bool keepSourceNodesActive;
    bool writeData;
    bool printToTerminal;
    StreamingType streamingType;
    InitialDataType initialDataType;
};
#endif //__INCLUDE_GUARD_Config_hh__