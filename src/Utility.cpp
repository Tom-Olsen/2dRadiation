#include "Utility.h"



template <typename T>
int sgn(T val)
{return (T(0) < val) - (val < T(0));}
template int sgn<int>(int val);
template int sgn<float>(float val);
template int sgn<double>(double val);



void Marker(std::string name, bool newline)
{
    static int i = 0;
    std::cout << name << ": " << i << std::endl;
    if(newline)
        std::cout << std::endl;
    i++;
}



// Returns number of digits of given number.
int numDigits(int number)
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
std::string FrameNumber(unsigned int f)
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
std::string Format(const double n, const int precision)
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



// Print double.
void PrintDouble(const double d, std::string name, bool newline, const int precision)
{
    std::cout << name << " = "<< Format(d,precision) << std::endl;
    if(newline)
        std::cout << std::endl;
}



// Sleep for given amount of seconds.
void sleep(double seconds)
{
    int microsecondsInSecond = 1000000;
    usleep(seconds*microsecondsInSecond);
}



// Exit Programm with Error Message.
void exit_on_error(const char* const msg)
{
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(errno);
}



// Minimum value of given array
double MinValue(const double* const array, int length)
{
    double min = 1e16;
    for(int i=0; i<length; i++)
        if(array[i]<min)
            min = array[i];
    return min;
}
// Maximum value of given array
double MaxValue(const double* const array, int length)
{
    double max = -1e16;
    for(int i=0; i<length; i++)
        if(max<array[i])
            max = array[i];
    return max;
}


// Inverse lerp.
double Map01(double x, double min, double max)
{
    return std::clamp((x - min) / (max - min), 0.0, 1.0);
}



// Debugging:
void PrintMat(const Eigen::MatrixXd& A)
{
    for(int i=0; i<A.cols(); i++)
    {
        std::cout << "(" << Format(A(i,0));
        for(int j=1; j<A.rows(); j++)
            std::cout << ", " << Format(A(i,j));
        std::cout << ")" << std::endl;
    }
}



// File management:
bool CreateDirectory(std::string path)
{
	namespace fs = std::filesystem;
    fs::path currentPath = fs::current_path();
	std::string directoryPath = currentPath/path;
    return fs::create_directories(directoryPath);
}



// Timer class:
template<int ID>
void Timer<ID>::Start()
{
    p = std::chrono::steady_clock::now();
}
template<int ID>
void Timer<ID>::End()
{
    using duration = std::chrono::duration<double>;
    auto d = std::chrono::steady_clock::now() - p;
    time += std::chrono::duration_cast<duration>(d).count();
}
template<int ID>
void Timer<ID>::Reset(std::string name_)
{
    name = name_;
    time = 0;
    p = std::chrono::steady_clock::now();
}
template<int ID>
void Timer<ID>::Print()
{
    std::cout << "Timer, " << name << "(" << ID << "): " << time << "s" << std::endl;
}

template<int ID>
double Timer<ID>::time;
template<int ID>
std::chrono::steady_clock::time_point Timer<ID>::p;
template<int ID>
std::string Timer<ID>::name;

template class Timer< 0>;   template class Timer<10>;   template class Timer<20>;   template class Timer<30>;   template class Timer<40>;
template class Timer< 1>;   template class Timer<11>;   template class Timer<21>;   template class Timer<31>;   template class Timer<41>;
template class Timer< 2>;   template class Timer<12>;   template class Timer<22>;   template class Timer<32>;   template class Timer<42>;
template class Timer< 3>;   template class Timer<13>;   template class Timer<23>;   template class Timer<33>;   template class Timer<43>;
template class Timer< 4>;   template class Timer<14>;   template class Timer<24>;   template class Timer<34>;   template class Timer<44>;
template class Timer< 5>;   template class Timer<15>;   template class Timer<25>;   template class Timer<35>;   template class Timer<45>;
template class Timer< 6>;   template class Timer<16>;   template class Timer<26>;   template class Timer<36>;   template class Timer<46>;
template class Timer< 7>;   template class Timer<17>;   template class Timer<27>;   template class Timer<37>;   template class Timer<47>;
template class Timer< 8>;   template class Timer<18>;   template class Timer<28>;   template class Timer<38>;   template class Timer<48>;
template class Timer< 9>;   template class Timer<19>;   template class Timer<29>;   template class Timer<39>;   template class Timer<49>;

template class Timer<50>;   template class Timer<60>;   template class Timer<70>;   template class Timer<80>;   template class Timer<90>;
template class Timer<51>;   template class Timer<61>;   template class Timer<71>;   template class Timer<81>;   template class Timer<91>;
template class Timer<52>;   template class Timer<62>;   template class Timer<72>;   template class Timer<82>;   template class Timer<92>;
template class Timer<53>;   template class Timer<63>;   template class Timer<73>;   template class Timer<83>;   template class Timer<93>;
template class Timer<54>;   template class Timer<64>;   template class Timer<74>;   template class Timer<84>;   template class Timer<94>;
template class Timer<55>;   template class Timer<65>;   template class Timer<75>;   template class Timer<85>;   template class Timer<95>;
template class Timer<56>;   template class Timer<66>;   template class Timer<76>;   template class Timer<86>;   template class Timer<96>;
template class Timer<57>;   template class Timer<67>;   template class Timer<77>;   template class Timer<87>;   template class Timer<97>;
template class Timer<58>;   template class Timer<68>;   template class Timer<78>;   template class Timer<88>;   template class Timer<98>;
template class Timer<59>;   template class Timer<69>;   template class Timer<79>;   template class Timer<89>;   template class Timer<99>;


template class Timer<100>;   template class Timer<110>;   template class Timer<120>;   template class Timer<130>;   template class Timer<140>;
template class Timer<101>;   template class Timer<111>;   template class Timer<121>;   template class Timer<131>;   template class Timer<141>;
template class Timer<102>;   template class Timer<112>;   template class Timer<122>;   template class Timer<132>;   template class Timer<142>;
template class Timer<103>;   template class Timer<113>;   template class Timer<123>;   template class Timer<133>;   template class Timer<143>;
template class Timer<104>;   template class Timer<114>;   template class Timer<124>;   template class Timer<134>;   template class Timer<144>;
template class Timer<105>;   template class Timer<115>;   template class Timer<125>;   template class Timer<135>;   template class Timer<145>;
template class Timer<106>;   template class Timer<116>;   template class Timer<126>;   template class Timer<136>;   template class Timer<146>;
template class Timer<107>;   template class Timer<117>;   template class Timer<127>;   template class Timer<137>;   template class Timer<147>;
template class Timer<108>;   template class Timer<118>;   template class Timer<128>;   template class Timer<138>;   template class Timer<148>;
template class Timer<109>;   template class Timer<119>;   template class Timer<129>;   template class Timer<139>;   template class Timer<149>;

template class Timer<150>;   template class Timer<160>;   template class Timer<170>;   template class Timer<180>;   template class Timer<190>;
template class Timer<151>;   template class Timer<161>;   template class Timer<171>;   template class Timer<181>;   template class Timer<191>;
template class Timer<152>;   template class Timer<162>;   template class Timer<172>;   template class Timer<182>;   template class Timer<192>;
template class Timer<153>;   template class Timer<163>;   template class Timer<173>;   template class Timer<183>;   template class Timer<193>;
template class Timer<154>;   template class Timer<164>;   template class Timer<174>;   template class Timer<184>;   template class Timer<194>;
template class Timer<155>;   template class Timer<165>;   template class Timer<175>;   template class Timer<185>;   template class Timer<195>;
template class Timer<156>;   template class Timer<166>;   template class Timer<176>;   template class Timer<186>;   template class Timer<196>;
template class Timer<157>;   template class Timer<167>;   template class Timer<177>;   template class Timer<187>;   template class Timer<197>;
template class Timer<158>;   template class Timer<168>;   template class Timer<178>;   template class Timer<188>;   template class Timer<198>;
template class Timer<159>;   template class Timer<169>;   template class Timer<179>;   template class Timer<189>;   template class Timer<199>;


template class Timer<200>;   template class Timer<210>;   template class Timer<220>;   template class Timer<230>;   template class Timer<240>;
template class Timer<201>;   template class Timer<211>;   template class Timer<221>;   template class Timer<231>;   template class Timer<241>;
template class Timer<202>;   template class Timer<212>;   template class Timer<222>;   template class Timer<232>;   template class Timer<242>;
template class Timer<203>;   template class Timer<213>;   template class Timer<223>;   template class Timer<233>;   template class Timer<243>;
template class Timer<204>;   template class Timer<214>;   template class Timer<224>;   template class Timer<234>;   template class Timer<244>;
template class Timer<205>;   template class Timer<215>;   template class Timer<225>;   template class Timer<235>;   template class Timer<245>;
template class Timer<206>;   template class Timer<216>;   template class Timer<226>;   template class Timer<236>;   template class Timer<246>;
template class Timer<207>;   template class Timer<217>;   template class Timer<227>;   template class Timer<237>;   template class Timer<247>;
template class Timer<208>;   template class Timer<218>;   template class Timer<228>;   template class Timer<238>;   template class Timer<248>;
template class Timer<209>;   template class Timer<219>;   template class Timer<229>;   template class Timer<239>;   template class Timer<249>;

template class Timer<250>;   template class Timer<260>;   template class Timer<270>;   template class Timer<280>;   template class Timer<290>;
template class Timer<251>;   template class Timer<261>;   template class Timer<271>;   template class Timer<281>;   template class Timer<291>;
template class Timer<252>;   template class Timer<262>;   template class Timer<272>;   template class Timer<282>;   template class Timer<292>;
template class Timer<253>;   template class Timer<263>;   template class Timer<273>;   template class Timer<283>;   template class Timer<293>;
template class Timer<254>;   template class Timer<264>;   template class Timer<274>;   template class Timer<284>;   template class Timer<294>;
template class Timer<255>;   template class Timer<265>;   template class Timer<275>;   template class Timer<285>;   template class Timer<295>;
template class Timer<256>;   template class Timer<266>;   template class Timer<276>;   template class Timer<286>;   template class Timer<296>;
template class Timer<257>;   template class Timer<267>;   template class Timer<277>;   template class Timer<287>;   template class Timer<297>;
template class Timer<258>;   template class Timer<268>;   template class Timer<278>;   template class Timer<288>;   template class Timer<298>;
template class Timer<259>;   template class Timer<269>;   template class Timer<279>;   template class Timer<289>;   template class Timer<299>;