//
// Basic instrumentation profiler by Cherno

// Usage: include this header file somewhere in your code (eg. precompiled header), and then use like:
//
// Profiler::Get().BeginSession("Session Name");        // Begin session 
// {
//     ProfileTimer timer("Profiled Scope Name");   // Place code like this in scopes you'd like to include in profiling
//     // Code
// }
// Profiler::Get().EndSession();                        // End Session
//
// You will probably want to macro-fy this, to switch on/off easily and use things like __FUNCSIG__ for the profile name.
//
#ifndef __INCLUDE_INSTRUMENTOR_hh__
#define __INCLUDE_INSTRUMENTOR_hh__
#include <string>
#include <chrono>
#include <algorithm>
#include <fstream>
// #include <mutex>
// #include <thread>



struct ProfileResult
{
    std::string Name;
    long long Start, End;
    uint32_t ThreadID;
};



struct ProfileSession
{
    std::string Name;
};



class Profiler
{
private:
    ProfileSession* m_CurrentSession;
    std::ofstream m_OutputStream;
    int m_ProfileCount;
    // std::mutex m_lock;

public:
    Profiler() : m_CurrentSession(nullptr), m_ProfileCount(0) {}

    void BeginSession(const std::string& name, const std::string& filepath = "output/results.json")
    {
        m_OutputStream.open(filepath);
        WriteHeader();
        m_CurrentSession = new ProfileSession{ name };
    }

    void EndSession()
    {
        WriteFooter();
        m_OutputStream.close();
        delete m_CurrentSession;
        m_CurrentSession = nullptr;
        m_ProfileCount = 0;
    }

    void WriteProfile(const ProfileResult& result)
    {
        // std::lock_guard<std::mutex> lock(m_lock);

        if (m_ProfileCount++ > 0)
            m_OutputStream << ",\n";

        std::string name = result.Name;
        std::replace(name.begin(), name.end(), '"', '\'');

        m_OutputStream << "\t\t{\n";
        m_OutputStream << "\t\t\t\"cat\":\"function\",\n";
        m_OutputStream << "\t\t\t\"dur\":" << (result.End - result.Start) << ",\n";
        m_OutputStream << "\t\t\t\"name\":\"" << name << "\",\n";
        m_OutputStream << "\t\t\t\"ph\":\"X\",\n";
        m_OutputStream << "\t\t\t\"pid\":0,\n";
        m_OutputStream << "\t\t\t\"tid\":" << result.ThreadID << ",\n";
        m_OutputStream << "\t\t\t\"ts\":" << result.Start << "\n";
        m_OutputStream << "\t\t}";

        m_OutputStream.flush();
    }

    void WriteHeader()
    {
        m_OutputStream << "{\n";
        m_OutputStream << "\t\"otherData\": {},\n";
        m_OutputStream << "\t\"traceEvents\":\n";
        m_OutputStream << "\t[\n";
        m_OutputStream.flush();
    }

    void WriteFooter()
    {
        m_OutputStream << "\n\t]\n";
        m_OutputStream << "}";
        m_OutputStream.flush();
    }

    static Profiler& Get()
    {
        static Profiler instance;
        return instance;
    }
};



class ProfileTimer
{
public:
    ProfileTimer(const char* name) : m_Name(name), m_Stopped(false)
    { m_StartTimepoint = std::chrono::steady_clock::now(); }

    ~ProfileTimer()
    {
        if (!m_Stopped)
            Stop();
    }

    void Stop()
    {
        auto endTimepoint = std::chrono::steady_clock::now();

        long long start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch().count();
        long long end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch().count();

        uint32_t threadID = 0;//std::hash<std::thread::id>{}(std::this_thread::get_id());
        Profiler::Get().WriteProfile({ m_Name, start, end, threadID });

        m_Stopped = true;
    }
private:
    const char* m_Name;
    std::chrono::time_point<std::chrono::steady_clock> m_StartTimepoint;
    bool m_Stopped;
};



#define PROFILING 0
#if PROFILING
    #define PROFILE_SCOPE(name) ProfileTimer timer##__LINE__(name)
    #define PROFILE_FUNCTION() PROFILE_SCOPE(__PRETTY_FUNCTION__)
#else
    #define PROFILE_FUNCTION()
#endif

#endif //__INCLUDE_INSTRUMENTOR_hh__