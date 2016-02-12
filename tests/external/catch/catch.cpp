#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

// On Windows, disable the "Application error" dialog box, because it
// requires an human intervention, and there is no one on Appveyor.
//
// On UNIX, does nothing
void silent_crash_handlers();

int main(int argc, char* argv[]) {
    silent_crash_handlers();
    return Catch::Session().run(argc, argv);
}

#if (defined(WIN32) || defined(WIN64))
    #include <windows.h>
    void silent_crash_handlers() {
        SetErrorMode(GetErrorMode() | SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
    }
#else
    void silent_crash_handlers() {}
#endif
