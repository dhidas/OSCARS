#ifndef GUARD_Version_h
#define GUARD_Version_h

#include <string>

// OSCARS Version information
// You must also edit by hand the version in setup.py to match these numbers

#define OSCARS_VMAJOR 1
#define OSCARS_VMINOR 36
#define OSCARS_REVISION 6

namespace OSCARS {

  std::string GetVersionString ()
  {
    // Get the version string based off of the compiler defines
    char ver[200];
    sprintf(ver, "%i.%02i.%02i", OSCARS_VMAJOR, OSCARS_VMINOR, OSCARS_REVISION);
    return std::string(ver);
  }

}




#endif
