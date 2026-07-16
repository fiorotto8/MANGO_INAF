// Compatibility entry point. The maintained reducer is RecoTrack.C.
#include "RecoTrack.C"

void RecoTrack_faster(std::string filename,
                      std::string particle = "e-",
                      bool primaryOnly = true)
{
    RecoTrack(filename, particle, primaryOnly);
}
