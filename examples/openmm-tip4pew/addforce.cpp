#include "OpenMM.h"

// System::addForce takes ownership. The Python SWIG wrapper sometimes
// refuses a PlumedForce built in a second module; the pointer is still
// an OpenMM::Force.
extern "C" void dseams_openmm_add_force(OpenMM::System *system, OpenMM::Force *force) {
  system->addForce(force);
}
