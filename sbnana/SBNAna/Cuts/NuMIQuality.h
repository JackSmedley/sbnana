#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


namespace ana
{

  class NuMIQuality
  {
  public:
    NuMIQuality();
    ~NuMIQuality();
    static NuMIQuality& Instance();

    bool IsBadTriggeredSpill(unsigned int run, unsigned int subrun, unsigned int event) const;

  protected:
    std::string filename_bad_triggered_spill;
  };

  extern const SpillCut kNuMINotBadTriggeredSpill;

  extern const SpillCut kNuMIGoodRun;

}
