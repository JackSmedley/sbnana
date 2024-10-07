#include "sbnana/SBNAna/Cuts/NuMIQuality.h"

namespace ana {

  //// ----------------------------------------------
  NuMIQuality::NuMIQuality()
  {
    const char* sbndata = std::getenv("SBNDATA_DIR");
    if (!sbndata) {
      std::cout << "NuMIQuality: $SBNDATA_DIR environment variable not set. Please setup "
                   "the sbndata product."
                << std::endl;
      std::abort();
    }

    filename_bad_triggered_spill = std::string(sbndata) +
                   "beamData/NuMIdata/bad_triggered_spills.csv";

    std::cout << "[NuMIQuality::NuMIQuality] filename_bad_triggered_spill = " << filename_bad_triggered_spill << std::endl;

  }

  NuMIQuality::~NuMIQuality()
  {
  }

  NuMIQuality& NuMIQuality::Instance()
  {
    static NuMIQuality numiQual;
    return numiQual;
  }


  bool NuMIQuality::IsBadTriggeredSpill(unsigned int run, unsigned int subrun, unsigned int event) const {

    std::ifstream csv_bad_triggered_spill(filename_bad_triggered_spill);
    std::string line;
    std::getline(csv_bad_triggered_spill, line); // skip first line

    while (std::getline(csv_bad_triggered_spill, line)) {
      std::istringstream iss(line);
      unsigned int entry, currentRun, currentSubrun, currentEvt;
      char comma;  // to read the commas in the CSV
      if (!(iss >> entry >> comma >> currentRun >> comma >> currentSubrun >> comma >> currentEvt)) {
        std::cerr << "Error reading line from file: " << filename_bad_triggered_spill << std::endl;
        abort();
      }

      // Check if the pair matches
      if (currentRun == run && currentSubrun == subrun && currentEvt == event) {
        printf("[BAD Trigger] (run, subrun, event) = (%d, %d, %d)\n", run, subrun, event);
        return true;  // Pair found
      }
    }
    return false;

  }


  const SpillCut kNuMINotBadTriggeredSpill( [](const caf::SRSpillProxy *sr) {

    if(sr->hdr.ismc) return true;

    const NuMIQuality& numiQual = NuMIQuality::Instance();

    if( numiQual.IsBadTriggeredSpill(sr->hdr.run, sr->hdr.subrun, sr->hdr.evt) ) return false;
    else return true;

  });

  /// \ref Good run
  const SpillCut kNuMIGoodRun( [](const caf::SRSpillProxy *sr) {

    if(sr->hdr.ismc) return true;

/*
    // used before 05/07/2024
    static std::vector<unsigned int> GoodRunList = {
8460, 8461, 8462, 8468, 8469, 8470, 8471,
8505, 8506, 8507, 8513, 8514, 8515, 8517, 8518, 8521, 8522, 8525, 8527, 8528, 8529, 8530, 8531, 8552, 8553,
9593, 9594, 9595, 9597, 9599,
9602, 9610, 9642, 9646, 9648, 9649, 9688, 9690, 9691, 9692, 9693, 9694, 9695, 9696, 9699,
9700, 9704, 9715, 9716, 9717, 9721, 9723, 9725, 9726, 9728, 9729, 9730, 9731, 9732, 9733,
9735, 9743, 9744, 9745, 9746, 9747, 9750, 9752, 9753, 9755, 9762, 9763, 9764, 9765, 9781,
9791, 9792, 9794, 9796, 9807, 9834, 9835, 9837, 9838, 9840, 9844, 9847, 9849, 9851, 9854,
9855, 9860, 9862, 9867, 9869, 9870, 9892, 9894, 9896, 9897, 9914, 9919, 9921, 9922, 9924,
9925, 9926, 9941, 9943, 9944, 9945, 9949, 9950, 9951, 9953, 9954, 9956, 9959, 9960, 9972,
9974, 9977, 9979, 9981, 10054, 10059, 10060, 10061, 10064, 10065, 10066, 10084, 10085, 10096, 10097
    };
*/
    // updated 05/07/2024
    static std::vector<unsigned int> GoodRunList = {
      8461, 8462, 8468, 8469, 8470, 8471,
      8505, 8506, 8507, 8513, 8514, 8515, 8521, 8522, 8527, 8528, 8529, 8530, 8531, 8552, 8553,
      9593, 9594, 9595, 9597, 9599,
      9602, 9610, 9642, 9646, 9648, 9649, 9688, 9690, 9691, 9692, 9693, 9694, 9695, 9696, 9699,
      9700, 9704, 9715, 9716, 9717, 9721, 9723, 9725, 9726, 9728, 9729, 9730, 9731, 9732, 9733,
        9735, 9743, 9744, 9745, 9746, 9747, 9750, 9752, 9753, 9755, 9762, 9763, 9764, 9765, 9781,
        9791, 9792, 9794, 9796,
      9807, 9834, 9835, 9837, 9838, 9840, 9844, 9847, 9849, 9851, 9854, 9855, 9860, 9862, 9867,
        9869, 9870, 9892, 9894, 9896, 9897, 
      9914, 9919, 9921, 9922, 9924, 9925, 9926, 9944, 9945, 9949, 9950, 9951, 9953, 9954, 9956,
        9959, 9960, 9972, 9974, 9977, 9979, 9981, 
      10054, 10059, 10060, 10061, 10064, 10065, 10066, 10084, 10085, 10096, 10097
    };
    unsigned int RunNum = sr->hdr.run;
    auto it = std::find(GoodRunList.begin(), GoodRunList.end(), RunNum);
    if(it != GoodRunList.end()){
      //std::cout << "RunNum is in GoodRunList" << std::endl;
      return true;
    }
    else{
      //std::cout << "RunNum is not in GoodRunList" << std::endl;
      return false;
    }

  });

}
