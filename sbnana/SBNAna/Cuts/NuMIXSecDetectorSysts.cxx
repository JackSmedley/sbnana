#include "sbnana/SBNAna/Cuts/NuMIXSecDetectorSysts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <iostream>

namespace ana {

  NuMIXSecDetectorSysts::NuMIXSecDetectorSysts(DetSystType detsyst_type, const std::string& name, const std::string& latexName):
    ISyst(name, latexName),
    kDetSystType(detsyst_type)
  {

  }
  void NuMIXSecDetectorSysts::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {

    int RecoMuonIdx = kRecoMuonIdx(sr);
    if(RecoMuonIdx<0) return;
    int RecoProtonIdx = kRecoProtonIdx(sr);
    int ScndProtonIdx = kScndProtonIdx(sr);
    if(RecoProtonIdx<0) return;

    //Take subleading proton momentum when available
    double RecoProtonP = ( ScndProtonIdx >=0 ? kScndProtonP(sr) : kRecoProtonP(sr) );

    if(kDetSystType==kFrontIndPlaneGain){

      // pre-reprocessing
      //weight *= 1. + sigma * ( RecoProtonP<0.5 ? 0.10 : 0.05 );

      // post-reprocessing; flat 3.8%
      weight *= 1. + sigma * (0.038) ;

    }
    else if(kDetSystType==kFrontIndPlaneNoise){

      // negative is mirrored
      double this_sigma = abs(sigma);
      
      // pre-reprocessing
      //weight *= 1. + this_sigma * (-0.10);

      // post-reprocessing; flat 10.5%
      weight *= 1. + this_sigma * (-0.105) ;

    }
    else if(kDetSystType==kFrontIndPlaneSignalShape){
      // negative is mirrored
      double this_sigma = abs(sigma);

      // pre-reprocessing
      weight *= 1. + this_sigma * (-0.10);
    }
    else if(kDetSystType==kFrontIndPlaneSignalShapeFitted){

      // negative is mirrored
      double this_sigma = abs(sigma);

      double this_stepfunc = GetSmoothStepFunction(
        RecoProtonP,
        -6.58186996e+01,
        6.76128082e-01,
        -1.99887838e+03,
        1.99984206e+03
      );
      double this_onesig = abs(1.-this_stepfunc) * (-1.);
      weight *= 1. + this_sigma * this_onesig;
    }
    else if(kDetSystType==kMiddleIndPlaneTransparency){

      // pre-reprocessing
      //weight *= 1. + sigma * 0.05;

      // post-reprocessing; flat 7.9%
      weight *= 1. + sigma * 0.079;

    }
    else if(kDetSystType==kCaloGain){

      // post-reprocessing; flat 8.0%
      weight *= 1. + sigma * 0.03;

    }
    else if(kDetSystType==kElectronLifetime){

      weight *= 1. + sigma * 0.06;

    }

    else{

    }


  }

  void NuMIXSecDetectorSysts::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

  }

  double NuMIXSecDetectorSysts::GetSmoothStepFunction(double x, double x_min, double x_max, double offset, double scale) const {

    double X = (x - x_min) / (x_max - x_min);
    if(X<0) X = 0.;
    if(X>1.0) X = 1.0;

    double out = (3.-2.*X)*(X*X);
    return offset + scale * out;

  }

  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneGainSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneGain,
    "NuMIXSecFrontIndPlaneGainSyst",
    "Front ind. plane gain #pm10%"
  );
  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneNoiseSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneNoise,
    "NuMIXSecFrontIndPlaneNoiseSyst",
    "Front ind. plane noise +10%"
  );
  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneSignalShapeSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneSignalShape,
    "NuMIXSecFrontIndPlaneSignalShapeSyst",
    "Front ind. plane signal shape"
  );
  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneSignalShapeFittedSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneSignalShapeFitted,
    "NuMIXSecFrontIndPlaneSignalShapeFittedSyst",
    "Front ind. plane signal shape"
  );
  const NuMIXSecDetectorSysts kNuMIXSecMiddleIndPlaneTransparencySyst(
    NuMIXSecDetectorSysts::kMiddleIndPlaneTransparency,
    "NuMIXSecMiddleIndPlaneTransparencySyst",
    "Middle ind. plane transparency"
  );
  const NuMIXSecDetectorSysts kNuMIXSeCaloGainSyst(
    NuMIXSecDetectorSysts::kMiddleIndPlaneTransparency,
    "NuMIXSecCaloGainSyst",
    "Calorimetry gain"
  );
  const NuMIXSecDetectorSysts kNuMIXSeElectronLifetimeSyst(
    NuMIXSecDetectorSysts::kMiddleIndPlaneTransparency,
    "NuMIXSecElectronLifetimeSyst",
    "Electron lifetime 9.2ms"
  );

} // end namespace ana
