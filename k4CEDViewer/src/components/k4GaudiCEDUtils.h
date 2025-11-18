/***********************************************************************************************
 Some utility functions ported from DDMarlinCED for key4hep and Gaudi

 @author F.Gaede, DESY
 @date   Nov 2025 
***********************************************************************************************/

#ifndef k4GaudiCEDUtils_h
#define k4GaudiCEDUtils_h 1

// --- DD4hep ---
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"
#include "DD4hep/DD4hepUnits.h" 
#include "DDRec/DetectorData.h"


// workaround for Gaudi logging
#define endmsg std::endl

inline std::ostream& info()   { return std::cout ; }
inline std::ostream& debug()  { return std::cout ; }
inline std::ostream& warning(){ return std::cout ; }
inline std::ostream& error()  { return std::cerr ; }


namespace k4ced {

//structure for calculating the track length of given particles
  struct CalorimeterDrawParams {
    double r_inner, delta_r, z_0, delta_z;
  };
  
  
/**Helper struct for drawing collections*/
  struct DrawParameters{
    DrawParameters(const std::string& colName, int size, int marker, int layer ) :
      ColName( colName ),
      Size( size ),
      Marker( marker ),
      Layer( layer ) {
    }
    std::string ColName ;
    int Size ;
    int Marker ;
    int Layer ;
  };


  CalorimeterDrawParams getCalorimeterParameters(dd4hep::Detector& theDetector, std::string name, bool selfCall=false );
 
  double calculateTrackLength(std::string barrelName, std::string endcapName, dd4hep::Detector& theDetector, double x, double y, double z, double px, double py, double pz, double rel_X0= 0.5) ;


//get the outer extents of the tracker
  double* getTrackerExtent(dd4hep::Detector& theDetector) ;

//get the outer extents of the yoke
  double* getYokeExtent(dd4hep::Detector& theDetector) ;

} /// end namespace

#endif

