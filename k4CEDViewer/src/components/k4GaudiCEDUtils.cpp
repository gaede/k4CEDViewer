#include "k4GaudiCEDUtils.h"
#define endmsg std::endl

namespace k4ced {


  std::ostream& info()     {
    if( k4ced::GlobalLog::instance().level() <= LogLevel::Info ){
      std::cout <<  k4ced::GlobalLog::instance().name() << "  INFO "  ; 
      return std::cout ;
    }
    return ( null_stream ) ;
  }
  std::ostream& debug()     {
    if( k4ced::GlobalLog::instance().level() <= LogLevel::Debug ){
      std::cout <<  k4ced::GlobalLog::instance().name() << "  DEBUG "  ; 
      return std::cout ;
    }
    return ( null_stream ) ;
  }
  std::ostream& verbose()     {
    if( k4ced::GlobalLog::instance().level() <= LogLevel::Verbose ){
      std::cout <<  k4ced::GlobalLog::instance().name() << "  VERBOSE "  ; 
      return std::cout ;
    }
    return ( null_stream ) ;
  }
  std::ostream& error()     {
    if( k4ced::GlobalLog::instance().level() <= LogLevel::Error ){
      std::cout <<  k4ced::GlobalLog::instance().name() << "  ERROR "  ; 
      return std::cout ;
    }
    return ( null_stream ) ;
  }
  std::ostream& warning()     {
    if( k4ced::GlobalLog::instance().level() <= LogLevel::Warning ){
      std::cout <<  k4ced::GlobalLog::instance().name() << "  WARNING "  ; 
      return std::cout ;
    }
    return ( null_stream ) ;
  }
 
  
  CalorimeterDrawParams getCalorimeterParameters(dd4hep::Detector& theDetector, std::string name, bool selfCall ){ 

    CalorimeterDrawParams params;
    if (selfCall)
      name[1] = tolower(name[1]);
    const std::vector< dd4hep::DetElement>& calorimeters  = theDetector.detectors( "calorimeter" ) ;
    dd4hep::DetElement calo;
    for( unsigned i=0,n=calorimeters.size() ; i<n ; ++i ){
      if ((std::string)calorimeters[i].name() == name){
	calo = calorimeters[i];
	break;
      }
    }
    dd4hep::rec::LayeredCalorimeterData* caloGeo;
    try{
      caloGeo = calo.extension<dd4hep::rec::LayeredCalorimeterData>() ;
    } catch(std::runtime_error& e){
      if (!selfCall)
	return getCalorimeterParameters(theDetector, name, true);
      else{
	info() <<   "Cannot find detector data for  " << name << " -> cannot draw particles"<<endmsg;
	params.delta_z = -1;    //no spatial extension --> no drawing
	return params;
      }
    }
    if (caloGeo->layoutType == dd4hep::rec::LayeredCalorimeterData::BarrelLayout){
      params.r_inner = caloGeo->extent[0]/dd4hep::mm;
      params.delta_r = caloGeo->extent[1]/dd4hep::mm - params.r_inner;
      params.z_0 = 0.;
      params.delta_z = caloGeo->extent[3]/dd4hep::mm;
    }else{
      params.r_inner = caloGeo->extent[0]/dd4hep::mm;
      params.delta_r = caloGeo->extent[1]/dd4hep::mm - params.r_inner;
      params.z_0 = caloGeo->extent[2]/dd4hep::mm;
      params.delta_z = (caloGeo->extent[3]/dd4hep::mm - params.z_0);  //we are interested in the full length!
      //CEDGeoTube only requires half length as an argument
    }

    return params;
  }


  double calculateTrackLength(const CalorimeterDrawParams& barrel, const CalorimeterDrawParams& endcap, double x, double y, double z, double px, double py, double pz, double rel_X0){


    if (barrel.delta_z == -1 || endcap.delta_z == -1) return 0;   //the case if the parameters could not be loaded properly

    double length;
    double pt = sqrt(px*px + py*py);
    double pt_over_pz = pt/fabs(pz);

    double r = sqrt(x*x+y*y);
    if (r > barrel.r_inner || r > endcap.r_inner) return 0;
    double p = 2 * (px * x + py * y)/pt;

    double q = r*r - barrel.r_inner * barrel.r_inner;
    double distance_to_barrel_r = -p/2 + sqrt(p*p/4 - q);

    q = r*r - endcap.r_inner * endcap.r_inner;
    double distance_to_endcap_r = -p/2 + sqrt(p*p/4 - q);

    double sign_pz = (pz >= 0) ? 1. : -1.;
    double distance_to_barrel_z = barrel.delta_z - sign_pz*z;
    double distance_to_endcap_z = endcap.z_0 - sign_pz*z;

    //case 1: barrel only
    if (pt_over_pz > (distance_to_barrel_r + barrel.delta_r)/distance_to_barrel_z){
      length = rel_X0 * barrel.delta_r + distance_to_barrel_r * sqrt(1. + pow(1./pt_over_pz,2));
    }
    //case 2: barrel + endcap hit
    else if(pt_over_pz > distance_to_barrel_r/distance_to_barrel_z){
      //x == path traversed in barrel, rotation symmetry is still assumed at this point which is a valid approximation most of the times
      double X = (distance_to_barrel_z - distance_to_barrel_r/pt_over_pz)*sqrt(1+pow(pt_over_pz,2));
      //case 2a: traversed path in the barrel is larger than defined interaction path --> case 1
      if (X>=rel_X0*barrel.delta_r)
	length = rel_X0 * barrel.delta_r + distance_to_barrel_r * sqrt(1. + pow(1./pt_over_pz,2));
      //case 2b: particle is not absorbed in barrel but reaches the endcap --> length as in case 3 minus x
      else{
	length = (rel_X0 - X / barrel.delta_r) * endcap.delta_z + distance_to_endcap_z * sqrt(1. + pow(pt_over_pz,2));
      }
      //case 2c: distance from z-axis exceeds endcap extension (e.g. if particle travels through gab)
      double length_r = length * pt_over_pz/sqrt(1. + pow(pt_over_pz,2));
      if (length_r > (distance_to_endcap_r + endcap.delta_r)){
	length = length * (distance_to_endcap_r + endcap.delta_r)/length_r;
      }
    }
    //case 3: endcap only
    else if(pt_over_pz > distance_to_endcap_r/distance_to_endcap_z){
      length = rel_X0 * endcap.delta_z + distance_to_endcap_z * sqrt(1. + pow(pt_over_pz,2));
    }
    //case 4: part of endcap hit
    else if(pt_over_pz > distance_to_endcap_r/(distance_to_endcap_z + endcap.delta_z)){
      double X = (distance_to_endcap_r/pt_over_pz - distance_to_endcap_z)*sqrt(1+pow(pt_over_pz,2));
      //case 4a: traversed path in endcap is larger than defined interaction path --> case 3
      if (X>=rel_X0*endcap.delta_z)
	length = rel_X0 * endcap.delta_z + distance_to_endcap_z * sqrt(1. + pow(pt_over_pz,2));
      //case 4b: particle is not fully absorbed the endcap --> draw up to the yoke
      else
	length = (distance_to_endcap_z + endcap.delta_z) * sqrt(1. + pow(pt_over_pz,2));
    }
    //case 5: neither the endcap nor a barrel is hit
    else{
      length = (distance_to_endcap_z + endcap.delta_z) * sqrt(1. + pow(pt_over_pz,2));
    }
    return fabs(length);
  }


//get the outer extents of the tracker
  std::array<double, 2> getTrackerExtent(dd4hep::Detector& theDetector){
    
    double extent0 =theDetector.constant<double>("tracker_region_rmax")/dd4hep::mm;
    double extent1 = theDetector.constant<double>("tracker_region_zmax")/dd4hep::mm;
    return { extent0, extent1 } ;
  }

  std::array<double, 2> getYokeExtent(dd4hep::Detector& theDetector) {
    double extent0(0.),extent1(0.) ;
    const std::vector< dd4hep::DetElement>& calorimeters     = theDetector.detectors( "calorimeter" ) ;
    dd4hep::DetElement yoke;
    for( unsigned i=0,n=calorimeters.size() ; i<n ; ++i ){
      std::string detName = calorimeters[i].name();
      bool isYokeBarrel = (detName == "YokeBarrel") ;
      bool isYokeEndcap = (detName == "YokeEndcap") ;
      if (!isYokeBarrel && !isYokeEndcap)
	continue;
      yoke = calorimeters[i];
      dd4hep::rec::LayeredCalorimeterData* yokeGeo;
      try{
	yokeGeo = yoke.extension<dd4hep::rec::LayeredCalorimeterData>() ;
      } catch(std::runtime_error& e){
	debug() <<  " cannot get detector data for  " << detName << " cannot draw ;-( "<< endmsg ;
	extent0 = extent1 = 0;
	return { extent0, extent1 } ;
      }
      if (isYokeBarrel)
	extent0 = yokeGeo->extent[1]/dd4hep::mm;
      if (isYokeEndcap)
	extent1 = yokeGeo->extent[3]/dd4hep::mm;
    }
    return  { extent0, extent1 } ;
  }


} // end namespace
