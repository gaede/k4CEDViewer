#include "k4GaudiCED.h"

#include "k4GaudiCEDUtils.h"

//#include <LCGeometryTypes.h>
#include "ced_cli.h"
 
#include <signal.h>


#include <unistd.h>//

//SJA:FIXED:added to make gcc4.3 compliant
#include <cstdlib>
#include <cstdio>

//hauke
#include <ctime>
#include <time.h>
#include <termios.h>
#include <poll.h>

//for kbhit
#include <sys/select.h>
#include <termios.h>



//for detector drawing (Thorben Quast)
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h" 
#include "DDRec/DetectorData.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"
#include "TColor.h"
#include "TROOT.h"

#define endmsg std::endl

using namespace dd4hep;

namespace k4ced {
  
  k4GaudiCED* k4GaudiCED::_me = 0;

  
  std::vector<std::string> k4GaudiCED::_descs(CED_MAX_LAYER, ""); //layer descriptions


//--------------------------------------------------------------------------------------------------------


  int PickingHandler::kbhit(void) {
    //http://stackoverflow.com/questions/448944/c-non-blocking-keyboard-input#448982
    struct timeval tv = { 0L, 0L };
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(0, &fds);
    return select(1, &fds, NULL, NULL, &tv);
  }


  void k4GaudiCED::add_layer_description(const std::string &desc, int layerID){
    std::string tmp;
    if(layerID > CED_MAX_LAYER || layerID < 0){return;}
    if( _descs.at(layerID).find(desc.c_str()) == std::string::npos){
      tmp=_descs.at(layerID);
      if(! tmp.empty()){
	tmp.append(", ");
      }
      tmp.append(desc);
      _descs.at(layerID)=tmp;
    }else{
    }
  }

  void k4GaudiCED::set_layer_description(const std::string &desc, int layerID){
    if(layerID > CED_MAX_LAYER || layerID < 0){return;}
    _descs.at(layerID)=desc;
    
  }

  void k4GaudiCED::write_layer_description(void){
    //std::cout<<"LAYER: write all layer in ced" << std::endl;
    unsigned int i;
    //for(i=0;i<25;i++){
    for(i=0; i<_descs.size(); i++){
      ced_describe_layer(_descs.at(i).c_str(), i);
    }
  }

  k4GaudiCED* k4GaudiCED::instance() {
    if( _me == 0 )
      _me = new k4GaudiCED ;
    return _me ;
  }


  void k4GaudiCED::init(const void* proc ) {
    
    if( instance()->_first == 0 ){
        
      instance()->_first = proc ;
        
      char *port, *host;
      port = getenv("CED_PORT");
      host = getenv("CED_HOST");
      if((port == NULL || port[0] == 0) && (host == NULL || host[0] == 0)){
	ced_client_init("localhost",7286);
      }else if(port == NULL || port[0] == 0){
	info()<< "Use user defined host " << host << endmsg;
	ced_client_init(host,7286);
      }else if(host == NULL|| host[0] == 0){
	info()<< "Use user defined port " << port << endmsg;
	ced_client_init("localhost",atoi(port));
      }else{
	info()<< "Use user defined host " << host << ", port " <<  port << endmsg;
	ced_client_init(host,atoi(port));
      }
        
      ced_register_elements();
    }
    
    instance()->_last = proc ;
  }


  void k4GaudiCED::newEvent(const void* proc ) {
    if( proc == instance()->_first ) {
      ced_new_event();

      PickingHandler::instance().clear() ;
    }
  }

//hauke hoelbe modify 08.02.2010
  void k4GaudiCED::draw(const void* proc , int waitForKeyboard ) {
    int i=0;
    
    
    if( proc == instance()->_last ) {
      //    ced_draw_event();
      k4GaudiCED::write_layer_description();
      //ced_picking_text("test1 test2 test3");
        
      ced_send_event();
      if ( waitForKeyboard == 1 ) {
	info() << "Double click for picking. Press <ENTER> for the next event." << endmsg;
	//test:
            
	signal(SIGWINCH,SIG_IGN);
            
	while(!PickingHandler::kbhit()){
	  //            while(!poll(pfd,1,0)){
                
	  usleep(100000); //micro seconds
                
	  int id = ced_selected_id_noblock();

	  if(id>=0) {

	    debug() << "debug: got id: " << id <<endmsg;

	    if(id == 0){

	      warning() << "Picking nothing, or an object with ID 0!" << endmsg;

	    }else{

	      PickingHandler::instance().printObject(id);

	      ced_picking_text("test1 test2 test3",i++);
	      ced_send_event();
	    }
	  }
	}
            
	signal(SIGWINCH,SIG_IGN);
	int c = getchar();
	if(c=='q'||c=='Q'||c==3){ //quit if the user pressed q or strg+c (3 = strg+c)
	  exit(0);
	}
	info() << "--------- END ---------------\n";
      }
    }
  }

/**
 * Improved drawHelix() method. Draws straight lines as well.
 */
//SM-H: Added id to drawHelix (default zero), which allows for implementation of picking
  void k4GaudiCED::drawHelix(float b, float charge, float x, float y, float z,
			     float px, float py, float pz, int marker, int size, unsigned int col,
			     float rmin, float rmax, float zmax, unsigned int id)  {
    // FIXME : check for zmin as well, i.e. cylindrical coordinates
    
    double cFactor = 2.9979251e-4;
    const double high_pt = 100.0;//Transverse momentum high enough for the particle not to curve noticeably
    double pt = sqrt(px*px + py*py);
    
    // FIXME: use a parameter for this cut or better this should be a function of the B field, charge and momentum 2006/07/04 OW
    
    // SD: FIXME: Adaptive step-number (or get rid of it!) and adaptive draw step!
    if ( (pt >= 0.01) && (pt <= high_pt && charge!=0) ) {
      double r =  pt / ( cFactor * b * std::abs( charge )  ) ;
      double sign =  charge > 0 ? 1 : -1 ;
        
      sign = - sign  ; // FIXME: need to check the convention - but this works !?
        
      double phi = std::atan2( py , px ) + ( 2. + sign ) * M_PI / 2. ;
      //center of helix
      double cx = x - ( sign * py * r / pt ) ;
      double cy = y + ( sign * px * r / pt ) ;
      double cz = z ;
        
      double x1 =  x ;
      double y1 =  y ;
      double z1 =  z ;
      double step = 0.05;  // initial 0.05
        
      // FIX ME: do the adaptive step number...
        
      // cheap adaptive algorithms
      if (px>1 || py >1 || px <-1 || py <-1 ){
	step = 0.005;
	if (px>5 || py >5 || px <-5 || py <-5){
	  step = 0.001;
	}
      }
        
      int nSteps = int(100/step); //hauke
        
        
      int count_lines=0;
      for (int j = 0; j < nSteps ; j++) {
            
	double alpha0 = step*j ;
            
	double x2 = cx + r * cos( phi + sign * alpha0 ) ;
	double y2 = cy + r * sin( phi + sign * alpha0 ) ;
	double z2 = cz + r * alpha0 * pz / pt ;
            
	double r_current  = sqrt(x2*x2 + y2*y2); // hypot( x2, y2 )
            
	/*
	 *  interpolation and loop break
	 */
	if( std::abs(z2) > zmax || r_current > rmax  ) {
                
	  double alpha = step*(j+0.5);
                
	  x2 = cx + r * cos( phi + sign * alpha ) ;
	  y2 = cy + r * sin( phi + sign * alpha ) ;
	  z2 = cz + r * alpha * pz / pt ;
	  break ;
	}
            
	if( r_current >= (rmin+step)) {
	  count_lines++;
	  ced_line_ID( x1, y1, z1, x2, y2, z2 , marker , size, col, id);
	}
	x1 = x2;
	y1 = y2;
	z1 = z2;
            
      }
    }
    //For high momentum tracks, just draw straight line
    else if (pt > high_pt) {
      debug() << "pt = " << pt << endmsg;
      float absP =sqrt(px*px + py*py + pz*pz);
      float k = 0.0;
      float kr = 0.0;
      float kz = 0.0;
      float summand = 0.0;
      float radicant = 0.0;
        
      // find intersection with rmax
      summand = (-1)*( absP*(px*x + py*y)/(pow(px,2) + pow(py,2)) );
      radicant = summand*summand - ( (pow(absP,2)*(pow(x,2)+pow(y,2)-pow(rmax,2)))/(pow(px,2) + pow(py,2)) );
        
      if (radicant < 0) {
	error() << "Error in 'k4GaudiCED::drawHelix()': Startpoint beyond (rmax,zmax)" << endmsg;
	return;
      }
        
      kr = summand + sqrt(radicant);
      kz = ((zmax-z)*absP)/pz;
        
      // this has been improved
        
      if (z + (kr*pz)/absP > zmax || z + (kr*pz)/absP < -zmax){
	k = kz;
      }
      else k = kr;
        
      if (k < 0.0 ) {
	debug() << "k4GaudiCED::drawHelix(): negative intersection parameter - will revert sign ... "
		<< endmsg;
	//fg: k cannot be negativ ( particle is moving along its 3-momentum ....)
	k = -k ;
      }
        
      float xEnd = x + (k*px)/absP;
      float yEnd = y + (k*py)/absP;
      float zEnd = z + (k*pz)/absP;
        
      if (rmin != 0){
	debug() << "FIX ME: Inner cylinder not taken into account!" << endmsg;
	return;
      }
        
        
      debug() << "k4GaudiCED::drawHelix()' - pt : " << pt << " |p| = " << absP
	      << ", x " << x
	      << ", y " << y
	      << ", z " << z
	      << ", px " << px
	      << ", py " << py
	      << ", pz " << pz
	      << ", xEnd " << xEnd
	      << ", yEnd " << yEnd
	      << ", zEnd " << zEnd
	      << endmsg ;
        
        
      ced_line_ID(x, y, z, xEnd, yEnd, zEnd , marker , size, col, id);
        
    }
    else {
      debug() << "Low momentum particle given point instead of helix" << endmsg;
      const double delta = 0.0001;
      ced_line_ID(x, y, z, x+delta, y+delta, z+delta, marker , size, col, id);
    }
  }


  void k4GaudiCED::drawDD4hepDetector( dd4hep::Detector& lcdd, bool _surfaces, StringVec _detailled){
    typedef std::vector< dd4hep::DetElement> DetVec ;
    // get DetElements for the main sub detectors from dd4hep 
    const DetVec& trackers     = lcdd.detectors( "tracker" ) ;
    const DetVec& calorimeters = lcdd.detectors( "calorimeter" ) ;
    const DetVec& passiveDets  = lcdd.detectors( "passive" ) ;


    //allocate reference to the surface manager
    // some models might not have a SurfaceManager extension:
    dd4hep::rec::SurfaceManager* sM = 0 ;
    try{  sM = lcdd.extension<dd4hep::rec::SurfaceManager>();
    } catch(const std::runtime_error& ) {
      lcdd.apply( "InstallSurfaceManager",0,0);
      sM = lcdd.extension<dd4hep::rec::SurfaceManager>();
    }
    const dd4hep::rec::SurfaceManager& surfMan = *sM;

    //some temporary parameters for visualization
    unsigned color; bool visible;
    //temporary objects
    DetElement det; std::string detName;

    std::vector<CEDGeoTube> gTV ; 
    int detLayer = NUMBER_DATA_LAYER; 


    //--- loop over all calorimeters
    for( unsigned i=0,n=calorimeters.size() ; i<n ; ++i ){
      det = calorimeters[i] ;
      detName = det.name() ;
      dd4hep::rec::LayeredCalorimeterData* calo = nullptr;
      debug( ) << " ......processing " << detName << endmsg;  
      //try to get the appropriate extension
      
      try{ 
	calo = det.extension<dd4hep::rec::LayeredCalorimeterData>();
      } catch(const std::runtime_error& e){
	debug( ) <<  detName 
		 << " has no extension of type LayeredCalorimeterData. "
		 <<   endmsg;  
      }

      //get the visAttributes of the detElement's volume or (if not existing) some default values
      getVisAttributes(det, color, visible);
      int layer = detLayer++;
      bool isDrawn = _surfaces && DrawSurfaces(surfMan, detName, color, layer);

      //draw if the object exists and its visibility is set true
      if(visible && !isDrawn) {                    
	if (calo != 0){
	  //get the required parameter set for drawing a CEDGeoTube
	  CEDGeoTubeParams params;    
	  params = CalorimeterParameterConversion(calo);
	  //consistently allocated this helper for linking the element correctly to the GUI control via ID
	  gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
	  if (!params.isBarrel){
	    //place the second one symmetric to the z-axis. An additional shift by the width of the layer is needed since the appropriate argument represents the left handed start of the geometry
	    gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  - (params.z0+2.0*params.delta_z), color , layer  ,1,1 ) ) ; 
	  }
	  isDrawn = true;
	}
	if (!calo)
	  isDrawn = DrawSurfaces(surfMan, detName, color, layer);
      }
      if(isDrawn){
	k4GaudiCED::set_layer_description( detName , layer );
	debug( ) << detName << " has successfully been processed."<< endmsg;
      }
      std::cout<<std::endl;
    }
    //--- draw trackers 
    //mostly repetition of the calorimeter 
    for( unsigned i=0,n=trackers.size() ; i<n ; ++i ){
      det = trackers[i] ;
      detName = trackers[i].name() ;
      dd4hep::rec::ZPlanarData* trkPlanar = nullptr;
      dd4hep::rec::ZDiskPetalsData* trkDisk = nullptr;
      dd4hep::rec::FixedPadSizeTPCData* trkTPC = nullptr;
      debug( ) << " ......processing" <<  detName << endmsg; 
    
      try{ 
	trkPlanar = det.extension<dd4hep::rec::ZPlanarData>(false);
      } catch(const std::runtime_error&){
	try{
	  trkDisk = det.extension<dd4hep::rec::ZDiskPetalsData>(false);
	}catch(const std::runtime_error&){
	  try{
	    trkTPC = det.extension<dd4hep::rec::FixedPadSizeTPCData>(false);
	  } catch(const std::runtime_error&){
	    debug( ) <<  detName 
		     << " has no extension of type ZPlanarData/ZDiskPetalsData. "
		     <<   endmsg;           
	  }
	}
      }

      int layer = detLayer++;  
      getVisAttributes(det, color, visible);
      bool isDrawn = _surfaces && DrawSurfaces(surfMan, detName, color, layer);
      //get the visAttributes of the detElement's volume or (if not existing) some default values
      if (!isDrawn && visible){
	//the following if statements are exclusive, i.e. only one may apply
	if(trkPlanar){
	  for (std::vector<dd4hep::rec::ZPlanarData::LayerLayout>::iterator thisLayer = trkPlanar->layers.begin(); thisLayer != trkPlanar->layers.end(); thisLayer++){
	    LayerGeometry Geo;
	    Geo = TrackerLayerParameterConversion(thisLayer);
	    if (detailledDrawing(_detailled, detName)){
	      for( unsigned stave_i=0; stave_i<Geo.staves.size() ; ++stave_i ){
		ced_geobox_r_ID( Geo.staves[stave_i].sizes, Geo.staves[stave_i].center, Geo.staves[stave_i].rotate, color, layer,0);
		ced_geobox_r_solid( Geo.staves[stave_i].sizes, Geo.staves[stave_i].center, Geo.staves[stave_i].rotate, color, layer);
	      }
	    }
	    else{ 
	      gTV.push_back( CEDGeoTube( Geo.tube.Rmax, Geo.tube.Rmin, Geo.tube.inner_symmetry, Geo.tube.outer_symmetry, Geo.tube.phi0, Geo.tube.delta_phi, Geo.tube.delta_z,  Geo.tube.z0, color , layer  ,1,1 ) );
	    }
	  }
	  isDrawn = true;
	}

	if(trkDisk){ 
	  if (detailledDrawing(_detailled, detName)){
	    debug( )<<detName<<": Not drawn for now (appropriate geometry does not exist)"<<endmsg;
	  }
	  else{
	    for (std::vector<dd4hep::rec::ZDiskPetalsData::LayerLayout>::iterator thisLayer = trkDisk->layers.begin(); thisLayer != trkDisk->layers.end(); thisLayer++){
	      CEDGeoTubeParams params;    
	      params = PetalParameterConversion(thisLayer);
	      gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
	      //place the second one symmetric to the z-axis. An additional shift by the width of the layer is needed since the appropriate argument represents the left handed start of the geometry
	      gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  -(params.z0+2*params.delta_z), color , layer  ,1,1 ) ) ; 
	      k4GaudiCED::set_layer_description( detName , layer );
	    }
	  }
	  isDrawn = true;
	}
      
	if(trkTPC) {
	  CEDGeoTubeParams params;    
	  params = TPCParameterConversion(trkTPC);
	  gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
	  isDrawn = true;
	}

	if(!trkTPC && !trkDisk && !trkPlanar)//if no simplified geometry is given, try surfaces
	  isDrawn = DrawSurfaces(surfMan, detName, color, layer);
      
      }
      if(isDrawn){
	k4GaudiCED::set_layer_description( detName , layer );
	debug( ) << detName << " has successfully been processed."<< endmsg;
      }
      std::cout<<std::endl;
    }


    //--- draw passive 
    for( unsigned i=0,n=passiveDets.size() ; i<n ; ++i ){
      det = passiveDets[i] ;
      detName = passiveDets[i].name() ;
      dd4hep::rec::ConicalSupportData* passiveConical = nullptr;
      dd4hep::rec::LayeredCalorimeterData* passiveCalo = nullptr;
      debug( ) << " ......processing " <<  detName << endmsg;  
      try{ 
        passiveConical = det.extension<dd4hep::rec::ConicalSupportData>(false);
      } catch(const std::runtime_error&){
	try{
	  passiveCalo = det.extension<dd4hep::rec::LayeredCalorimeterData>(false);
	} catch(const std::runtime_error&){
	  debug( ) <<  detName 
		  << " has no extension of type ConicalSupportData/LayeredCalorimeterData. "
		  <<   endmsg;  
	}
      }
      //get the visAttributes of the detElement's volume or (if not existing) some default values
      getVisAttributes(det, color, visible);
      int layer = detLayer++;
      bool isDrawn = _surfaces && DrawSurfaces(surfMan, detName, color, layer);
      if (!isDrawn && visible){
	if (passiveConical){
	  debug( )<<detName<<" is not drawn for now (not needed)."<<endmsg;
	  for (std::vector<dd4hep::rec::ConicalSupportData::Section>::iterator thisSection = passiveConical->sections.begin(); thisSection != passiveConical->sections.end(); thisSection++){
	  }
	  isDrawn = true;
	}
	if (passiveCalo){
	  for (std::vector<dd4hep::rec::LayeredCalorimeterData::Layer>::iterator thisLayer = passiveCalo->layers.begin(); thisLayer != passiveCalo->layers.end(); thisLayer++){
	    CEDGeoTubeParams params;
	    params = CalorimeterLayerParameterConversion(thisLayer);
	    gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry,  params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
	  }
	  isDrawn = true;
	}
	if (!passiveConical && !passiveCalo)
	  isDrawn = DrawSurfaces(surfMan, detName, color, layer);
      }
      if (isDrawn){
	k4GaudiCED::set_layer_description( detName , layer );
	debug( ) << detName << " has successfully been processed."<< endmsg;
      }
      std::cout<<std::endl;
    }
    // ========================================================================
    //Draw the tubes:
    ced_geotubes( gTV.size() ,  (CED_GeoTube*) &gTV[0] );

    // ========================================================================
  
    k4GaudiCED::write_layer_description();
  }

  void DDdraw_helix( float b, float charge, float x, float y, float z,
		     float px, float py, float pz,
		     int marker, int size, unsigned int col,
		     float rmin, float rmax, float zmax, unsigned int id){
    
    
    k4GaudiCED::drawHelix( b,  charge,  x,  y,  z, px,  py,  pz,  marker,  size, col, rmin,  rmax,  zmax,  id) ;
    
  }



/***detector draw helpers***/

  bool detailledDrawing(StringVec _detailled, std::string detName){
    if(_detailled.size() == 0) 
      return false;
    unsigned index = 0;
    while(index < _detailled.size()){
      if (detName.compare(_detailled[index])==0)
	return true;
      index++;
    }
    return false;
  }

//converts the parameters in LayeredCalorimeterData given by the appropriate drivers
//into those required by the CEDGeoTube
  CEDGeoTubeParams CalorimeterParameterConversion (dd4hep::rec::LayeredCalorimeterData *calo){
    //get all the information from the lcdd class
    double rMin = calo->extent[0]/dd4hep::mm ;
    double rMax = calo->extent[1]/dd4hep::mm ;
    //attention! the meaning of these two variables depends on their context 
    //(see e.g. ECalEndcap_o2_v01_geo.cpp lines 78/79 vs. ECalBarrel_o1_v01.cpp lines 74/75)
    double zMin = calo->extent[2]/dd4hep::mm ; 
    double zMax = calo->extent[3]/dd4hep::mm ;
  
    double inner_symmetry = calo->inner_symmetry;
    double outer_symmetry = calo->outer_symmetry;
    //Note that the phi0's are given in degrees rather than rad in this implementation
    double inner_phi0 = calo->inner_phi0*180/M_PI;
    double outer_phi0 = calo->outer_phi0*180/M_PI;
    //new convention to prevent weird overlaps
    int NCircle = 36;
    //correct for the case of circle approximation, the given is conventional
    if (outer_symmetry < 1)  outer_symmetry = NCircle;
    if (inner_symmetry < 1)  inner_symmetry = NCircle;

    CEDGeoTubeParams returnParams;
    //given: distance middle point - center of edge section
    //required: distance middle point - edge intersection (corner) to prevent weird overlaps
    returnParams.Rmax = rMax/cos(M_PI/outer_symmetry); 
    returnParams.Rmin = rMin/cos(M_PI/inner_symmetry);
  
    //nothing to convert
    returnParams.inner_symmetry = inner_symmetry;
    returnParams.outer_symmetry = outer_symmetry;
  
    //by default, CED draws tubes with the corner facing downward
    //what we want: phi0 = 0 <=> straight edge parallel to x-y plane
    //therefore, the tube must be rotated by 360 - (90 + 360./(2*number of sides)) degrees since 
    //respecting the rotation symmetry in 360/n, the minimal phi0 is calculated to prevent interference with phi cuts in the CEDViewer
    //phi0 == 0 implies a normal vector of the first cell parallel to the +x-axis 
    //but CED starts drawing symmetrically from the +y-axis
    returnParams.phi0 = outer_phi0 + 270. - 180./outer_symmetry;
    returnParams.phi0 = returnParams.phi0 - (360./outer_symmetry)*(int (returnParams.phi0/(360./outer_symmetry))); 
    //in ILD and CLIC, both inner and outer shapes agree. For a generic solution, another parameter in LayeredCalorimeterData (e.g. phi_inner) should be introduced
    //i.e. delta_phi that is the rotational angle of the inner angle with respect to the outer layer
    returnParams.delta_phi = -returnParams.phi0 + inner_phi0 + 270. - 180./inner_symmetry;
    returnParams.delta_phi = returnParams.delta_phi - (360./inner_symmetry)*(int (returnParams.delta_phi/(360./inner_symmetry)));

  

    //endcaps and barrels take the same parameters in CED but the interpretation of the z-coordinates differs:
    returnParams.isBarrel = calo->layoutType == dd4hep::rec::LayeredCalorimeterData::BarrelLayout;
  
    //type specific conversions
    if (returnParams.isBarrel){
      //barrels are drawn centered at z=0, i.e. they start at -zMax (z0) and their half length (delta_z) is zMax
      returnParams.delta_z = zMax;
      returnParams.z0 = -zMax;
    }
    else {
      //In the drivers, zMin is given to be the start of one disk. (Symmetry at the z-axis is assumed during placement.)
      //zMax is passed as zMin + thickness such that half the distance is correctly obtained by the line below.
      returnParams.delta_z = 0.5*(zMax - zMin);
      returnParams.z0 = zMin;
    }
    return returnParams;
  }
//converts the parameters in ZDiskPetalsData given by the appropriate drivers
//into those required by the CEDGeoTube
  CEDGeoTubeParams PetalParameterConversion (std::vector<dd4hep::rec::ZDiskPetalsData::LayerLayout>::iterator thisLayer){

    double phi0 = thisLayer->phi0*180/M_PI;
    double distanceSensitive = thisLayer->distanceSensitive/dd4hep::mm;
    //Supposedly, this is the actual width of the ring. The expansion along the z-axis is hard coded to 0.3.
    double lengthSensitive = thisLayer->lengthSensitive/dd4hep::mm;
    double zPosition = thisLayer->zPosition/dd4hep::mm;
    double thicknessSupport = thisLayer->thicknessSupport/dd4hep::mm;
    int petalNumber = thisLayer->petalNumber;

    CEDGeoTubeParams returnParams;
    /*   
    //Old implementation with GEAR
    returnParams.Rmax = distanceSensitive;
    returnParams.Rmin = distanceSensitive-1.0*lengthSensitive;
    */
    returnParams.Rmax = distanceSensitive+1.0*lengthSensitive;
    returnParams.Rmin = distanceSensitive;

    //(see comment in line 151)
    returnParams.Rmax = returnParams.Rmax/cos(M_PI/petalNumber);
    returnParams.Rmin = returnParams.Rmin/cos(M_PI/petalNumber);
    //Number of edges = number of petals
    returnParams.inner_symmetry = petalNumber;
    returnParams.outer_symmetry = petalNumber;
    //(see comment in line 159)
    returnParams.phi0 = phi0 + 270. - 180./petalNumber - (360./petalNumber)*(int ((270.-180./petalNumber)/(360./petalNumber)));
    returnParams.delta_phi = 0.0;
  
    //thicknessSupport is negligibly small
    returnParams.delta_z = thicknessSupport; 
    //Again: z0 is the left handed starting point for drawing
    returnParams.z0 = zPosition-returnParams.delta_z;

    returnParams.isBarrel = false;

    return returnParams;
  }

//converts the parameters from a LayeredCalorimeterData layer given by the appropriate drivers
//into those required by the CEDGeoTube
  CEDGeoTubeParams CalorimeterLayerParameterConversion(std::vector<dd4hep::rec::LayeredCalorimeterData::Layer>::iterator thisLayer){
    double distance = thisLayer->distance/dd4hep::mm;
    double thickness = thisLayer->inner_thickness/dd4hep::mm + thisLayer->outer_thickness/dd4hep::mm ;
    double cellSize0 = thisLayer->cellSize0/dd4hep::mm;
    double cellSize1 = thisLayer->cellSize1/dd4hep::mm;
  
    int NCircle = 36;   //hard coded number of edges to form a circle

    CEDGeoTubeParams returnParams;
    returnParams.Rmax = distance + 1.0*thickness;
    returnParams.Rmin = distance;
    //(see comment in line 151)
    returnParams.Rmax = returnParams.Rmax/cos(M_PI/NCircle);
    returnParams.Rmin = returnParams.Rmin/cos(M_PI/NCircle);
    //assume round edges
    returnParams.inner_symmetry = NCircle;
    returnParams.outer_symmetry = NCircle;
    returnParams.phi0 = 270. - 180./NCircle - (360./NCircle)*(int ((270.-180./NCircle)/(360./NCircle)));
    returnParams.delta_phi = 0.0;

    //cellSize1 is half the length along the z-axis (see e.g. Solenoid_o1_v01_gep.cpp line 58/59)
    returnParams.delta_z = cellSize1;
    //cellSize0 is defined to be the middle point of the geometry 
    returnParams.z0 = -cellSize1 + cellSize0;

    returnParams.isBarrel = true;

    return returnParams;
  }

//converts the parameters from a FixedPadSizeTPCData given by the appropriate drivers
//into those required by the CEDGeoTube
  CEDGeoTubeParams TPCParameterConversion(dd4hep::rec::FixedPadSizeTPCData *tpc){
    double zHalf = tpc->zHalf/dd4hep::mm;
    //these radii include the insensitive space!
    double rMin = tpc->rMin/dd4hep::mm;
    double rMax = tpc->rMax/dd4hep::mm;
  
    int NCircle = 36;   //hard coded number of edges to form a circle

    CEDGeoTubeParams returnParams;

    returnParams.Rmax = rMax;
    returnParams.Rmin = rMin;
    //(see comment in line 151)
    returnParams.Rmax = returnParams.Rmax/cos(M_PI/NCircle);
    returnParams.Rmin = returnParams.Rmin/cos(M_PI/NCircle);
    //assume round edges
    returnParams.inner_symmetry = NCircle;
    returnParams.outer_symmetry = NCircle;
    returnParams.phi0 = 270. - 180./NCircle - (360./NCircle)*(int ((270.-180./NCircle)/(360./NCircle)));
    returnParams.delta_phi = 0.0;

  
    returnParams.delta_z = zHalf;
    returnParams.z0 = -zHalf;

    returnParams.isBarrel = true;

    return returnParams;
  }

//converts the parameters from a ZPlanarData::LayerLayout layer given by the appropriate drivers
//into those required by the CEDGeoBox (for drawing of staves) or by CEDGeoTube (for approximation of the set of staves into tubes)
  LayerGeometry TrackerLayerParameterConversion(std::vector<dd4hep::rec::ZPlanarData::LayerLayout>::iterator thisLayer){
    int nLadders = thisLayer->ladderNumber;
    double phi0 = thisLayer->phi0*180/M_PI;
  
    double distance_sensitive = thisLayer->distanceSensitive/dd4hep::mm ;
    double thickness_sensitive = thisLayer->thicknessSensitive/dd4hep::mm ;
    double width_sensitive = thisLayer->widthSensitive/dd4hep::mm ;
    double offset_sensitive = thisLayer->offsetSensitive/dd4hep::mm ;
    double zHalf_sensitive = thisLayer->zHalfSensitive/dd4hep::mm ;  
  
    LayerGeometry Geometry;
    //two possibilites to draw: 
    //1) Draw layer consisting of staves
    double currentPhi; double radius;  double deltaPhi = 360./nLadders;
    for (int i=0; i<nLadders; i++){
      //placement shall begin along -z axis like for all other geometries
      currentPhi = phi0 + i*deltaPhi;
    
      //distance_sensitive is passed to be the minimal radius by the driver;
      //for drawing, it should be converted to the radius of the middle line
      radius = distance_sensitive+0.5* thickness_sensitive;

      CEDGeoBox stave;
      //place the center of the box at the appropriate coordinates
      //offset_sensitive is an additonal shift in placement orthogonal to the original 2D radius vector
      stave.center[0] = (radius*cos(currentPhi*M_PI/180) - offset_sensitive*sin(currentPhi*M_PI/180));
      stave.center[1] = (radius*sin(currentPhi*M_PI/180) + offset_sensitive*cos(currentPhi*M_PI/180));
      stave.center[2] = 0.0;  //placed z=0
      //dimensions are straight forward in an xyz coordinate system
      stave.sizes[0]  = thickness_sensitive;
      stave.sizes[1]  = width_sensitive;
      stave.sizes[2]  = zHalf_sensitive * 2;
      //the individual staves are finally rotated in the coordinate system to their appropriate position
      //herby, the rotation is performed around the z-axis
      stave.rotate[0] = 0.0;
      stave.rotate[1] = 0.0;
      stave.rotate[2] = currentPhi;

      Geometry.staves.push_back(stave);
    }


    //2) Summarize the set of staves into a geotube
    //analogous conversions as for the CalorimeterLayerParameterConversion (see above)
    Geometry.tube.Rmax = (distance_sensitive+thickness_sensitive)/cos(M_PI/nLadders); 
    Geometry.tube.Rmin = distance_sensitive/cos(M_PI/nLadders);
    Geometry.tube.inner_symmetry = nLadders;
    Geometry.tube.outer_symmetry = nLadders;
    Geometry.tube.phi0 = phi0 + 270. - 180./nLadders - (360./nLadders)*(int ((270.-180./nLadders)/(360./nLadders)));
    Geometry.tube.delta_phi = 0.0;
    Geometry.tube.delta_z = zHalf_sensitive;
    Geometry.tube.z0 = - zHalf_sensitive;

    //return both possible sets of parameters describing the geometry. The choice for either is implemented in the main draw routine
    return Geometry;
  }

//draws the given surfaces as a set of individual lines from indicated start- to the endpoint
  bool DrawSurfaces(const dd4hep::rec::SurfaceManager &surfMan, std::string detName, unsigned color, int layer){
    typedef dd4hep::rec::SurfaceMap SMap;
    const SMap* sMap = surfMan.map(detName);
    int lineCounter = 0;
    if(sMap) {
      for (SMap::const_iterator it = sMap->begin(); it != sMap->end(); ++it){
	dd4hep::rec::Surface* surf = dynamic_cast<dd4hep::rec::Surface*> (it->second);
	if (!surf) continue;
	if (!(surf->type().isVisible())) continue;
	const std::vector<std::pair<dd4hep::rec::Vector3D,dd4hep::rec::Vector3D> > lines = surf->getLines();
	if (lines.empty()){
	  debug( )<<" **** drawSurfaces(): empty lines vector for surface "<< *surf <<endmsg;
	  continue;
	}
	for(unsigned i = 0; i <lines.size(); i++){
	  unsigned default_width = 2;
	  unsigned default_type = layer;
	  ced_line(lines[i].first.x()/dd4hep::mm ,lines[i].first.y()/dd4hep::mm ,lines[i].first.z()/dd4hep::mm ,
		   lines[i].second.x()/dd4hep::mm ,lines[i].second.y()/dd4hep::mm ,lines[i].second.z()/dd4hep::mm ,
		   default_type,default_width, color);
	  ced_line_ID(lines[i].first.x()/dd4hep::mm ,lines[i].first.y()/dd4hep::mm ,lines[i].first.z()/dd4hep::mm ,
		      lines[i].second.x()/dd4hep::mm ,lines[i].second.y()/dd4hep::mm ,lines[i].second.z()/dd4hep::mm ,
		      default_type,default_width, color, layer);
	  lineCounter++; 
	}
      }
    } 
    if (lineCounter > 0)
      debug( )<<"Surfaces have been used."<<endmsg;
    return lineCounter > 0; //at least one line must have been drawn, otherwise the surface is considered to be undrawn
  }



  void getVisAttributes(dd4hep::DetElement det, unsigned &color, bool &visible) {
    dd4hep::VisAttr thisVisAttribute = det.volume().visAttributes();
    if (thisVisAttribute.isValid()){
      TColor* c = gROOT->GetColor( thisVisAttribute.color() );
      if ( c != nullptr ){
	//convert the given TColor into a hexadecimal whereby R_i, G_i, B_i integers
	//color = 0x00|R_1 R_2|G_1 G_2|B_1 B_2|
	color = ((int(255*c->GetRed())<<16) |  //multiply with 2^16 to get the last two bits
		 (int(255*c->GetGreen())<<8)|          //multiply with 2^8 to get the middle bits
		 (int(255*c->GetBlue())<<0));          //multiply with 2^0 to get the first two bits
	//the | operator is a bitwise addition of the number
      }
      else{
	color = ((int(255)<<16) | (int(255)<<8) | (int(255)<<0));
      }
    }
    else{
      debug()<<"color: pointer does not exist"<<endmsg;
      color  = 8947848; //== 0xff999999
      visible = true;   
    }
    //TODO: Hard coded to make every element visible for now
    visible = true;
  }


} /// end namespace
