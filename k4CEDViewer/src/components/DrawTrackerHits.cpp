/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/MCParticle.h"
#include "k4FWCore/Consumer.h"
#include "k4GaudiCED.h"
#include "k4GaudiCEDUtils.h"
#include "k4CEDColors.h"

#include "DD4hep/DetType.h"

#include <string>



using namespace k4ced ;


struct DrawTrackerHits final : k4FWCore::Consumer<void(const std::vector<const edm4hep::TrackerHit3DCollection*>&,
						       const std::vector<const edm4hep::TrackerHitPlaneCollection*>&)> {
  DrawTrackerHits(const std::string& name, ISvcLocator* svcLoc)
    : Consumer(name, svcLoc,  { KeyValues("colNamesTH3D", {"TPCTrackerHits"} ), KeyValues("colNamesTHPlane", {"VertexBarrelTrackerHits","SETTrackerHits"} ) }) {
    
    k4GaudiCED::init(this) ;
  }
  
  
  Gaudi::Property<int> layer{ this, "layer" , 11 , "layer to draw TrackerHits " };
  Gaudi::Property<int> size{  this, "size"  , 2 , "size for drawning  TrackerHits " };
  Gaudi::Property<int> marker{  this, "marker"  , 0 , "marker for drawning  TrackerHits " };
  Gaudi::Property<int> color{  this, "color"  , 0xee0000 , "color for drawning  TrackerHits (default: 0xee0000)" };

  
//===========================================================================================

  void operator()(const std::vector<const edm4hep::TrackerHit3DCollection*>& cols3D, const std::vector<const edm4hep::TrackerHitPlaneCollection*>& colsPlane) const override {
    


    k4ced::GlobalLog::instance().level()  = msgSvc()->outputLevel() ;
    k4ced::GlobalLog::instance().name()   = name() ;
    
    k4GaudiCED::newEvent(this) ;

    std::stringstream sstr ;
    info()  <<  " +++++++  drawing TrackerHits from collections: \n" ;
    for(unsigned i=0 ; i< inputLocations(0).size() ; ++i){

      sstr <<  inputLocations(0)[i] << ", " ;
      info() << "    " << inputLocations(0)[i] << "\n" ; 
    }
    for(unsigned i=0 ; i< inputLocations(1).size() ; ++i){

      sstr <<  inputLocations(1)[i] << ", " ;
      info() << "    " << inputLocations(1)[i] << "\n" ; 
    }

    info() << endmsg ;

    
    //-----------------------

    k4GaudiCED::add_layer_description( sstr.str(), layer);

    for( const auto* col : cols3D ) {

//      debug() << "  will draw these hits : " << *col << endmsg ;

      unsigned myColID = PickingHandler::instance().colID() * k4ced::IDFactor ;
      printfun f =  PrintEDM4hep<edm4hep::TrackerHit3DCollection>( *col )  ;
      PickingHandler::instance().registerFunctor( myColID/IDFactor , f ) ;
      
      for( int i=0, n=col->size(); i<n ; i++ ){
	
	auto h = col->at(i) ;

	int id =  myColID + h.id().index ;
	  
	  ced_hit_ID( h.getPosition()[0],
		      h.getPosition()[1],
		      h.getPosition()[2],
		      marker,layer, size , color, id ) ;
      }
    }

    for( const auto* col : colsPlane ) {

//      debug() << "  will draw these hits : " << *col << endmsg ;

      unsigned myColID = PickingHandler::instance().colID() * k4ced::IDFactor ;
      printfun f =  PrintEDM4hep<edm4hep::TrackerHitPlaneCollection>( *col )  ;
      PickingHandler::instance().registerFunctor( myColID/IDFactor , f ) ;
      
      for( int i=0, n=col->size(); i<n ; i++ ){
	
	auto h = col->at(i) ;
	
	int id =  myColID + h.id().index ;
	
	ced_hit_ID( h.getPosition()[0],
		    h.getPosition()[1],
		    h.getPosition()[2],
		    marker,layer, size , color, id ) ;
      }
    }
    
    k4GaudiCED::draw(this, 1 );
  }
      
};

DECLARE_COMPONENT(DrawTrackerHits)
