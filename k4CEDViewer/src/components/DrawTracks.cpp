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


#include "edm4hep/TrackCollection.h"
#include "k4FWCore/Consumer.h"
#include "k4GaudiCED.h"
#include "k4GaudiCEDUtils.h"
#include "k4CEDColors.h"

#include "DD4hep/DetType.h"

#include <string>


typedef std::vector< edm4hep::TrackerHit > TrackerHitVec;

namespace k4ced{


  // helper function to get track states
  edm4hep::TrackState getTrackStateAt( const edm4hep::Track& trk , int location ){

    edm4hep::TrackState ts0 ;

    for( auto ts : trk.getTrackStates() ){

      if( ts.location == location)
	return  ts ;
    }
    return ts0 ;
  } 

}

using namespace k4ced ;


struct DrawTracks final : k4FWCore::Consumer<void(const edm4hep::TrackCollection&)> {
  DrawTracks(const std::string& name, ISvcLocator* svcLoc)
    : Consumer(name, svcLoc,  KeyValue("colName", {"Track"}) ) {

    k4GaudiCED::init(this) ;
  }

  
  Gaudi::Property<int> layer{ this, "layer" , 5 , "layer to draw Tracks " };
  Gaudi::Property<int> size{  this, "size"  , 3 , "size for drawning  Tracks " };
  Gaudi::Property<int> marker{  this, "marker"  , 0 , "marker for drawning  Tracks " };
  Gaudi::Property<int> drawHelixForTracks{  this, "drawHelixForTracks"  , 1 , "draw a helix for Track objects: -1: none, 0: default, 1: atIP, 2: atFirstHit, 3: atLastHit, 4: atCalorimeter" };
  
  Gaudi::Property<float> lineThickness{  this, "lineThickness"  , 1. , "line thickness for drawning helix  " };
  Gaudi::Property<bool> useColorForHelixTracks{  this, "useColorForHelixTracks"  , 0 , "draw helices in the color of the track/PFO: 0: disabled (lightgrey), 1: enabled " };
  Gaudi::Property<unsigned> colorScheme{  this, "colorScheme"  , 0 ,
					  "Red:0,Orange:1,Plum:2,Violet:3,Blue:4,LightBlue:5,Aquamarine:6,Green:7,Olive:8,Yellow:9,Dark:10,Light:11,Classic:12" };
  
//===========================================================================================

  void operator()(const edm4hep::TrackCollection& col) const override {



    k4ced::GlobalLog::instance().level()  = msgSvc()->outputLevel() ;
    k4ced::GlobalLog::instance().name()   = name() ;
    
    k4GaudiCED::newEvent(this) ;


    info()  <<  " +++++++  drawing Track collection with " <<  col.size()  << " particles "
	    << " outputLevel = " <<   k4ced::GlobalLog::instance().level()
	    << endmsg ;

    unsigned myColID = PickingHandler::instance().colID() * k4ced::IDFactor ;

    printfun f =  PrintEDM4hep<edm4hep::TrackCollection>( col )  ;
    PickingHandler::instance().registerFunctor( myColID/IDFactor , f ) ;

    dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();


    //-----------------------

    // fixme ...
    double _helix_max_r = 2000. ;
    double _helix_max_z = 2300. ;


    Colors colors(colorScheme) ; 

    //------------------------
    
    for(unsigned i=0 ; i< col.size() ; ++i ){
      edm4hep::Track trk = col[ i ] ;

      // -- collect hits from all track segments
      TrackerHitVec tHV ;
      
      debug() << " -- track has "<< trk.getTracks().size() << " subtracks - will use these for displaying hits "
	      << endmsg ;

      std::copy( trk.getTrackerHits().begin() , trk.getTrackerHits().end() , std::back_inserter(  tHV ) ) ;

      for( unsigned j=0 ,N = trk.getTracks().size() ; j<N ; ++j ){
	auto t = trk.getTracks()[j] ;

	debug() << " -- track j= "<< j << " has " <<  t.getTrackerHits().size() << " hits  "
		<< endmsg ;
	
	std::copy( t.getTrackerHits().begin() , t.getTrackerHits().end() , std::back_inserter(  tHV ) ) ;

      }
      const TrackerHitVec& hits = tHV ;

      int color = colors.current()[ i % colors.size() ] ;

      k4GaudiCED::add_layer_description( inputLocations(0)[0] , layer);

      for( TrackerHitVec::const_iterator it = hits.begin();  it != hits.end() ; it++ ) {

	ced_hit_ID( (*it).getPosition()[0],
		    (*it).getPosition()[1],
		    (*it).getPosition()[2],
		    marker, layer , size , color , myColID + trk.id().index ) ;

      } 
      edm4hep::TrackState ts ;

      switch( drawHelixForTracks ){
	// Case 0 is the default value which takes the whatever first track state
	// (e.g. InteractionPoint for ILD, or the simple track for test beam data)
      case 0:
	if (trk.getTrackStates().size() > 0) {
	  ts = trk.getTrackStates().at(0) ;
	}
	break ;

      case 1:
	ts = getTrackStateAt( trk , edm4hep::TrackState::AtIP          ) ;
	break ;

      case 2:
	ts = getTrackStateAt( trk, edm4hep::TrackState::AtFirstHit    ) ;
	break ;

      case 3: // AtLastHit

	ts = (  trk.getTracks().empty() ?  getTrackStateAt( trk, edm4hep::TrackState::AtLastHit )
		:       getTrackStateAt( trk.getTracks().back(), edm4hep::TrackState::AtLastHit ) )   ;
	break ;

      case 4:  // AtCalo

	ts = (  trk.getTracks().empty() ?  getTrackStateAt( trk, edm4hep::TrackState::AtCalorimeter )
		:       getTrackStateAt( trk.getTracks().back(), edm4hep::TrackState::AtCalorimeter ) )   ;


	ced_hit_ID( ts.referencePoint[0],
		    ts.referencePoint[1],
		    ts.referencePoint[2],
		    1 , layer , size*20 , color , myColID + trk.id().index ) ;

	break ;
      }


      std::vector<double> bFieldVector(3) ;
      theDetector.field().combinedMagnetic(dd4hep::Position(0,0,0), &

					   bFieldVector[0] ) ;
      double bField = bFieldVector[2] / dd4hep::tesla;

      double pt;
      if (bField != 0.0 && std::abs(ts.omega) > 0.00001 )
	pt = bField * 3e-4 / std::abs( ts.omega ) ;
      else
	pt = 1.e10;
      double charge = ( ts.omega > 0. ?  1. : -1. ) ;
      double px = pt * std::cos(  ts.phi ) ;
      double py = pt * std::sin(  ts.phi ) ;
      double pz = pt * ts.tanLambda ;
      
      // start point for drawing ( PCA to reference point )
      double xs = ts.referencePoint[0] -  ts.D0 * sin( ts.phi ) ;
      double ys = ts.referencePoint[1] +  ts.D0 * cos( ts.phi ) ;
      double zs = ts.referencePoint[2] +  ts.Z0 ;
      
      
      if( drawHelixForTracks >= 0 && pt > 0.01 ) {
	
	int helixColor = ( useColorForHelixTracks ? color : 0xdddddd ) ;
	
	k4GaudiCED::drawHelix( bField , charge, xs, ys, zs ,
			       px, py, pz, layer ,  lineThickness , helixColor  ,
			       0.0, _helix_max_r,
			       _helix_max_z, myColID + trk.id().index ) ;
	
      }
    }
    
    k4GaudiCED::draw(this, 1 );
  }




};

DECLARE_COMPONENT(DrawTracks)
