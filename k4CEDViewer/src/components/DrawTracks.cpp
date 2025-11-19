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

#include "DD4hep/DetType.h"

#include <string>


typedef std::vector< edm4hep::TrackerHit > TrackerHitVec;

namespace k4ced{

  static std::vector<unsigned> colors = {
    0x510000,0x660202,0x720202,0x840202,0x960303,0xa80303,0xbc0303,
    0xce0404,0xe00404,0xef0404,0xff0000,0xfc0505,0xf91111,0xf92222,
    0xfc4949,0xf97777,0xf98e8e,0xf9aeae,0xf7b2b2,0xf7cdcd,
    
    0x512401,0x662d03,0x7a3602,0x934204,0xa54a04,0xb75103,0xc65803,
    0xd35e04,0xe56604,0xff6e00,0xf4710c,0xf4791a,0xf9842a,0xf98f3e,
    0xf99d57,0xf9a768,0xf9b47f,0xf9bf93,0xf9c9a4,0xf9d2b3,
    
    0x48014c,0x610266,0x700277,0x93039b,0xb103ba,0xc904d3,0xda04e5,
    0xe604f2,0xfa00ff,0xf00ffc,0xec1df7,0xef2ff9,0xeb42f4,0xec58f4,
    0xed6bf4,0xf486f9,0xf59af9,0xf8b5fc,0xf4c3f7,0xf8d6f9,
    
    0x2d0251,0x3d026d,0x4a0284,0x550399,0x5f03aa,0x6903bc,0x7102cc,
    0x8004e5,0x9800ff,0x8e0ef7,0x9922f9,0xa134f9,0xa845f9,0xb057f9,
    0xbc70f9,0xbf77f9,0xc98ef9,0xd3a4f9,0xddbbf9,0xecd6ff,
    
    0x00004f,0x020268,0x03037f,0x030399,0x0303b2,0x0404cc,0x0404e0,
    0x0404ef,0x0000ff,0x0c0cf9,0x1b1bf9,0x2a2af9,0x3939f9,0x4d4df9,
    0x6363f9,0x7272f9,0x8484f9,0x9898f9,0xb3b3f9,0xcacaf9,
    
    0x01434c,0x025c68,0x027382,0x04899b,0x0397aa,0x03abc1,0x04b9d1,
    0x04cbe5,0x00d8ff,0x0cdef9,0x1bddf7,0x2ae1f9,0x3ddff4,0x59e4f7,
    0x6be9f9,0x7cebf9,0x94f0fc,0xa3edf7,0xb6f2f9,0xc7f4f9,
    
    0x014443,0x01605f,0x027f7d,0x039996,0x03b2af,0x01c6c3,0x04ddda,
    0x04efeb,0x00ffff,0x09f9f5,0x0ef9f5,0x20f9f6,0x32fcf9,0x3ef9f6,
    0x52f9f7,0x6bf9f7,0x7ff9f7,0x95f9f8,0xb1f9f8,0xcaf9f9,
    
    0x016001,0x027002,0x027f02,0x029102,0x05aa05,0x05bf05,0x06d306,
    0x04e504,0x00ff00,0x09f909,0x18f918,0x2cf92c,0x43f943,0x52f952,
    0x63f963,0x77f977,0x8bf98b,0x9ff99f,0xb3f9b3,0xcaf9ca,
    
    0x344701,0x4b6603,0x608202,0x739b04,0x83b203,0x96cc04,0xa7e204,
    0xb1ef04,0xb6ff00,0xbaf713,0xbff725,0xc5f73b,0xcbf751,0xd3f968,
    0xd7f97a,0xd8f48b,0xe2f9a2,0xe1f4ad,0xe7f7bb,0xe9f4c8,
    
    0x565501,0x727002,0x898702,0xa5a303,0xb7b403,0xd1cd04,0xe2df04,
    0xefeb04,0xffff00,0xf9f509,0xf9f618,0xf9f62a,0xf7f43b,0xf9f64a,
    0xf9f759,0xf9f76b,0xf9f77c,0xf9f88b,0xf9f89f,0xfcfbbd
  };
  

  // helper function to get track states
  edm4hep::TrackState getTrackStateAt( const edm4hep::Track& trk , int location ){

    edm4hep::TrackState ts0 ;
    
    auto trkStates = trk.getTrackStates() ; 

    for( auto ts : trkStates){

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


//  Gaudi::Property<double> mcpECut{  this, "ECut"  , 0.01 , "energy cut for Tracks (in GeV) " };

  
  
  
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

      int color = colors[ i % colors.size() ] ;

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
