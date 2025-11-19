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


#include "edm4hep/MCParticleCollection.h"
#include "k4FWCore/Consumer.h"
#include "k4GaudiCED.h"
#include "k4GaudiCEDUtils.h"


#include "DD4hep/DetType.h"

using namespace k4ced ;

#include <string>

struct DrawMCParticles final : k4FWCore::Consumer<void(const edm4hep::MCParticleCollection&)> {
  DrawMCParticles(const std::string& name, ISvcLocator* svcLoc)
    : Consumer(name, svcLoc,  KeyValue("colName", {"MCParticle"}) ) {

    k4GaudiCED::init(this) ;
  }

  std::vector<DrawParameters> _drawParameters ;
  
  Gaudi::Property<int> layer{ this, "layer" , 0 , "layer to draw MCParticles " };
  Gaudi::Property<int> size{  this, "size"  , 3 , "size for drawning  MCParticles " };
  Gaudi::Property<double> mcpECut{  this, "ECut"  , 0.01 , "energy cut for MCParticles (in GeV) " };
  Gaudi::Property<bool> usingParticleGun{  this, "isParticleGun"  , false , " should be true of file contains events produced with a gun" };
  Gaudi::Property<bool> drawMCParticlesCreatedInSimulation{  this, "drawMCParticlesCreatedInSimulation"  , false , "draw particles created in simulation ? " };

  
  Gaudi::Property<std::string> m_ecalBarrelName{this, "EcalBarrelName", "EcalBarrel" , "name of ecal barrel detector "};
  Gaudi::Property<std::string> m_ecalEndcapName{this, "EcalEndcapName", "EcalEndcap" , "name of ecal endcap detector "};
  Gaudi::Property<std::string> m_hcalBarrelName{this, "HcalBarrelName", "HcalBarrel" , "name of hcal barrel detector "};
  Gaudi::Property<std::string> m_hcalEndcapName{this, "HcalEndcapName", "HcalEndcap" , "name of hcal endcap detector "};


  
//===========================================================================================

  void operator()(const edm4hep::MCParticleCollection& col) const override {



    k4ced::GlobalLog::instance().level()  = msgSvc()->outputLevel() ;
    k4ced::GlobalLog::instance().name()   = name() ;
    
    k4GaudiCED::newEvent(this) ;


    info()  <<  " +++++++  drawing MCParticle collection with " <<  col.size()  << " particles "
	    << " outputLevel = " <<   k4ced::GlobalLog::instance().level()
	    << endmsg ;

    unsigned myColID = PickingHandler::instance().colID() * k4ced::IDFactor ;

    printfun f =  PrintEDM4hep<edm4hep::MCParticleCollection>( col )  ;
    PickingHandler::instance().registerFunctor( myColID/IDFactor , f ) ;

    dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

    

    bool _surfaces = false ;

    //------------------------

    k4GaudiCED::drawDD4hepDetector(theDetector, _surfaces, {} ) ;

    //------------------------
    // get some calorimeter parameters for drawing the lines and helices

    // could directly get the calorimeter in question - yet this would need some error handling ....
    // dd4hep::DetElement ecalBarrel =  theDetector.detectors( (dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL |  dd4hep::DetType::ELECTROMAGNETIC ), dd4hep::DetType::AUXILIARY )[0] ;
    // dd4hep::DetElement ecalEndcap =  theDetector.detectors( (dd4hep::DetType::CALORIMETER | dd4hep::DetType::ENDCAP |  dd4hep::DetType::ELECTROMAGNETIC ), dd4hep::DetType::AUXILIARY )[0] ;
    // dd4hep::DetElement hcalBarrel =  theDetector.detectors( (dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL |  dd4hep::DetType::HADRONIC ), dd4hep::DetType::AUXILIARY )[0] ;
    // dd4hep::DetElement hcalEndcap =  theDetector.detectors( (dd4hep::DetType::CALORIMETER | dd4hep::DetType::ENDCAP |  dd4hep::DetType::HADRONIC ), dd4hep::DetType::AUXILIARY )[0] ;

						    
     
    CalorimeterDrawParams ecalBarrelParams =  getCalorimeterParameters(theDetector, m_ecalBarrelName) ;
    CalorimeterDrawParams ecalEndcapParams =  getCalorimeterParameters(theDetector, m_ecalEndcapName) ;
    CalorimeterDrawParams hcalBarrelParams =  getCalorimeterParameters(theDetector, m_hcalBarrelName) ;
    CalorimeterDrawParams hcalEndcapParams =  getCalorimeterParameters(theDetector, m_hcalEndcapName) ;
     
    

    //-----------------------
    for( edm4hep::MCParticle mcp : col ){

      float charge = mcp.getCharge ();

      if(mcp.getGeneratorStatus() != 1 // do not draw unstable particles

	 && usingParticleGun == false  // unless we are in particleGun mode

	 && !(drawMCParticlesCreatedInSimulation && mcp.isCreatedInSimulation()) // or we want to draw particles created in the simulation
        ) {
	continue;
      }

      if ( mcp.getEnergy() < mcpECut )
	continue ;         
      
      debug()  << "  drawing MCParticle pdg "
	       << mcp.getPDG()
	       << " genstat: " << mcp.getGeneratorStatus()
	       << endmsg ;
      
      
      k4GaudiCED::add_layer_description( inputLocations(0)[0] , layer );
      
      double px = mcp.getMomentum()[0];
      double py = mcp.getMomentum()[1];
      double pz = mcp.getMomentum()[2];
      double pt = sqrt(pow(px,2)+pow(py,2));
      double p = sqrt(pow(pt,2)+pow(pz,2));
      double x = mcp.getVertex()[0] ;
      double y = mcp.getVertex()[1] ;
      double z = mcp.getVertex()[2] ;
      
      if( std::fabs( charge ) > 0.0001  ) {
	  
	std::vector<double> bFieldVector(3) ;
	theDetector.field().combinedMagnetic(dd4hep::Position(0,0,0), &bFieldVector[0]) ;
	double bField = bFieldVector[2]  / dd4hep::tesla;
	  
	debug() << "  drawing MCParticle helix for p_t "
		<< sqrt(px*px+py*py)
		<< endmsg ;
	  
	  
	//maximal extension of all charged tracks - use ecal for all but muons
	  
	double _hmr = ecalBarrelParams.r_inner + ecalBarrelParams.delta_r ; 
	double _hmz = ecalEndcapParams.z_0 + ecalEndcapParams.delta_z ;
	  
	if( std::abs(mcp.getPDG()) == 13 ){
	    
	  _hmr = hcalBarrelParams.r_inner + hcalBarrelParams.delta_r ;
	  _hmz = hcalEndcapParams.z_0 + hcalEndcapParams.delta_z;
	    
	}
	  
	k4GaudiCED::drawHelix( bField , charge, x, y, z,
			       px, py, pz, layer , size , 0x7af774  ,
			       0.0,  _hmr ,
			       _hmz,  myColID + mcp.id().index   );
	  
      } else { // neutral
	int color  ;
	double length, yokeR, yokeZ;
	  
	auto ext =  getYokeExtent(theDetector) ;
	  
	switch(  std::abs(mcp.getPDG() )  ){
	  //refactored length calculation (T. Quast 7 Aug 15)
	case 22:
	  color = 0xf9f920;          // photon
	  length = calculateTrackLength( ecalBarrelParams, ecalEndcapParams, x, y, z, px, py, pz);
	  break ;
	case 12:  case 14: case 16: // neutrino
	  color =  0xdddddd  ;
	  yokeR = ext[0];
	  yokeZ = ext[1];
	  length = (fabs(pt/pz) > yokeR/yokeZ) ?
	    yokeR * sqrt(1. + pow(pz/pt,2)):
	    yokeZ * sqrt(1. + pow(pt/pz,2));
	  break ;
	default:
	  color = 0xb900de  ;        // neutral hadron
	  length = calculateTrackLength( hcalBarrelParams, hcalEndcapParams, x, y, z, px, py, pz);
	}
	//tracks with vertex outside the according calorimeter are not drawn, length is passed as 0
	ced_line_ID( x , y , z ,
		     x + length*px/p ,  y + length*py/p ,  z + length*pz/p ,
		     layer  , size, color,   myColID + mcp.id().index  );
      }
    }

    
    k4GaudiCED::draw(this, 1 );
  }




};

DECLARE_COMPONENT(DrawMCParticles)
