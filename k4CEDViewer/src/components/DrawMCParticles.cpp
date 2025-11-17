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

#define MCPARTICLE_LAYER 7

#include "edm4hep/MCParticleCollection.h"
#include "k4FWCore/Consumer.h"
#include "k4GaudiCED.h"
#include "k4GaudiCEDUtils.h"

using namespace k4ced ;

#include <string>

struct DrawMCParticles final : k4FWCore::Consumer<void(const edm4hep::MCParticleCollection&)> {
  DrawMCParticles(const std::string& name, ISvcLocator* svcLoc)
    : Consumer(name, svcLoc,  KeyValue("colName", {"MCParticle"}) ) {

    k4GaudiCED::init(this) ;
  }
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
  std::vector<DrawParameters> _drawParameters ;
  
  Gaudi::Property<int> m_layer{this, "layer", 0 , "layer to draw MCParticles "};

  Gaudi::Property<std::string> m_ecalBarrelName{this, "EcalBarrelName", "EcalBarrel" , "name of ecal barrel detector "};
  Gaudi::Property<std::string> m_ecalEndcapName{this, "EcalEndcapName", "EcalEndcap" , "name of ecal endcap detector "};
  Gaudi::Property<std::string> m_hcalBarrelName{this, "HcalBarrelName", "HcalBarrel" , "name of hcal barrel detector "};
  Gaudi::Property<std::string> m_hcalEndcapName{this, "HcalEndcapName", "HcalEndcap" , "name of hcal endcap detector "};


  Gaudi::Property<std::string> m_colName{this, "colName", "MyMCParticles" , "name of the MCPArticle collection "};
  
  
  void operator()(const edm4hep::MCParticleCollection& col) const override {


//    info() <<  " +++++++  drawing MCParticle collection : " <<  col << endmsg;
    std::cout  <<  " +++++++  drawing MCParticle collection with " <<  col.size()  << " particles "<< std::endl;

    k4GaudiCED::newEvent(this) ;

    dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
    
    // define some consts that should be made parameters ...
    bool usingParticleGun = false; 
    int layer = m_layer ;
    bool drawMCParticlesCreatedInSimulation = true ;
    double mcpECut = 0.01 ;
    double _helix_max_r = 2000. ;
    double _helix_max_z = 2500. ;

    int marker = 3 ;
    int size = 3 ;
    float scaleLineThickness =1. ;
      
    for(unsigned i=0; i< col.size() ; i++){

      auto mcp = col[i]  ;
      
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
	

//	layer = ( layer > -1 ? layer : MCPARTICLE_LAYER ) ;
//        this->drawParameters[np].Layer = layer ;
	
        k4GaudiCED::add_layer_description( m_colName.value() , layer );

        double px = mcp.getMomentum()[0];
        double py = mcp.getMomentum()[1];
        double pz = mcp.getMomentum()[2];
        double pt = sqrt(pow(px,2)+pow(py,2));
        double p = sqrt(pow(pt,2)+pow(pz,2));
        double x = mcp.getVertex()[0] ;
        double y = mcp.getVertex()[1] ;
        double z = mcp.getVertex()[2] ;

        if( std::fabs( charge ) > 0.0001  ) {

	  std::vector<double> bFieldVector = { 0., 0., 3.5 } ;   // FIXME: get from dd4hep
	  // theDetector.field().combinedMagnetic(Position(0,0,0), bFieldVector) ;

	  
	  double bField = bFieldVector[2] ;    /// FIXME:  / dd4hep::tesla;

	  debug() << "  drawing MCParticle helix for p_t "
		    << sqrt(px*px+py*py)
		    << endmsg ;
            const int ml = marker | ( layer << CED_LAYER_SHIFT ) ;
            //maximal extension of all charged tracks
            double _hmr, _hmz;
            switch(std::abs(mcp.getPDG() ) ){
                case 13:
                    _hmr = getCalorimeterParameters(theDetector, m_hcalBarrelName).r_inner + getCalorimeterParameters(theDetector, m_hcalBarrelName).delta_r;
                    _hmz = getCalorimeterParameters(theDetector, m_hcalEndcapName).z_0 + getCalorimeterParameters(theDetector, m_hcalEndcapName).delta_z;
                    break;
                default:
                    _hmr = _helix_max_r;
                    _hmz = _helix_max_z;
            }
            k4GaudiCED::drawHelix( bField , charge, x, y, z,
				   px, py, pz, ml , size*scaleLineThickness , 0x7af774  ,
				   0.0,  _hmr ,
				   _hmz, mcp.id().index ) ;

        } else { // neutral
            int color  ;
            double length, yokeR, yokeZ;
            switch(  std::abs(mcp.getPDG() )  ){
                //refactored length calculation (T. Quast 7 Aug 15)
                case 22:
                    color = 0xf9f920;          // photon
                    length = calculateTrackLength(m_ecalBarrelName, m_ecalEndcapName, theDetector, x, y, z, px, py, pz);
                    break ;
                case 12:  case 14: case 16: // neutrino
                    color =  0xdddddd  ;
                    yokeR = getYokeExtent(theDetector)[0];
                    yokeZ = getYokeExtent(theDetector)[1];
                    length = (fabs(pt/pz) > yokeR/yokeZ) ?
                            yokeR * sqrt(1. + pow(pz/pt,2)):
                            yokeZ * sqrt(1. + pow(pt/pz,2));
                    break ;
                default:
                    color = 0xb900de  ;        // neutral hadron
                    length = calculateTrackLength(m_hcalBarrelName, m_hcalEndcapName, theDetector, x, y, z, px, py, pz);
            }
            //tracks with vertex outside the according calorimeter are not drawn, length is passed as 0
            ced_line_ID( x , y , z ,
                        x + length*px/p ,  y + length*py/p ,  z + length*pz/p ,
			 layer  , size, color, mcp.id().index );
        }
    }

    
    k4GaudiCED::draw(this, 1 );
  }




};

DECLARE_COMPONENT(DrawMCParticles)
