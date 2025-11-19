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

struct DrawDetector final : k4FWCore::Consumer<void(const edm4hep::MCParticleCollection&)> {
  DrawDetector(const std::string& name, ISvcLocator* svcLoc)
    : Consumer(name, svcLoc,  KeyValue("colName", {"MCParticles"}) ) {

    k4GaudiCED::init(this) ;
  }

  Gaudi::Property<bool> drawSurfaces{  this, "drawSurfaces"  , false , " draw detector surfaces if available" };
  Gaudi::Property<std::vector<std::string>> drawDetailed{  this, "drawDetailed"  , {} , " draw these detectors with more details (e.g. staves" };

  
//===========================================================================================

  void operator()(const edm4hep::MCParticleCollection& ) const override {



    k4ced::GlobalLog::instance().level()  = msgSvc()->outputLevel() ;
    k4ced::GlobalLog::instance().name()   = name() ;
    
    k4GaudiCED::newEvent(this) ;


    info()  <<  " +++++++  drawing the detector  "
	    << " outputLevel = " <<   k4ced::GlobalLog::instance().level()
	    << endmsg ;


    dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();   

    //------------------------

    k4GaudiCED::drawDD4hepDetector(theDetector, drawSurfaces, drawDetailed ) ;

    //------------------------

    
    k4GaudiCED::draw(this, 1 );
  }




};

DECLARE_COMPONENT(DrawDetector)
