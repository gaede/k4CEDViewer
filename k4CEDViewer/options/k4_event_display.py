#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from Gaudi.Configuration import (DEBUG, INFO , WARNING)
from k4FWCore import ApplicationMgr
from k4FWCore.parseArgs import parser
from k4FWCore import IOSvc
#from Configurables import ExampleConsumer, DrawMCParticles
from Configurables import (GeoSvc,
                           DrawMCParticles,
                           EventDataSvc
                           )

algList = []
svcList = []


parser.add_argument(
    "--compactFile", help="Compact detector file to use", type=str, default=""
)

reco_args = parser.parse_known_args()[0]
compact_file = reco_args.compactFile   #// or get_compact_file_path(det_model)


evtsvc = EventDataSvc("EventDataSvc")
svcList.append(evtsvc)



iosvc = IOSvc("IOSvc")
iosvc.Input = "input_edm4hep.root"
#svcList.append(iosvc)


geoSvc = GeoSvc("GeoSvc")
geoSvc.detectors = [compact_file]
geoSvc.OutputLevel = INFO
geoSvc.EnableGeant4Geo = False
svcList.append(geoSvc)


algList.append(  DrawMCParticles("draw_mcps", colName = "MCParticles" , layer=0 , size=3 ) )


ApplicationMgr(TopAlg=algList,
               EvtSel="NONE",
               EvtMax=100,
               ExtSvc=svcList,
               OutputLevel=WARNING,
               )
