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


#include "edm4hep/ClusterCollection.h"
#include "k4FWCore/Consumer.h"
#include "k4GaudiCED.h"
#include "k4GaudiCEDUtils.h"
#include "k4CEDColors.h"
#include "ColorMap.h"

#include "DD4hep/DetType.h"

#include <string>


//typedef std::vector< edm4hep::Cluster >
typedef podio::RelationRange<edm4hep::Cluster> ClusterVec;
typedef podio::RelationRange<edm4hep::CalorimeterHit> CalorimeterHitVec;

using namespace k4ced ;


struct DrawClusters final : k4FWCore::Consumer<void(const edm4hep::ClusterCollection&)> {
  DrawClusters(const std::string& name, ISvcLocator* svcLoc)
    : Consumer(name, svcLoc,  KeyValue("colName", {"PandoraClusters"}) ) {

    k4GaudiCED::init(this) ;
  }

  
  Gaudi::Property<int> layer{ this, "layer" , 8 , "layer to draw Clusters " };
  Gaudi::Property<int> size{  this, "size"  , 2 , "size for drawning  Clusters " };
  Gaudi::Property<int> marker{  this, "marker"  , 0 , "marker for drawning  Clusters " };

  Gaudi::Property<unsigned> colorScheme{  this, "colorScheme"  , 10 ,
					  "Red:0,Orange:1,Plum:2,Violet:3,Blue:4,LightBlue:5,Aquamarine:6,Green:7,Olive:8,Yellow:9,Dark:10,Light:11,Classic:12" };


  Gaudi::Property<bool> drawEllipsoidForPFOClusters{ this, "DrawEllipsoidForPFOClusters" , false, 
						      "draw ellipsoids for clusters" };

  Gaudi::Property<bool> colorEnergy{ this, "ColorByEnergy" , false, 
				      "color recunstructed particle by energy" };

  Gaudi::Property<bool> colorEnergyAuto{ this, "ColorByEnergyAuto" , false, 
				      "Automatically adjust event by event the blue to min energy and red to max energy of event" };
  
  Gaudi::Property<double> colorEnergyMin{ this, "ColorByEnergyMin" , 0.,
					   "Minimal value for energy which will be represented as blue" };
  
  Gaudi::Property<double> colorEnergyMax{ this, "ColorByEnergyMax" , 10., 
					   "Maximal value for energy which will be represented as red" };
  
  Gaudi::Property<double> colorEnergySaturation{ this, "ColorByEnergySaturation" , 0.8, 
						  "Hue value that will be used to determine the pallete" };
  
  
  Gaudi::Property<double> colorEnergyValue{ this, "ColorByEnergyBrightnessn" , 0.8, 
					     "Brigtness value that will be used to determine the pallete" };

  
//===========================================================================================

  void operator()(const edm4hep::ClusterCollection& col) const override {



    k4ced::GlobalLog::instance().level()  = msgSvc()->outputLevel() ;
    k4ced::GlobalLog::instance().name()   = name() ;
    
    k4GaudiCED::newEvent(this) ;

    
    info()  <<  " +++++++  drawing Cluster collection with "
	    <<  col.size()  << " particles "
	    << " outputLevel = " <<   k4ced::GlobalLog::instance().level()
	    << endmsg ;

    unsigned myColID = PickingHandler::instance().colID() * k4ced::IDFactor ;

    printfun f =  PrintEDM4hep<edm4hep::ClusterCollection>( col )  ;
    PickingHandler::instance().registerFunctor( myColID/IDFactor , f ) ;

    dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

    Colors colors(colorScheme) ; 

    //------------------------
    
    //Determine the maximal and minimal cluster energy depositions in the event for color scaling (-->when drawing ellipsoids/cylinders).
    double eMin = 99999.; double eMax = 0;

    for( auto clu : col ){ 
      double e = clu.getEnergy();
      eMin = fmin(eMin, e);
      eMax = fmax(eMax, e);
    }

    unsigned ic = 0 ;
    for( auto cluster : col ){ 


      unsigned color = colors.current()[ ic++
					 % colors.size() ] ;


      double ene = cluster.getEnergy();
      
      if( colorEnergy ) {
	if( colorEnergyAuto ) {
	  color = ColorMap::NumberToTemperature(ene,eMin,eMax,colorEnergySaturation,colorEnergyValue);
	} else {
	  color = ColorMap::NumberToTemperature(ene,colorEnergyMin,colorEnergyMax,colorEnergySaturation,colorEnergyValue);
	}
      }
      CalorimeterHitVec hitvec = cluster.getHits();
      int nHits = (int)hitvec.size();

      for (int iHit = 0; iHit < nHits; ++iHit) {

	edm4hep::CalorimeterHit hit = hitvec[iHit];

	float x = hit.getPosition()[0];
	float y = hit.getPosition()[1];
	float z = hit.getPosition()[2];

	ced_hit_ID(x,y,z,marker, layer ,size,color,  myColID + cluster.id().index );
      }
    }

    //fg: this code below needs some work: the ellipsoids are flat and they should use the cluster
    //    shower parameters really, also this should probably be drawn in a different layer ....
    if( drawEllipsoidForPFOClusters ) {
      //refactored Cluster drawing as ellipsoids
      //by Thorben Quast, CERN Summer Student 2015

      for( auto cluster : col ){ 

	//Energy clusters are drawn as ellipsoids.
	//For each cluster, it's (energy weighted) central position, the deposited energy and the intrinsic direction in terms of sperical angles are given.
	//The minimal and maximal deposited energies among all clusters in the displayed event have been determined previously and will be needed for coloring.

	double cluster_center[] = {cluster.getPosition()[0], cluster.getPosition()[1], cluster.getPosition()[2]};
	double phi = cluster.getPhi();
	double theta = cluster.getITheta();
	//Use the direction of the cluster center w.r.t. the origin if no intrinsic angles are given.
	if( phi ==0. && theta==0.){
	  theta = atan( sqrt( cluster_center[0]*cluster_center[0] + cluster_center[1]*cluster_center[1] ) / cluster_center[2]  ) ;
	  phi = atan2( cluster_center[1] , cluster_center[0] ) ;
              }
              //Energy weighted moments of inertia are calculated. Ultimately, the eigenvalues of the 3x3 matrix will be a measure of the ellipsoids' extensions.
              CalorimeterHitVec hitvec = cluster.getHits();
              int nHits = (int)hitvec.size();
              double Etot = 0;
              double I[3][3]; for(int i=0; i<3; i++) for(int j=0; j<3; j++) I[i][j] = 0;
              //The angles theta and phi are used to transform the coordinates of each hit into a c.s. in which the x'-axis is parrallel to
              //the ellisoid's intrinsic direction. Note that this does not describe an unambigious system as any rotation along the x'-axis does not touch this constraint.
              //The transformation is achieved by a typical multiplication of rotation matrices: R(theta, phi) = R_y(Pi/2 - theta)*R_z(phi)
              double R[3][3];
              R[0][0] = cos(phi) * sin(theta); R[1][0] = -sin(phi); R[2][0] = -cos(phi)*cos(theta); R[0][1] = sin(phi)*sin(theta);
              R[1][1] = cos(phi); R[2][1] = -sin(phi)*cos(theta); R[0][2] = cos(theta); R[1][2] = 0; R[2][2] = sin(theta);
              double tot_x =0;
              double tot_y =0;
              double tot_z =0;
              for (int q = 0; q < nHits; q++){
		edm4hep::CalorimeterHit  hit = hitvec[q];
                float x = hit.getPosition()[0];
                float y = hit.getPosition()[1];
                float z = hit.getPosition()[2];
                float e = hit.getEnergy();
                //ced_hit_ID(x,y,z,marker, layer ,size,color,myColID + part.id().index);   //this line draws the indivdual hits within a cluster
                //translation and rotation of the coordinates
                x -= cluster_center[0];
                y -= cluster_center[1];
                z -= cluster_center[2];
                tot_x += x*e;
                tot_y += y*e;
                tot_z += z*e;

                double new_x = x * R[0][0] + y * R[0][1] + z * R[0][2];
                double new_y = x * R[1][0] + y * R[1][1] + z * R[1][2];
                double new_z = x * R[2][0] + y * R[2][1] + z * R[2][2];
                x = new_x; y = new_y; z = new_z;
                //calculate moments of inertia
                I[0][0] += x*x*e; I[1][1] += y*y*e; I[2][2] += z*z*e;
                I[0][1] = I[1][0] += x*y*e; I[0][2] = I[2][0] += x*z*e; I[1][2] = I[2][1] += y*z*e;
                Etot += e;
              }
              //The result of the rotation by the matrix R, the following coordinates correspond with each other:
              //  (component in coordinate system with x' || intrinsic direction)       (system in which ellipsoids are initially placed)
              //                          x'                                      <-->        z
              //                          y'                                      <-->        y
              //                          z'                                      <-->       -x
              //These assignments are corrected for by a modified rotation of the ellipsoid along its y-axis (see declaration of double rotate[])

              //I is not diagonal yet as only one axis was fixed when applying the rotation R.
              //The remaining lengths are determined by the solution of the 2x2 Eigenvalues (p-q formula).
              double lambda1 = 0.5*(I[2][2]+I[1][1]) + sqrt( pow(0.5*(I[2][2]+I[1][1]),2 ) + pow(I[2][1],2)-I[2][2]*I[1][1]);
              double lambda2 = 0.5*(I[2][2]+I[1][1]) - sqrt( pow(0.5*(I[2][2]+I[1][1]),2 ) + pow(I[2][1],2)-I[2][2]*I[1][1]);
              double sizes[3];
              sizes[0] = I[0][0]; sizes[1] = lambda1; sizes[2] = lambda2;
              //Remaining: (more or less) Arbitrary rescaling and transformation to a physical length (sqrt + energy division)
              for (int i=0; i<3; i++)  sizes[i] = sqrt(17.727)*sqrt(sizes[i])/Etot;

              double alpha = 0.5*asin(2*I[1][2]/fabs(lambda1-lambda2)) * (I[1][1]-I[2][2])/fabs(I[1][1]-I[2][2]);
              //We want to rotate the ellipsoid w.r.t. to the y-axis by 90deg - theta. Due to the rotation by R that maps x' <--> z, it is now
              //upside down such that an additional sign is needed.
              double rotate[] = {alpha, -(90-theta*180/M_PI), phi*180/M_PI};

              //The colors (blue and red) are set according the deposited energy in the cluster by comparison to other clusters in the event
              int ellipsoid_color = returnRGBClusterColor(cluster.getEnergy(), eMin, eMax, 256, 'a', 3);

              //Draw the ellipsoids, uncommenting the line with cylinders works as well.
              ced_ellipsoid_r(sizes, cluster_center, rotate, layer, ellipsoid_color);
              //ced_geocylinder_r(0.25*(sizes[0]+sizes[1]), sizes[2], cluster_center, rotate, 36, ellipsoid_color, layer);
      }
    }

    
    k4GaudiCED::draw(this, 1 );
  }


 int returnRGBClusterColor(float eneCluster, float cutoff_min, float cutoff_max, int color_steps, char scale, int colorMap) const {
    int color = 0x000000; //default colour: black
    int color_delta = 0; //colour step in the [0, color_steps] spectrum
    unsigned int rgb[] = {0, 0, 0}; //array of RGB to be returned as one 0x000000 HEX value

    /**
    * Check the input values for sanity */
    if (cutoff_min > cutoff_max) {
        std::cout << "Error in 'DSTViewer::returnRGBClusterColor': cutoff_min < cutoff_max" << std::endl;
    }
    if (eneCluster < 0.0) {
        std::cout << "Error in 'DSTViewer::returnRGBClusterColor': eneCluster is negative!" << std::endl;
    }
    if (cutoff_min < 0.0) {
        std::cout << "Error in 'DSTViewer::returnRGBClusterColor': eneCluster is negative!" << std::endl;
    }
    if (colorMap < 0 || colorMap > 6) {
        std::cout << "Error in 'DSTViewer::returnRGBClusterColor': wrong colorMap param!" << std::endl;
    }
    // Input values in log-scale
    float log_ene = std::log(eneCluster+1);
    float log_min = std::log(cutoff_min+1);
    float log_max = std::log(cutoff_max+1);
    float log_delta = log_max - log_min;
    float log_step = log_delta/(float)color_steps;

    switch(scale){
        case 'a': default: //log
            color_delta = (int) ((log_ene-log_min)/log_step); // which colour bin does the value go to? We have [0x00,0xFF] bins
            break;
        case 'b': //linear
            color_delta = (int)((eneCluster - cutoff_min)/(cutoff_max - cutoff_min)*color_steps);
            break;
    }


    if (color_delta >= color_steps){
        color_delta = color_steps;
    }
    if (color_delta < 0){
        color_delta = 0;
    }

    ColorMap::selectColorMap(colorMap)(rgb, color_delta, 0, color_steps);
    color = ColorMap::RGB2HEX(rgb[0],rgb[1],rgb[2]);

    return color;
} 



};

DECLARE_COMPONENT(DrawClusters)
