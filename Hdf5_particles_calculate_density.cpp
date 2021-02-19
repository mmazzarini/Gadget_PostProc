#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "library_MMazzarini_density.h"
#include "library_MMazzarini_structs.h"

/*
*
* This function calculates the center of mass one time, then repeats iteratively within a sphere of given radius that is shrinked
* every time
*
*/

void density_center_calc_MW(struct mass_center *dcms, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3,struct header_h5 *hp){

  float mySoft = 0.; /* Average softening */
  float myResol = 0.; /* Softening resolution */
  float myMass = 0.; /* Total mass */
  struct mass_center ocm; //old center of mass
  struct mass_center ncm; //new center of mass
  float maxRad = 0.;
  int flag = 0;	
  int nparts = 0;

  dcms_initialize(&ncm);  /* Initialize ncm to vector 0 */

  dcms_update(&ocm,&ncm);  /* Update ocm to ncm values */

  flag = 0;

  dcms_calculate_MW(&ncm,&ocm,pp1,pp2,pp3,hp,&maxRad,&nparts,flag);  /* Calculate center of mass and update to ncm */

  dcms_update(&ocm,&ncm); 
 
  /* Calculate average softening for resolution in density center calculation */	
 
  int i = 0;
  for (i = 0; i< hp->nbody1; i++){
      mySoft += (pp1->at(i).eps)*(pp1->at(i).mass);
      myMass += pp1->at(i).mass;
  }
  i = 0;
  for (i = 0; i< hp->nbody2; i++){
      mySoft += (pp2->at(i).eps)*(pp2->at(i).mass);
      myMass += pp2->at(i).mass;
  }
  i = 0;
  for (i = 0; i< hp->nbody3; i++){
      mySoft += (pp3->at(i).eps)*(pp3->at(i).mass);
      myMass += pp3->at(i).mass;
  }

  mySoft = mySoft/myMass;
  myResol = mySoft;
  std::cout << " " << std::endl;
  std::cout << "Working on MW with average softening resolution: " << myResol << std::endl;
  std::cout << " " << std::endl;
  
  int iter = 1;
  flag = 1;

  std::cout << "DEBUG ::: ocm is :  " << ocm.pos[0] << "\t" << ocm.pos[1] << "\t" << ocm.pos[2] << "\n" << ocm.vel[0] << "\t" << ocm.vel[1] << "\t"<< ocm.vel[2] << "\n" << std::endl;	

  while(( maxRad > myResol ) and (nparts >= 100)){

     std::cout << " Iteration number " << iter << std::endl;
     maxRad = 0.75*maxRad;
     dcms_initialize(&ncm);
 
     dcms_calculate_MW(&ncm,&ocm,pp1,pp2,pp3,hp,&maxRad,&nparts,flag);
     
     dcms_update(&ocm,&ncm);
        
     iter++;
     std::cout << " " << std::endl;
 
  } 
      
  std::cout << "End of iterations " << std::endl;
  std::cout << "Total number of iterations: " << iter << std::endl;
  std::cout << "Final particles within radius " << maxRad << " : " << nparts << std::endl; 

  dcms_determine(dcms,&ncm);  
 	
}

//The same but for satellite particles

void density_center_calc_Sat(std::vector<struct mass_center> *dcms_sat, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5,struct header_h5 *hp){

    int initPartDark = 0;
    int endPartDark = 0;
    int initPartStar = 0;
    int endPartStar = 0; 	    	  

    for(int ind_sat = 0;ind_sat<hp->sat_numb;ind_sat++){

  	if (ind_sat == 0){

		initPartDark = 0;
		endPartDark = hp->myNumbDark[0];
		initPartStar = 0;
		endPartStar = hp->myNumbStar[0];

	}else{
		
		initPartDark += hp->myNumbDark[ind_sat-1];
		endPartDark += hp->myNumbDark[ind_sat];
		initPartStar += hp->myNumbStar[ind_sat-1];
		endPartStar += hp->myNumbStar[ind_sat];	
 
	}

  	  std::cout << "Starting with dark particle " << initPartDark << "..." << std::endl;	
          std::cout << "Ending with dark particle " << endPartDark << "..." << std::endl;
  	  std::cout << "Starting with star particle " << initPartStar << "..." << std::endl;
          std::cout << "Ending with star particle " << endPartStar << "..." << std::endl;

	  float mySoft = 0.; /* Average softening */
    	  float myResol = 0.; /* Softening resolution */
    	  float myMass = 0.; /* Total mass */
    	  struct mass_center ocm; //old center of mass
    	  struct mass_center ncm; //new center of mass
    	  float maxRad = 0.;
    	  int flag = 0;	
    	  int nparts = 0;

	  dcms_initialize(&ncm);  /* Initialize ncm to vector 0 */

          dcms_update(&ocm,&ncm);  /* Update ocm to ncm values */
  
          flag = 0;

          dcms_calculate_Sat(&ncm,&ocm,pp4,pp5,hp,&maxRad,&nparts,flag,&ind_sat,&initPartDark,&endPartDark, &initPartStar,&endPartStar);  /* Calculate center of mass and update to ncm */

          dcms_update(&ocm,&ncm); 
 
          /* Calculate average softening for resolution in density center calculation */	

          int i;

          i = 0;
          for (i = initPartDark; i< endPartDark; i++){
          mySoft += (pp4->at(i).eps)*(pp4->at(i).mass);
          myMass += pp4->at(i).mass;
         }
 	 i = 0;

  	for (i = initPartStar; i< endPartStar; i++){
      		mySoft += (pp5->at(i).eps)*(pp5->at(i).mass);
      		myMass += pp5->at(i).mass;
  	}
  

  	mySoft = mySoft/myMass;
  	myResol = mySoft;
  	std::cout << " " << std::endl;
  	std::cout << std::endl << "Working on satellite " << hp->sat_numb+1 << " with average softening resolution: " << myResol << std::endl;
  	std::cout << " " << std::endl;
  
 	int iter = 1;
  	flag = 1;

  	std::cout << "DEBUG ::: ocm is :  " << ocm.pos[0] << "\t" << ocm.pos[1] << "\t" << ocm.pos[2] << "\n" << ocm.vel[0] << "\t" << ocm.vel[1] << "\t"<< ocm.vel[2] << "\n" << std::endl;	

  	while(( maxRad > myResol ) and (nparts >= 100)){

	std::cout << " Iteration number " << iter << std::endl;
     	maxRad = 0.75*maxRad;
     	dcms_initialize(&ncm);
 
     	dcms_calculate_Sat(&ncm,&ocm,pp4,pp5,hp,&maxRad,&nparts,flag,&ind_sat,&initPartDark,&endPartDark,&initPartStar,&endPartStar );
     
     	dcms_update(&ocm,&ncm);
        
     	iter++;
     	std::cout << " " << std::endl;
 
  	} 
      
  	std::cout << "End of iterations " << std::endl;
  	std::cout << "Total number of iterations: " << iter << std::endl;
  	std::cout << "Final particles within radius " << maxRad << " : " << nparts << std::endl; 

  	dcms_determine(&(dcms_sat->at(ind_sat)),&ncm);  
 
  } //End for loop
	
}

//For satellite Baryons --in the case of Superbox-original data one must know in advance, when launching the code, that its only ptype4 that gives the whole density center 

/*

void density_center_calc_SatelliteBaryon(struct mass_center *dcms_sat, std::vector<struct h5_particle> *pp4, struct header_h5 *hp){

  float mySoft = 0.; // Average softening 
  float myResol = 0.; // Softening resolution 
  float myMass = 0.; // Total mass
  struct mass_center ocm; //old center of mass
  struct mass_center ncm; //new center of mass
  float maxRad = 0.;
  int flag = 0;	
  int nparts = 0;

  dcms_initialize(&ncm);  // Initialize ncm to vector 0 

  dcms_update(&ocm,&ncm);  // Update ocm to ncm values 

  flag = 0;

  dcms_calculate_SatelliteBaryon(&ncm,&ocm,pp4,hp,&maxRad,&nparts,flag);  // Calculate center of mass and update to ncm 

  dcms_update(&ocm,&ncm); 
 
  // Calculate average softening for resolution in density center calculation 	

  i = 0;
  for (i = 0; i< hp->nbody4; i++){
      mySoft += (pp4->at(i).eps)*(pp4->at(i).mass);
      myMass += pp4->at(i).mass;
  }

  mySoft = mySoft/myMass;
  myResol = 2.*mySoft;
  std::cout << " " << std::endl;
  std::cout << "Working on satellite with average softening resolution: " << myResol << std::endl;
  std::cout << " " << std::endl;
  
  int iter = 1;
  flag = 1;

  std::cout << "DEBUG ::: ocm is :  " << ocm.pos[0] << "\t" << ocm.pos[1] << "\t" << ocm.pos[2] << "\n" << ocm.vel[0] << "\t" << ocm.vel[1] << "\t"<< ocm.vel[2] << "\n" << std::endl;	

  while(( maxRad > myResol ) and (nparts >= 100)){

     std::cout << " Iteration number " << iter << std::endl;
     maxRad = 0.75*maxRad;
     dcms_initialize(&ncm);
 
     dcms_calculate_SatelliteBaryon(&ncm,&ocm,pp4,hp,&maxRad,&nparts,flag);
     
     dcms_update(&ocm,&ncm);
        
     iter++;
     std::cout << " " << std::endl;
 
  } 
      
  std::cout << "End of iterations " << std::endl;
  std::cout << "Total number of iterations: " << iter << std::endl;
  std::cout << "Final particles within radius " << maxRad << " : " << nparts << std::endl; 

  dcms_determine(dcms_Baryon,&ncm);  
 	
}

*/

//For satellite DM --in the case of Superbox one must know in advance, when launching the code, that its only ptype4 that gives the whole bunch of particles

/*

void density_center_calc_SatelliteDark(struct mass_center *dcms_sat, std::vector<struct h5_particle> *pp5,struct header_h5 *hp){

  float mySoft = 0.; // Average softening 
  float myResol = 0.; // Softening resolution 
  float myMass = 0.; // Total mass 
  struct mass_center ocm; //old center of mass
  struct mass_center ncm; //new center of mass
  float maxRad = 0.;
  int flag = 0;	
  int nparts = 0;

  dcms_initialize(&ncm);  // Initialize ncm to vector 0 

  dcms_update(&ocm,&ncm);  // Update ocm to ncm values

  flag = 0;

  dcms_calculate_SatelliteDark(&ncm,&ocm,pp5,hp,&maxRad,&nparts,flag);  // Calculate center of mass and update to ncm

  dcms_update(&ocm,&ncm); 
 
  // Calculate average softening for resolution in density center calculation 	

  i = 0;
  for (i = 0; i< hp->nbody5; i++){
      mySoft += (pp5->at(i).eps)*(pp5->at(i).mass);
      myMass += pp5->at(i).mass;
  }

  mySoft = mySoft/myMass;
  myResol = 2.*mySoft;
  std::cout << " " << std::endl;
  std::cout << "Working on satellite with average softening resolution: " << myResol << std::endl;
  std::cout << " " << std::endl;
  
  int iter = 1;
  flag = 1;

  std::cout << "DEBUG ::: ocm is :  " << ocm.pos[0] << "\t" << ocm.pos[1] << "\t" << ocm.pos[2] << "\n" << ocm.vel[0] << "\t" << ocm.vel[1] << "\t"<< ocm.vel[2] << "\n" << std::endl;	

  while(( maxRad > myResol ) and (nparts >= 100)){

     std::cout << " Iteration number " << iter << std::endl;
     maxRad = 0.75*maxRad;
     dcms_initialize(&ncm);
 
     dcms_calculate_Dark(&ncm,&ocm,pp5,hp,&maxRad,&nparts,flag);
     
     dcms_update(&ocm,&ncm);
        
     iter++;
     std::cout << " " << std::endl;
 
  } 
      
  std::cout << "End of iterations " << std::endl;
  std::cout << "Total number of iterations: " << iter << std::endl;
  std::cout << "Final particles within radius " << maxRad << " : " << nparts << std::endl; 

  dcms_determine(dcms_Dark,&ncm);  
 	
}

*/




/*
*
* This function initializes the center of mass
*
*/


void dcms_intialize(struct mass_center *cm){

  int n;
  cm->mass   = 0.;
  for (n=0; n<3;n++){
  cm->pos[n] = 0.;
  cm->vel[n] = 0.;
  }
}


/*
*
* This function updates the old center of mass to the new one
*
*/


void dcms_update(struct mass_center *ocm, struct mass_center *ncm){

  int n;
  ocm->mass   = ncm->mass;
  for (n=0; n<3;n++){
  ocm->pos[n] = ncm->pos[n];
  ocm->vel[n] = ncm->vel[n];
  } 
}

/*
*
* This function is the density center core-calculation routine
*
*/

void dcms_calculate_MW(struct mass_center *ncm,struct mass_center *ocm, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3,struct header_h5 *hp,float *maxRad,int *nparts,const int flag){

    int i1 = 0;
    int i2 = 0;
    int i3 = 0;

    std::cout << "Debug max Radius: " << *maxRad << std::endl;	

    if (flag == 0)
    {
      for (i1 = 0; i1< hp->nbody1;  i1++){
          float myRad = std::sqrt( std::pow(pp1->at(i1).pos[0],2) + std::pow(pp1->at(i1).pos[1],2) + std::pow(pp1->at(i1).pos[2],2) );
          if( myRad > (*maxRad) ){ 
            *maxRad = myRad; 
          }		
          center_mass_increment_h5Particle(ncm,&(pp1->at(i1)));	  
          (*nparts)++;   
       } 

	std::cout << "Debug i1:  " <<  i1 << std::endl;

       for (i2 = 0; i2< hp->nbody2;  i2++){
          float myRad = std::sqrt( std::pow(pp2->at(i2).pos[0],2) + std::pow(pp2->at(i2).pos[1],2) + std::pow(pp2->at(i2).pos[2],2) );
          if( myRad > (*maxRad) ){ 
            *maxRad = myRad; 
          }			  
          center_mass_increment_h5Particle(ncm,&(pp2->at(i2)));	  
          (*nparts)++;     
       } 

	std::cout << "Debug i2:  " <<  i2 << std::endl;

       for (i3 = 0; i3< hp->nbody3;  i3++){
          float myRad = std::sqrt( std::pow(pp3->at(i3).pos[0],2) + std::pow(pp3->at(i3).pos[1],2) + std::pow(pp3->at(i3).pos[2],2) );
          if( myRad > (*maxRad) ){ 
            *maxRad = myRad; 
          }			  
          center_mass_increment_h5Particle(ncm,&(pp3->at(i3)));	  
          (*nparts)++;   
       } 

	std::cout << "Debug i3:  " <<  i3 << std::endl;

       int n;
       for(n = 0; n<3; n++){
          ncm->pos[n] = (ncm->pos[n])/(ncm->mass);
          ncm->vel[n] = (ncm->vel[n])/(ncm->mass);
       } 
       std::cout << "Number of particles within a sphere of radius " << *maxRad << " :   " << *nparts << std::endl;
   }

    else if (flag == 1){
       (*nparts) = 0;
       for (i1 = 0; i1< hp->nbody1;  i1++){
          float myRad = std::sqrt( std::pow(pp1->at(i1).pos[0]-ocm->pos[0],2) + std::pow(pp1->at(i1).pos[1]-ocm->pos[1],2) + std::pow(pp1->at(i1).pos[2]-ocm->pos[2],2) );
          if( myRad <= (*maxRad) ){ 
            center_mass_increment_h5Particle(ncm,&(pp1->at(i1)));	  
            (*nparts)++;   
          }
       }  

       for (i2 = 0; i2< hp->nbody2;  i2++){
          float myRad = std::sqrt( std::pow(pp2->at(i2).pos[0]-ocm->pos[0],2) + std::pow(pp2->at(i2).pos[1]-ocm->pos[1],2) + std::pow(pp2->at(i2).pos[2]-ocm->pos[2],2) );
          if( myRad <= (*maxRad) ){ 
            center_mass_increment_h5Particle(ncm,&(pp2->at(i2)));	  
           (*nparts)++;     
          }
       } 

       for (i3 = 0; i3< hp->nbody3;  i3++){
          float myRad = std::sqrt( std::pow(pp3->at(i3).pos[0]-ocm->pos[0],2) + std::pow(pp3->at(i3).pos[1]-ocm->pos[1],2) + std::pow(pp3->at(i3).pos[2]-ocm->pos[2],2) );
          if( myRad <= (*maxRad) ){ 
            center_mass_increment_h5Particle(ncm,&(pp3->at(i3)));	  
            (*nparts)++;   
          }
       } 

       int n;
       for(n = 0; n<3; n++){
          ncm->pos[n] = (ncm->pos[n])/(ncm->mass);
          ncm->vel[n] = (ncm->vel[n])/(ncm->mass);
       }   
     
       std::cout << "Number of particles within a sphere of radius " << *maxRad << " :   " << *nparts << std::endl;       

    }

}


//Same for satellite

void dcms_calculate_Sat(struct mass_center *ncm,struct mass_center *ocm, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp,float *maxRad,int *nparts,const int flag,int *ind_sat, int *initPartDark, int *finPartDark, int *initPartStar, int *finPartStar ){

    int i4 = 0;
    int i5 = 0;

    std::cout << "Debug max Radius: " << *maxRad << std::endl;	

    if (flag == 0)
    {
      for (i4 = *initPartDark; i4< *finPartDark; i4++){
          float myRad = std::sqrt( std::pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2) + std::pow(pp4->at(i4).pos[2],2) );
          if( myRad > (*maxRad) ){ 
            *maxRad = myRad; 
          }		
          center_mass_increment_h5Particle(ncm,&(pp4->at(i4)));	  
          (*nparts)++;   
       } 

	std::cout << "Debug i4:  " <<  i4 << std::endl;

      for (i5 = *initPartStar; i5< *finPartStar;  i5++){
          float myRad = std::sqrt( std::pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2) + std::pow(pp5->at(i5).pos[2],2) );
          if( myRad > (*maxRad) ){ 
            *maxRad = myRad; 
          }			  
          center_mass_increment_h5Particle(ncm,&(pp5->at(i5)));	  
          (*nparts)++;   
       } 

	

       //std::cout << "Debug i5:  " <<  i5 << std::endl;

       int n;
       for(n = 0; n<3; n++){
          ncm->pos[n] = (ncm->pos[n])/(ncm->mass);
          ncm->vel[n] = (ncm->vel[n])/(ncm->mass);
       } 
       std::cout << "Number of particles within a sphere of radius " << *maxRad << ":   " << *nparts << std::endl;
   }

    else if (flag == 1){
       (*nparts) = 0;
       for (i4 = *initPartDark; i4< *finPartDark; i4++){
          float myRad = std::sqrt( std::pow(pp4->at(i4).pos[0]-ocm->pos[0],2) + std::pow(pp4->at(i4).pos[1]-ocm->pos[1],2) + std::pow(pp4->at(i4).pos[2]-ocm->pos[2],2) );
          if( myRad <= (*maxRad) ){ 
            center_mass_increment_h5Particle(ncm,&(pp4->at(i4)));	  
            (*nparts)++;   
          }
       }  

       	
	
       for (i5 = *initPartStar; i5< *finPartStar;  i5++){
          float myRad = std::sqrt( std::pow(pp5->at(i5).pos[0]-ocm->pos[0],2) + std::pow(pp5->at(i5).pos[1]-ocm->pos[1],2) + std::pow(pp5->at(i5).pos[2]-ocm->pos[2],2) );
          if( myRad <= (*maxRad) ){ 
            center_mass_increment_h5Particle(ncm,&(pp5->at(i5)));	  
           (*nparts)++;     
          }
       } 

       	

       int n;
       for(n = 0; n<3; n++){
          ncm->pos[n] = (ncm->pos[n])/(ncm->mass);
          ncm->vel[n] = (ncm->vel[n])/(ncm->mass);
       }   
     
       std::cout << "Number of particles within a sphere of radius " << *maxRad << ":   " << *nparts << std::endl;       

    }

}


/*

//Same for baryons

void dcms_calculate_SatelliteBaryons(struct mass_center *ncm,struct mass_center *ocm, std::vector<struct h5_particle> *pp4, struct header_h5 *hp,float *maxRad,int *nparts,const int flag){

    int i4 = 0;

    std::cout << "Debug max Radius: " << *maxRad << std::endl;	

    if (flag == 0)
    {
      for (i4 = 0; i4< hp->nbody4;  i4++){
          float myRad = std::sqrt( std::pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2) + std::pow(pp4->at(i4).pos[2],2) );
          if( myRad > (*maxRad) ){ 
            *maxRad = myRad; 
          }		
          center_mass_increment_h5Particle(ncm,&(pp4->at(i4)));	  
          (*nparts)++;   
       } 

	std::cout << "Debug i4:  " <<  i4 << std::endl;

       int n;
       for(n = 0; n<3; n++){
          ncm->pos[n] = (ncm->pos[n])/(ncm->mass);
          ncm->vel[n] = (ncm->vel[n])/(ncm->mass);
       } 
       std::cout << "Number of particles within a sphere of radius " << *maxRad << ":   " << *nparts << std::endl;
   }

    else if (flag == 1){
       (*nparts) = 0;
       for (i4 = 0; i4< hp->nsph;  i4++){
          float myRad = std::sqrt( std::pow(pt4->at(i4).pos[0]-ocm->pos[0],2) + std::pow(pt4->at(i4).pos[1]-ocm->pos[1],2) + std::pow(pt4->at(i4).pos[2]-ocm->pos[2],2) );
          if( myRad <= (*maxRad) ){ 
            center_mass_increment_h5Particle(ncm,&(pt4->at(i4)));	  
            (*nparts)++;   
          }
       }  

       int n;
       for(n = 0; n<3; n++){
          ncm->pos[n] = (ncm->pos[n])/(ncm->mass);
          ncm->vel[n] = (ncm->vel[n])/(ncm->mass);
       }   
     
       std::cout << "Number of particles within a sphere of radius " << *maxRad << ":   " << *nparts << std::endl;       

    }

}


//Same for dark


void dcms_calculate_SatelliteDark(struct mass_center *ncm,struct mass_center *ocm, std::vector<struct h5_particle> *pp5, struct header_h5 *hp,float *maxRad,int *nparts,const int flag){

    int i5 = 0;

    std::cout << "Debug max Radius: " << *maxRad << std::endl;	

    if (flag == 0)
    {

       for (i5 = 0; i5< hp->nbody5;  i5++){
          float myRad = std::sqrt( std::pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2) + std::pow(pp5->at(i5).pos[2],2) );
          if( myRad > (*maxRad) ){ 
            *maxRad = myRad; 
          }			  
          center_mass_increment_h5Particle(ncm,&(pp5->at(i5)));	  
          (*nparts)++;   
       } 

	std::cout << "Debug i5:  " <<  i5 << std::endl;

       int n;
       for(n = 0; n<3; n++){
          ncm->pos[n] = (ncm->pos[n])/(ncm->mass);
          ncm->vel[n] = (ncm->vel[n])/(ncm->mass);
       } 
       std::cout << "Number of particles within a sphere of radius " << *maxRad << ":   " << *nparts << std::endl;
   }

    else if (flag == 1){
       (*nparts) = 0;

       for (i5 = 0; i5< hp->ndark;  i5++){
          float myRad = std::sqrt( std::pow(pt5->at(i5).pos[0]-ocm->pos[0],2) + std::pow(pt5->at(i5).pos[1]-ocm->pos[1],2) + std::pow(pt5->at(i5).pos[2]-ocm->pos[2],2) );
          if( myRad <= (*maxRad) ){ 
            center_mass_increment_h5Particle(ncm,&(pt5->at(i5)));	  
           (*nparts)++;     
          }
       } 

       int n;
       for(n = 0; n<3; n++){
          ncm->pos[n] = (ncm->pos[n])/(ncm->mass);
          ncm->vel[n] = (ncm->vel[n])/(ncm->mass);
       }   
     
       std::cout << "Number of particles within a sphere of radius " << *maxRad << ":   " << *nparts << std::endl;       

    }

}


*/


/*
*
* This function updates the density center to its final value
*
*/    

void dcms_determine(struct mass_center *dcms, struct mass_center *cm){

  int n;
  dcms->mass   = cm->mass;
  for (n=0; n<3;n++){
  dcms->pos[n] = cm->pos[n];
  dcms->vel[n] = cm->vel[n];
  }
 
}

/*
*
* This function writes the density center in a file with root names as the TIPSY file and with suffix "_Density_Center.ascii"
*
*/    

void density_center_write(std::string myTipsyFile,struct mass_center *dcms){

  std::string myDensityName;
  myDensityName = myTipsyFile+"_Density_Center.ascii";

  std::cout << "Writing density center on file " << myDensityName << std::endl;

  std::ofstream myDensityFile;
  myDensityFile.open(myDensityName.c_str());
  myDensityFile << std::scientific << std::showpos << dcms->pos[0] << "\t" << dcms->pos[1] << "\t" << dcms->pos[2] << "\n" << dcms->vel[0] << "\t" << dcms->vel[1] << "\t" << dcms->vel[2] << std::endl;
  myDensityFile.close();

  std::cout << "Density center written on file " << myDensityName << std::endl;

}

/*
*
* Functions to increment center_of_mass mass, position and velocity with corresponding data from h5 particles
*
*/


void center_mass_increment_h5Particle(struct mass_center *cm, struct h5_particle *pp){   

	    int n;
            cm->mass += pp->mass;
	    for (n=0;n<3;n++){
               cm->pos[n] += pp->pos[n]*pp->mass;
               cm->vel[n] += pp->vel[n]*pp->mass;
            }

}

void dcms_initialize(struct mass_center *cm){

	int n;
	cm->mass   = 0.;
	for (n = 0; n<3; n++){
	   cm->pos[n] = 0.;
	   cm->vel[n] = 0.;	
	}

}







