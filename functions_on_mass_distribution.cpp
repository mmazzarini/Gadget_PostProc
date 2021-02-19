#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "library_MMazzarini_structs.h"
#include "library_MMazzarini_mass_distribution.h"
#include "library_MMazzarini_utilities.h"
//#include "externalModule_spline.h"

/* Function to calculate the fraction of bound satellite mass...*/

void bound_mass_calculate(std::vector<struct mass_center> *dcms_Sat, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp){

	std::cout << "I am doing this" << std::endl;
	std::cout << "\n\n\n\nNbodySat is " << std::endl;
	

	for (int i=0;i<hp->sat_numb; i++){
	
		std::cout << "Dark " << hp->myNumbDark[i] << " " << "Star " << hp->myNumbStar[i] << std::endl;	
	 	hp->nbodySat[i] = int(hp->myNumbDark[i]) + int(hp->myNumbStar[i]);
			std::cout << "Total " << hp->nbodySat[i] << std::endl;
	}
	
	std::cout << "\n\n\n\n\n" << std::endl;

	int debug_tot_part = 0;

	for (int i=0;i<hp->sat_numb; i++){

		debug_tot_part += hp->nbodySat[i];

	}

	std::cout << "tot particles from array " << debug_tot_part << std::endl; 	
	std::cout << "while my array of part4 has size " << pp4->size() << std::endl; 
	std::cout << "And my array of part5 has size " << pp5->size() << std::endl;
	for(int i=0; i<hp->sat_numb; i++){

		std::cout << "The index i is " << i << std::endl;

		float totMass = 0.; //Total sum
		float totMass_p4 = 0.; // only dark or p4 particles
		float totMass_p5 = 0.; //only baryons or p5 particles
		//float old_fracBound = 0.; //This is to update fraction of bound particles and check convergence
		float fracBound = 0.;
		float fracBound_p4 = 0.;
		float fracBound_p5 = 0.;
		int i4;
		int i5;
		//float err = 1000.;
		//int iter = 1;
		int initPP4, initPP5;
		int endPP4,endPP5;


		std::cout << "Qui1" << std::endl;
		if (i == 0){

			initPP4 = 0;
			endPP4 = hp->myNumbDark[0];
                        initPP5	 = 0;
                        endPP5 = hp->myNumbStar[0];
					
		}else{

			initPP4 += hp->myNumbDark[i-1]; 
			endPP4 += hp->myNumbDark[i]; 
                     	initPP5 += hp->myNumbStar[i-1];
                        endPP5 += hp->myNumbStar[i];
				
		}		
		

		std::cout << "Qui2" << std::endl;

		std::cout << "initPP " << initPP4 << "  endPP " << endPP4 << std::endl; 
                std::cout << "initPP " << initPP5 << "  endPP " << endPP5 << std::endl;
		
		std::cout << "size of satellite4 is " << endPP4-initPP4 << std::endl;
		std::cout << "size of satellite5 is " << endPP5-initPP5 << std::endl;
		std::cout << "\n\n\n the mass of particle " << initPP4 << " is " << pp4->at(initPP4).mass  << std::endl;
		std::cout << " the mass of particle " << endPP4-1 << " is " << pp4->at(endPP4-1).mass  << "\n\n\n " << std::endl;
                std::cout << "\n\n\n the mass of particle " << initPP5 << " is " << pp5->at(initPP5).mass  << std::endl;
                std::cout << " the mass of particle " << endPP5-1 << " is " << pp5->at(endPP5-1).mass  << "\n\n\n " << std::endl;

		for(i4=initPP4;i4<endPP4;i4++){

			pp4->at(i4).flag_bound = 1;
			totMass += pp4->at(i4).mass;
			totMass_p4 += pp4->at(i4).mass;

		}

		if(hp->wantBaryons == 1){

			for(i5=initPP5;i5<endPP5;i5++){

				pp5->at(i5).flag_bound = 1;
				totMass += pp5->at(i5).mass;
				totMass_p5 += pp5->at(i5).mass;

			}

		}

		
			
	/*

	std::cout << "   " << std::endl;
	std::cout << "Calculating Fraction of Bound Particles" << std::endl;
	std::cout << "Working with total inital mass: " << totMass << std::endl;
	std::cout << "	  And with total baryon mass: " << totMass_p4 << std::endl;
	std::cout << "	    And with total dark mass: " << totMass_p5 << std::endl;
 	std::cout << "   " << std::endl;
	

	//while(err > 0.01){

		//old_fracBound = fracBound;
		fracBound = 0.;
		fracBound_p4 = 0.;
		fracBound_p5 = 0.;

		std::cout << "   " << std::endl;		
	
	          std::cout << "Debug empty infinite loop on particles 4" << std::endl;	

	if(hp->wantBaryons == 0){
		for(i4=0;i4<hp->nbody4;i4++){

			if(i4%1000 == 0){

				std::cout << i4 << std::endl;

			}

			if(pp4->at(i4).flag_bound == 1){		

				float pot = 0.;
				float kin = 1./2.*(pp4->at(i4).mass)*( pow((pp4->at(i4).vel[0] - dcms_Sat->vel[0]),2) + pow((pp4->at(i4).vel[1] - dcms_Sat->vel[1]),2) + pow((pp4->at(i4).vel[2] - dcms_Sat->vel[2]),2) );		

				for(int j4=0;j4<hp->nbody4;j4++){
	
					if((j4 != i4) and (pp4->at(j4).flag_bound ==1)){
					
						float radius_i4j4 = sqrt( pow((pp4->at(i4).pos[0] - pp4->at(j4).pos[0]),2) + pow((pp4->at(i4).pos[1] - pp4->at(j4).pos[1]),2) + pow((pp4->at(i4).pos[2] - pp4->at(j4).pos[2]),2));
						pot += -1.*(pp4->at(i4).mass)*(pp4->at(j4).mass)/(radius_i4j4);		
					
					}

				}

				for(int j5=0;j5<hp->nbody5;j5++){

					if(pp5->at(j5).flag_bound ==1){
			
						float radius_i4j5 = sqrt( pow((pp4->at(i4).pos[0] - pp5->at(j5).pos[0]),2) + pow((pp4->at(i4).pos[1] - pp5->at(j5).pos[1]),2) + pow((pp4->at(i4).pos[2] - pp5->at(j5).pos[2]),2) );
						pot += -1.*(pp4->at(i4).mass)*(pp5->at(j5).mass)/(radius_i4j5);		
				
					}
	
				}

				pp4->at(i4).potKin = pot+kin;
						
				if ((pp4->at(i4).potKin) < 0){

					fracBound += pp4->at(i4).mass;
					fracBound_p4 += pp4->at(i4).mass;
			
				}

			}   //End if condition on particle i4 flag						

		}   //End loop on particles i4					

	}

	if(hp->wantBaryons == 1){

		for(i5=0;i5<hp->nbody5;i5++){

		if(i5%10 == 0){

				std::cout << i5 << std::endl;

			}

			if(pp5->at(i5).flag_bound == 1){		

				float pot = 0.;
				float kin = 1./2.*(pp5->at(i5).mass)*( pow((pp5->at(i5).vel[0] - dcms_Sat->vel[0]),2) + pow((pp5->at(i5).vel[1] - dcms_Sat->vel[1]),2) + pow((pp5->at(i5).vel[2] - dcms_Sat->vel[2]),2) );		

				for(int j4=0;j4<hp->nbody4;j4++){
	
					if((pp4->at(j4).flag_bound ==1)){
					
						float radius_i5j4 = sqrt( pow((pp5->at(i5).pos[0] - pp4->at(j4).pos[0]),2) + pow((pp5->at(i5).pos[1] - pp4->at(j4).pos[1]),2) + pow((pp5->at(i5).pos[2] - pp4->at(j4).pos[2]),2));
						pot += -1.*(pp5->at(i5).mass)*(pp4->at(j4).mass)/(radius_i5j4);		
					
					}

				}

				for(int j5=0;j5<hp->nbody5;j5++){

					if((j5 != i5) and (pp5->at(j5).flag_bound ==1)){
			
						float radius_i5j5 = sqrt( pow((pp5->at(i5).pos[0] - pp5->at(j5).pos[0]),2) + pow((pp5->at(i5).pos[1] - pp5->at(j5).pos[1]),2) + pow((pp5->at(i5).pos[2] - pp5->at(j5).pos[2]),2) );
						pot += -1.*(pp5->at(5).mass)*(pp5->at(j5).mass)/(radius_i5j5);		
				
					}
	
				}

				pp5->at(i5).potKin = pot+kin;
						
				if ((pp5->at(i5).potKin) < 0){

					fracBound += pp5->at(i5).mass;
					fracBound_p5 += pp5->at(i5).mass;
			
				}

			}   //End if condition on particle i5 flag						

		}   //End loop on particles i5					

	}	

		//Calculate error

		//fracBound /= oldMass;	

		//err = abs(old_fracBound/totMass-fracBound/totMass);
	
		//Flag particles
		
		std::cout << "   Bound mass: " << fracBound << std::endl;
		//std::cout << "   Convergence:    " << err << std::endl;
		std::cout << "   Bound baryons (or p4): " << fracBound_p4 << std::endl;
		std::cout << "   Bound dark mass (or p5): " << fracBound_p5 << std::endl;	
		std::cout << "   Total considered mass: " << totMass << std::endl;
		std::cout << "   Total considered baryons (or p4): " << totMass_p4 << std::endl;
		std::cout << "   Total considered dark mass (or p5): " << totMass_p5 << std::endl;
		std::cout << "   " << std::endl;


		/////
		if(hp->wantBaryons == 0){

			for(i4=0;i4<hp->nbody4;i4++){
				if((pp4->at(i4).potKin) >= 0){		
					pp4->at(i4).flag_bound = 0;
				}
			}
		}

		if(hp->wantBaryons == 1){
		
			for(i5=0;i5<hp->nbody5;i5++){
				if((pp5->at(i5).potKin) >= 0){		
					pp5->at(i5).flag_bound = 0;
				}
			}
		}

		///// 

		//iter++;

	//}  //End while loop	
	
	//std::cout << "Ended iterations " << std::endl;
	std::cout << " " << std::endl;
	std::cout << "Fraction of bound mass: " << fracBound/totMass << std::endl;
 	std::cout << "Fraction of bound baryons (p4): " << fracBound_p4/totMass_p4 << std::endl;	
	std::cout << "Fraction of bound dark (p5): " << fracBound_p5/totMass_p5 << std::endl;

	*/

		std::cout << "Debug calculation sat_totmass " << std::endl;

		hp->sat_massBound.at(i) = fracBound;	
		hp->sat_massBound_p4.at(i) = fracBound_p4;
		hp->sat_massBound_p5.at(i) = fracBound_p5;
		hp->sat_totMass.at(i) = totMass;
		hp->sat_totMass_p4.at(i) = totMass_p4;
		hp->sat_totMass_p5.at(i) = totMass_p5;

		std::cout << "the tot mass for satellite " << i << " is " << hp->sat_totMass.at(i) << std::endl;

	}

}

/* Integrate MW mass within distance of satellite density center to MW density center*/

void integrate_MW_mass(struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2, std::vector<struct h5_particle> *pp3, struct header_h5 *hp){
	
	hp->tot_halo_mass = 0.;

	std::cout << "Calculating total MW mass within sphere of radius = distance MW-center <-> satellite " << std::endl;
 
	int i1, i2, i3;
	
	for(int ind_sat = 0;ind_sat<hp->sat_numb;ind_sat++){

		hp->MW_integratedMass[ind_sat] = 0.;
		float rad = std::sqrt( std::pow((dcms_Sat->at(ind_sat).pos[0] - dcms->pos[0]),2) + std::pow((dcms_Sat->at(ind_sat).pos[1] - dcms->pos[1]),2) + std::pow((dcms_Sat->at(ind_sat).pos[2] - dcms->pos[2]),2) );

	//Disk

	for(i1=0;i1<hp->nbody1;i1++){

		if(std::sqrt( std::pow((pp1->at(i1).pos[0] - dcms->pos[0]),2) + std::pow((pp1->at(i1).pos[1] - dcms->pos[1]),2) + std::pow((pp1->at(i1).pos[2] - dcms->pos[2]),2) ) < rad ){

			hp->MW_integratedMass[ind_sat] += pp1->at(i1).mass;	

		}
		if(ind_sat == 0){
			hp->tot_halo_mass += pp1->at(i1).mass;
		}
	}

	//Bulge

	for(i2=0;i2<hp->nbody2;i2++){

		if(std::sqrt( std::pow((pp2->at(i2).pos[0] - dcms->pos[0]),2) + std::pow((pp2->at(i2).pos[1] - dcms->pos[1]),2) + std::pow((pp2->at(i2).pos[2] - dcms->pos[2]),2) ) < rad ){

			hp->MW_integratedMass[ind_sat] += pp2->at(i2).mass;	

		}

		if(ind_sat == 0){
			hp->tot_halo_mass += pp2->at(i2).mass;
		}

	}

	//Halo

	for(i3=0;i3<hp->nbody3;i3++){

		if(std::sqrt( std::pow((pp3->at(i3).pos[0] - dcms->pos[0]),2) + std::pow((pp3->at(i3).pos[1] - dcms->pos[1]),2) + std::pow((pp3->at(i3).pos[2] - dcms->pos[2]),2) ) < rad ){

			hp->MW_integratedMass[ind_sat] += pp3->at(i3).mass;	

		}
	
		if(ind_sat == 0){
			hp->tot_halo_mass += pp3->at(i3).mass;
		}
	
	}

	}

}


/* Tidal radius calculation */

void tidal_radius_calculate(struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, struct header_h5 *hp){

	float G_const = 1.; //4.30e-6; //4.30e-6; //const of gravitation 
	float R200 =  242.8; //assume c.a. the r_200 for the halo of the milky way; r_scale= 25 kpc
	float c_par = hp->MW_concentration;
	float r_scale = R200/c_par;
	float M_tot = hp->tot_halo_mass;
	std::cout << "m TOT: " << M_tot << std::endl;
	float a_halo = r_scale*std::sqrt( 2.*( log(1.+c_par)-c_par/(1.+c_par) ) );
	std::cout << "a_halo: " << a_halo << std::endl;

	for (int ind_sat=0; ind_sat<hp->sat_numb;ind_sat++){

		float sat_distance = std::sqrt( std::pow((dcms_Sat->at(ind_sat).pos[0] - dcms->pos[0]),2) + std::pow((dcms_Sat->at(ind_sat).pos[1] - dcms->pos[1]),2) + std::pow((dcms_Sat->at(ind_sat).pos[2] - dcms->pos[2]),2) );
		float sat_velocity = std::sqrt( std::pow((dcms_Sat->at(ind_sat).vel[0] - dcms->vel[0]),2) + std::pow((dcms_Sat->at(ind_sat).vel[1] - dcms->vel[1]),2) + std::pow((dcms_Sat->at(ind_sat).vel[2] - dcms->vel[2]),2) );
		float reduced_ang_m_x = (dcms_Sat->at(ind_sat).pos[1] - dcms->pos[1])*(dcms_Sat->at(ind_sat).vel[2] - dcms->vel[2]) - (dcms_Sat->at(ind_sat).pos[2] - dcms->pos[2])*(dcms_Sat->at(ind_sat).vel[1] - dcms->vel[1]);
		float reduced_ang_m_y = (dcms_Sat->at(ind_sat).pos[2] - dcms->pos[2])*(dcms_Sat->at(ind_sat).vel[0] - dcms->vel[0]) - (dcms_Sat->at(ind_sat).pos[0] - dcms->pos[0])*(dcms_Sat->at(ind_sat).vel[2] - dcms->vel[2]);
		float reduced_ang_m_z = (dcms_Sat->at(ind_sat).pos[0] - dcms->pos[0])*(dcms_Sat->at(ind_sat).vel[1] - dcms->vel[1]) - (dcms_Sat->at(ind_sat).pos[1] - dcms->pos[1])*(dcms_Sat->at(ind_sat).vel[0] - dcms->vel[0]);
		//float reduced_ang_m = std::sqrt( pow(reduced_ang_m_x,2) + pow(reduced_ang_m_y,2) + pow(reduced_ang_m_z,2) );
		float M_hern = M_tot*(sat_distance*sat_distance)/((sat_distance+a_halo)*(sat_distance+a_halo));
		float M_hern_1st = 2.*M_tot*(sat_distance)*a_halo/((sat_distance+a_halo)*(sat_distance+a_halo)*(sat_distance+a_halo));
		float M_hern_2nd = 2.*M_tot*a_halo*(a_halo-2.*sat_distance)/((sat_distance+a_halo)*(sat_distance+a_halo)*(sat_distance+a_halo));
		std::cout << "M: " << M_hern << " M': " << M_hern_1st << " M'': " << M_hern_2nd << std::endl;
		//float mu = (hp->sat_totMass)*(hp->MW_integratedMass)/((hp->sat_totMass) + (hp->MW_integratedMass));
		//float mfac = (hp->sat_totMass)*(hp->MW_integratedMass);

		std::cout << "  Satellite has dcms : \n" << "  " << dcms_Sat->at(ind_sat).pos[0] << " " <<  dcms_Sat->at(ind_sat).pos[1] << " "	<< 
		dcms_Sat->at(ind_sat).pos[2] << " "	<< dcms_Sat->at(ind_sat).vel[0] << " "	<< dcms_Sat->at(ind_sat).vel[1] << " "	<< 
		dcms_Sat->at(ind_sat).vel[2] << std::endl;
	
		float v_ang_2 = (sat_velocity*sat_velocity)/(sat_distance*sat_distance);
		std::cout << "Angular speed " << v_ang_2 << std::endl;
		float myLog = std::log(sat_distance/R200);
		//float d2Phi_dR2 = G_const*M_hern/(sat_distance*sat_distance*sat_distance)*(-3.+2.*myLog) + 2.*G_const*M_hern_1st/(sat_distance*sat_distance)*(1.-myLog) + G_const*M_hern_2nd/(sat_distance)*myLog;
		//float d2Phi_dR2 = - 
		float d2Phi_dR2 = -2.*G_const*M_tot/(std::pow((sat_distance+a_halo),3)); // For hernquist profile
		std::cout << "second derivative of potential: " << d2Phi_dR2 << std::endl;
		//float V_circ = std::sqrt(G_const*(hp->MW_integratedMass)/sat_distance);
		//float tot_En = -1.*G*mfac/sat_distance + 1./2.*mu*(sat_velocity*sat_velocity); 
		//float tot_En = 1./2.*(sat_velocity*sat_velocity) + V_circ*V_circ*std::log(sat_distance/R0);
		//float ecc = std::sqrt( 1. + 2.*tot_En*pow(reduced_ang_m,2)*pow(mu,2)/(mu*mfac*mfac*G_const*G_const) );
		//float M_grav_param = G_const*(hp->MW_integratedMass);
		//float a_ax = pow(reduced_ang_m,2)/(M_grav_param*(1 - ecc*ecc));
	
		//float r_p = a_ax*(1. - ecc);
		//float r_a = a_ax*(1. + ecc);
	 
		std::cout << "Calculating tidal radius for satellite " << ind_sat+1  << " ..." << std::endl;
		std::cout << "difference: " << v_ang_2 - d2Phi_dR2 << std::endl;
		std::cout << "The satellite mass used here is: " << hp->sat_totMass[ind_sat] << std::endl;
		hp->sat_tidRadius[ind_sat] = std::pow(((G_const*(hp->sat_totMass[ind_sat]))/(v_ang_2 - d2Phi_dR2)),1./3.);
		//hp->sat_tidRadius = pow(((hp->sat_massBound)/(3.*(hp->MW_integratedMass))),(1./3.))*sqrt( pow((dcms_Sat->pos[0] - dcms->pos[0]),2) + pow((dcms_Sat->pos[1] - dcms->pos[1]),2) + pow((dcms_Sat->pos[2] - dcms->pos[2]),2) );	

		std::cout << "Tidal radius is " << hp->sat_tidRadius[ind_sat] << std::endl; 
	
	}
	
}

/* Calculate tidal mass fraction */

void calculate_TidalMass_Fraction(std::vector<struct mass_center> *dcms_Sat, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, std::string filenumb){ 	

	//std::ofstream myFileDebrisStar((hp->myHDF5Root+"_040_starsDebris.txt").c_str());
	//std::ofstream myFileDebrisDark((hp->myHDF5Root+"_040_darkDebris.txt").c_str());	

	//myFileDebrisDark << 0 << std::endl;
	//myFileDebrisDark << 0 << std::endl;
	//myFileDebrisDark << std::scientific << std::setprecision(6) << 0 << std::endl;

	//myFileDebrisStar << 0 << std::endl;
	//myFileDebrisStar << 0 << std::endl;
	//myFileDebrisStar << std::scientific << std::setprecision(6) << 0 << std::endl;

	int initPart4, initPart5;
        int finPart4, finPart5;

	float x_vels_mean_p4[hp->sat_numb];	
	float x_vsigma_p4[hp->sat_numb];
	float y_vels_mean_p4[hp->sat_numb];	
	float y_vsigma_p4[hp->sat_numb];
	float z_vels_mean_p4[hp->sat_numb];	
	float z_vsigma_p4[hp->sat_numb];
	float num_part_p4[hp->sat_numb];
	float mass_vels_p4[hp->sat_numb];

	float x_vels_mean_p5[hp->sat_numb];	
	float x_vsigma_p5[hp->sat_numb];
	float y_vels_mean_p5[hp->sat_numb];	
	float y_vsigma_p5[hp->sat_numb];
	float z_vels_mean_p5[hp->sat_numb];	
	float z_vsigma_p5[hp->sat_numb];
	float num_part_p5[hp->sat_numb];
	float mass_vels_p5[hp->sat_numb];

	for (int ind_sat=0;ind_sat<hp->sat_numb;ind_sat++){

		std::cout << "Calculating tidal mass for satellite " << ind_sat+1 << " within tidal radius " <<  hp->sat_tidRadius[ind_sat] << std::endl; 

		int i4, i5;
		hp->sat_tidalMass[ind_sat] = 0.;
		hp->sat_tidalMass_p4[ind_sat] = 0.;
		hp->sat_tidalMass_p5[ind_sat] = 0.;

		if(ind_sat == 0){

			initPart4 = 0;
			finPart4 = hp->myNumbDark[0];
			initPart5 = 0;
			finPart5 = hp->myNumbStar[0];

		}else{

			initPart4 += hp->myNumbDark[ind_sat-1];
			finPart4 += hp->myNumbDark[ind_sat];
                      	initPart5 += hp->myNumbStar[ind_sat-1];
                        finPart5 += hp->myNumbStar[ind_sat];


		}

		x_vels_mean_p4[ind_sat] = 0.;	
		mass_vels_p4[ind_sat] = 0.;
		x_vsigma_p4[ind_sat] = 0.;
		y_vels_mean_p4[ind_sat] = 0.;	
		y_vsigma_p4[ind_sat] = 0.;
		z_vels_mean_p4[ind_sat] = 0.;	
		z_vsigma_p4[ind_sat] = 0.;
		num_part_p4[ind_sat] = 0;

		x_vels_mean_p5[ind_sat] = 0.;	
		mass_vels_p5[ind_sat] = 0.;
		x_vsigma_p5[ind_sat] = 0.;
		y_vels_mean_p5[ind_sat] = 0.;	
		y_vsigma_p5[ind_sat] = 0.;
		z_vels_mean_p5[ind_sat] = 0.;	
		z_vsigma_p5[ind_sat] = 0.;
		num_part_p5[ind_sat] = 0;

		for(i4=initPart4;i4<finPart4;i4++){

			if(std::sqrt(std::pow((pp4->at(i4).pos[0] - dcms_Sat->at(ind_sat).pos[0]),2) + std::pow((pp4->at(i4).pos[1] - dcms_Sat->at(ind_sat).pos[1]),2) + std::pow((pp4->at(i4).pos[2] - dcms_Sat->at(ind_sat).pos[2]),2) ) < hp->sat_tidRadius[ind_sat]){

				hp->sat_tidalMass[ind_sat] += pp4->at(i4).mass;
				hp->sat_tidalMass_p4[ind_sat] += pp4->at(i4).mass;	
				pp4->at(i4).flag_bound = 1;
				//calculate then  the velocities for the sigma calculation for DM
				x_vels_mean_p4[ind_sat] += pp4->at(i4).mass*(pp4->at(i4).vel[0]-dcms_Sat->at(ind_sat).vel[0]);
				y_vels_mean_p4[ind_sat] += pp4->at(i4).mass*(pp4->at(i4).vel[1]-dcms_Sat->at(ind_sat).vel[1]);		
				z_vels_mean_p4[ind_sat] += pp4->at(i4).mass*(pp4->at(i4).vel[2]-dcms_Sat->at(ind_sat).vel[2]);
				if(pp4->at(i4).mass != 0.){
					num_part_p4[ind_sat] += 1;
				}
	
				mass_vels_p4[ind_sat] += pp4->at(i4).mass;				
						
				
 			}else{

				pp4->at(i4).flag_bound = 0;
				hp->sat_outTidMass[ind_sat] += pp4->at(i4).mass;		
				hp->sat_outTidMass_p4[ind_sat] += pp4->at(i4).mass;			

				//myFileDebrisDark << std::scientific << std::setprecision(6) << pp4->at(i4).idf << "\t" << pp4->at(i4).mass << "\t"  << pp4->at(i4).pos[0] << "\t"  << pp4->at(i4).pos[1] << "\t" << pp4->at(i4).pos[2] << "\t" << pp4->at(i4).vel[0] << "\t" << pp4->at(i4).vel[1] << "\t" << pp4->at(i4).vel[2] << std::endl;

			}

		}
		
		if(mass_vels_p4[ind_sat] != 0.){
			x_vels_mean_p4[ind_sat] /= mass_vels_p4[ind_sat];
			y_vels_mean_p4[ind_sat] /= mass_vels_p4[ind_sat];
			z_vels_mean_p4[ind_sat] /= mass_vels_p4[ind_sat];
		}else{
			x_vels_mean_p4[ind_sat] = 0.;
			y_vels_mean_p4[ind_sat] = 0.;
			z_vels_mean_p4[ind_sat] = 0.;
		}

		for(i4=initPart4;i4<finPart4;i4++){
			if(pp4->at(i4).flag_bound == 1){
				x_vsigma_p4[ind_sat] += pp4->at(i4).mass*(pp4->at(i4).vel[0]-dcms_Sat->at(ind_sat).vel[0]-x_vels_mean_p4[ind_sat])*(pp4->at(i4).vel[0]-dcms_Sat->at(ind_sat).vel[0]-x_vels_mean_p4[ind_sat]);
				y_vsigma_p4[ind_sat] += pp4->at(i4).mass*(pp4->at(i4).vel[1]-dcms_Sat->at(ind_sat).vel[1]-y_vels_mean_p4[ind_sat])*(pp4->at(i4).vel[1]-dcms_Sat->at(ind_sat).vel[1]-y_vels_mean_p4[ind_sat]);
				z_vsigma_p4[ind_sat] += pp4->at(i4).mass*(pp4->at(i4).vel[2]-dcms_Sat->at(ind_sat).vel[2]-z_vels_mean_p4[ind_sat])*(pp4->at(i4).vel[2]-dcms_Sat->at(ind_sat).vel[2]-z_vels_mean_p4[ind_sat]);
			}
		}
		x_vsigma_p4[ind_sat] = std::sqrt( x_vsigma_p4[ind_sat]/( ((num_part_p4[ind_sat]-1)/(num_part_p4[ind_sat]))*mass_vels_p4[ind_sat] ) );
		y_vsigma_p4[ind_sat] = std::sqrt( y_vsigma_p4[ind_sat]/( ((num_part_p4[ind_sat]-1)/(num_part_p4[ind_sat]))*mass_vels_p4[ind_sat] ) );
		z_vsigma_p4[ind_sat] = std::sqrt( z_vsigma_p4[ind_sat]/( ((num_part_p4[ind_sat]-1)/(num_part_p4[ind_sat]))*mass_vels_p4[ind_sat] ) );

		

		std::cout << "for satellite " << ind_sat+1 << "there are " << hp->sat_tidalMass[ind_sat] << " unit masses inside and " << hp->sat_outTidMass[ind_sat] << " outside" << std::endl; 


		if(hp->wantBaryons == 1){	

			for(i5=initPart5;i5<finPart5;i5++){

				if(std::sqrt( std::pow((pp5->at(i5).pos[0] - dcms_Sat->at(ind_sat).pos[0]),2) + std::pow((pp5->at(i5).pos[1] - dcms_Sat->at(ind_sat).pos[1]),2) + std::pow((pp5->at(i5).pos[2] - dcms_Sat->at(ind_sat).pos[2]),2) ) < hp->sat_tidRadius[ind_sat]  ){

					hp->sat_tidalMass[ind_sat] += pp5->at(i5).mass;	
					hp->sat_tidalMass_p5[ind_sat] += pp5->at(i5).mass;
					pp5->at(i5).flag_bound = 1;
					x_vels_mean_p5[ind_sat] += pp5->at(i5).mass*(pp5->at(i5).vel[0]-dcms_Sat->at(ind_sat).vel[0]);
					y_vels_mean_p5[ind_sat] += pp5->at(i5).mass*(pp5->at(i5).vel[1]-dcms_Sat->at(ind_sat).vel[1]);		
					z_vels_mean_p5[ind_sat] += pp5->at(i5).mass*(pp5->at(i5).vel[2]-dcms_Sat->at(ind_sat).vel[2]);
					if(pp5->at(i5).mass != 0.){
						num_part_p5[ind_sat] += 1;
					}
					mass_vels_p5[ind_sat] += pp4->at(i5).mass;				


				}else{

					pp5->at(i5).flag_bound = 0;
					hp->sat_outTidMass[ind_sat] += pp5->at(i5).mass;		
					hp->sat_outTidMass_p5[ind_sat] += pp5->at(i5).mass;			

					//myFileDebrisStar << std::scientific << std::setprecision(6) << pp5->at(i5).idf << "\t" << pp5->at(i5).mass << "\t"  << pp5->at(i5).pos[0] << "\t"  << pp5->at(i5).pos[1] << "\t" << pp5->at(i5).pos[2] << "\t" << pp5->at(i5).vel[0] << "\t" << pp5->at(i5).vel[1] << "\t" << pp5->at(i5).vel[2] << std::endl;

				}

			}//for loop on p5

			if(mass_vels_p5[ind_sat] != 0.){
				x_vels_mean_p5[ind_sat] /= mass_vels_p5[ind_sat];
				y_vels_mean_p5[ind_sat] /= mass_vels_p5[ind_sat];
				z_vels_mean_p5[ind_sat] /= mass_vels_p5[ind_sat];
			}else{
				x_vels_mean_p5[ind_sat] = 0.;
				y_vels_mean_p5[ind_sat] = 0.;
				z_vels_mean_p5[ind_sat] = 0.;
			}

			for(i5=initPart5;i5<finPart5;i5++){
				if(pp5->at(i5).flag_bound == 1){
					x_vsigma_p5[ind_sat] += pp5->at(i5).mass*(pp5->at(i5).vel[0]-dcms_Sat->at(ind_sat).vel[0]-x_vels_mean_p5[ind_sat])*(pp5->at(i5).vel[0]-dcms_Sat->at(ind_sat).vel[0]-x_vels_mean_p5[ind_sat]);
					y_vsigma_p5[ind_sat] += pp5->at(i5).mass*(pp5->at(i5).vel[1]-dcms_Sat->at(ind_sat).vel[1]-y_vels_mean_p5[ind_sat])*(pp5->at(i5).vel[1]-dcms_Sat->at(ind_sat).vel[1]-y_vels_mean_p5[ind_sat]);
					z_vsigma_p5[ind_sat] += pp5->at(i5).mass*(pp5->at(i5).vel[2]-dcms_Sat->at(ind_sat).vel[2]-z_vels_mean_p5[ind_sat])*(pp5->at(i5).vel[2]-dcms_Sat->at(ind_sat).vel[2]-z_vels_mean_p5[ind_sat]);
				}
			}
			x_vsigma_p5[ind_sat] = std::sqrt( x_vsigma_p5[ind_sat]/( ((num_part_p5[ind_sat]-1)/(num_part_p5[ind_sat]))*mass_vels_p5[ind_sat] ) );
			y_vsigma_p5[ind_sat] = std::sqrt( y_vsigma_p5[ind_sat]/( ((num_part_p5[ind_sat]-1)/(num_part_p5[ind_sat]))*mass_vels_p5[ind_sat] ) );
			z_vsigma_p5[ind_sat] = std::sqrt( z_vsigma_p5[ind_sat]/( ((num_part_p5[ind_sat]-1)/(num_part_p5[ind_sat]))*mass_vels_p5[ind_sat] ) );

		}//want baryons
	
	}//satellite loop

	std::ofstream myFileSigmasPrint((hp->myHDF5Root+"_sigmas_sat_within_TidalRad_"+filenumb+".txt").c_str());

	for(int isat=0;isat<hp->sat_numb;isat++){

		myFileSigmasPrint << std::scientific << std::setprecision(6) << x_vsigma_p4[isat] << "\t" << y_vsigma_p4[isat] << "\t" << z_vsigma_p4[isat] << "\t" << num_part_p4[isat] << "\t"  << mass_vels_p4[isat] << "\t" << x_vsigma_p5[isat] << "\t"  << y_vsigma_p5[isat] << "\t"  << z_vsigma_p5[isat] << "\t" << num_part_p5[isat] <<  "\t"  << mass_vels_p4[isat] << std::endl; 

	}

	myFileSigmasPrint.close();

	//myFileDebrisStar.close();	
	//myFileDebrisDark.close();	

}

/* Function to relocate particle coordinates with respect to density center */

void density_center_relocate(struct mass_center *dcms, std::vector<struct h5_particle> *pp, int *nPart){

	int iP;
	
	for (iP =0;iP<*nPart;iP++){

		for (int ii=0;ii<3;ii++){

			pp->at(iP).pos[ii] -= dcms->pos[ii];
			pp->at(iP).vel[ii] -= dcms->vel[ii];

		}

	}

}

void calc_Energy_Lang_DcmsSat(struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp){
	
  int i4, i5, minIndex;
	int initPart4, initPart5;
        int finPart4, finPart5;
	float rad, minDist;
	std::vector<struct mass_center> sat_cms;
	std::vector<float> cms_mass = std::vector<float> (0);
	std::vector<float> aver_potential = std::vector<float> (0);
	std::vector<float> en_kin = std::vector<float> (0);
	std::vector<float> en_tot = std::vector<float> (0);
	std::vector<float> ang_mom_z = std::vector<float> (0);

	sat_cms.resize(hp->sat_numb);
	cms_mass.resize(hp->sat_numb);
	aver_potential.resize(hp->sat_numb);
	en_kin.resize(hp->sat_numb);
	en_tot.resize(hp->sat_numb);
	ang_mom_z.resize(hp->sat_numb);

	for (int ind_sat=0;ind_sat<hp->sat_numb;ind_sat++){

		cms_mass[ind_sat] = 0.;	
		aver_potential[ind_sat] = 0.;
		en_kin[ind_sat] = 0.;
		ang_mom_z[ind_sat] = 0.;	

	}

	for (int ind_sat=0;ind_sat<hp->sat_numb;ind_sat++){

		if(ind_sat == 0){

			initPart4 = 0;
			finPart4 = hp->myNumbDark[0];
			initPart5 = 0;
			finPart5 = hp->myNumbStar[0];

		}else{
	
			initPart4 += hp->myNumbDark[ind_sat-1];
			finPart4 += hp->myNumbDark[ind_sat];
	                initPart5 += hp->myNumbStar[ind_sat-1];
	                finPart5 += hp->myNumbStar[ind_sat];

		}
	
		aver_potential[ind_sat] = 0.;
		minDist = 1000.;

		for(i4=initPart4;i4<finPart4;i4++){

			rad = std::sqrt(std::pow(pp4->at(i4).pos[0]-dcms_Sat->at(ind_sat).pos[0],2) + std::pow(pp4->at(i4).pos[1]-dcms_Sat->at(ind_sat).pos[1],2) + std::pow(pp4->at(i4).pos[2]-dcms_Sat->at(ind_sat).pos[2],2) );

			if( (pp4->at(i4).flag_bound == 1) and (rad < hp->sat_tidRadius[ind_sat]/15.) ){
			
			        //aver_potential[ind_sat] += pp4->at(i4).phi*pp4->at(i4).mass;
				cms_mass[ind_sat] += pp4->at(i4).mass;
				
				for (int j=0;j<3;j++){	
				
					sat_cms[ind_sat].pos[j] += pp4->at(i4).pos[j]*pp4->at(i4).mass;
					sat_cms[ind_sat].vel[j] += pp4->at(i4).vel[j]*pp4->at(i4).mass;

				}
			
				if(rad < minDist){

                                        minDist = rad;  
					aver_potential[ind_sat] = pp4->at(i4).phi*pp4->at(i4).mass;
				}	
				  
			}

		}

		for(i5=initPart5;i5<finPart5;i5++){

			rad = std::sqrt(std::pow(pp5->at(i5).pos[0]-dcms_Sat->at(ind_sat).pos[0],2) + std::pow(pp5->at(i5).pos[1]-dcms_Sat->at(ind_sat).pos[1],2) + std::pow(pp5->at(i5).pos[2]-dcms_Sat->at(ind_sat).pos[2],2) );

			if( (pp5->at(i5).flag_bound == 1) and (rad < hp->sat_tidRadius[ind_sat]/15.) ){
			
			        //aver_potential[ind_sat] += pp5->at(i5).phi*pp5->at(i5).mass;
				cms_mass[ind_sat] += pp5->at(i5).mass;

				for (int j=0;j<3;j++){	
				
					sat_cms[ind_sat].pos[j] += pp5->at(i5).pos[j]*pp5->at(i5).mass;
					sat_cms[ind_sat].vel[j] += pp5->at(i5).vel[j]*pp5->at(i5).mass;

				}

				if(rad < minDist){

				  minDist = rad;
				  aver_potential[ind_sat] = pp5->at(i5).phi*pp5->at(i5).mass;

				}
			
			}

		}

		std::cout << "Average potential is: " << aver_potential[ind_sat] << " and mass \t" << cms_mass[ind_sat] << " and minDist is " << minDist << std::endl;	
			
		//aver_potential[ind_sat] /= cms_mass[ind_sat];

		for(int j=0;j<3;j++){

		  sat_cms[ind_sat].pos[j] /= cms_mass[ind_sat];
		  sat_cms[ind_sat].vel[j] /= cms_mass[ind_sat];

		}    

		en_kin[ind_sat] = 1./2.*( std::pow( sat_cms.at(ind_sat).vel[0],2 ) + std::pow( sat_cms.at(ind_sat).vel[1],2 ) + std::pow( sat_cms.at(ind_sat).vel[2],2 ) );   
		ang_mom_z[ind_sat] = (sat_cms.at(ind_sat).pos[0])*(sat_cms.at(ind_sat).vel[1]) - (sat_cms.at(ind_sat).pos[1])*(sat_cms.at(ind_sat).vel[0]);
		en_tot[ind_sat] = aver_potential[ind_sat] + en_kin[ind_sat];		
			
	}	

	
	std::ofstream EnAng_DcmsSat_OutFile((hp->myHDF5Root+"_enAng_DcmsSat.txt").c_str());	
	
	for (int ind_sat = 0 ;ind_sat < hp->sat_numb; ind_sat++){

		EnAng_DcmsSat_OutFile << std::scientific << std::setprecision(6) << en_kin[ind_sat] << "\t" << aver_potential[ind_sat] << "\t" << en_tot[ind_sat] << "\t" << ang_mom_z[ind_sat] << std::endl; 

		std::cout << "     and tot energy is " << en_tot[ind_sat] << std::endl;
								
	}

	EnAng_DcmsSat_OutFile.close();

}


/* Calculate angular momentum of Milky Way disk */

void angular_momentum_MW_calculate(struct mass_center *dcms, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3,struct header_h5 *hp){

	std::cout << "Calculate angular momentum for Milky Way" << std::endl;

	int i1, i2, i3;

	for(i1=0;i1<hp->nbody1;i1++){	

		dcms->angM[0] += pp1->at(i1).mass*( (pp1->at(i1).pos[1])*(pp1->at(i1).vel[2]) - (pp1->at(i1).pos[2])*(pp1->at(i1).vel[1]) );
		dcms->angM[1] += pp1->at(i1).mass*( (pp1->at(i1).pos[2])*(pp1->at(i1).vel[0]) - (pp1->at(i1).pos[0])*(pp1->at(i1).vel[2]) );
		dcms->angM[2] += pp1->at(i1).mass*( (pp1->at(i1).pos[0])*(pp1->at(i1).vel[1]) - (pp1->at(i1).pos[1])*(pp1->at(i1).vel[0]) );

	}

	/*

	for(i2=0;i2<hp->nbody2;i2++){

		dcms->angM[0] += pp2->at(i2).mass*( (pp2->at(i2).pos[1])*(pp2->at(i2).vel[2]) - (pp2->at(i2).pos[2])*(pp2->at(i2).vel[1]) );
		dcms->angM[1] += pp2->at(i2).mass*( (pp2->at(i2).pos[2])*(pp2->at(i2).vel[0]) - (pp2->at(i2).pos[0])*(pp2->at(i2).vel[2]) );
		dcms->angM[2] += pp2->at(i2).mass*( (pp2->at(i2).pos[0])*(pp2->at(i2).vel[1]) - (pp2->at(i2).pos[1])*(pp2->at(i2).vel[0]) );

	}

	for(i3=0;i3<hp->nbody3;i3++){

		dcms->angM[0] += pp3->at(i3).mass*( (pp3->at(i3).pos[1])*(pp3->at(i3).vel[2]) - (pp3->at(i3).pos[2])*(pp3->at(i3).vel[1]) );
		dcms->angM[1] += pp3->at(i3).mass*( (pp3->at(i3).pos[2])*(pp3->at(i3).vel[0]) - (pp3->at(i3).pos[0])*(pp3->at(i3).vel[2]) );
		dcms->angM[2] += pp3->at(i3).mass*( (pp3->at(i3).pos[0])*(pp3->at(i3).vel[1]) - (pp3->at(i3).pos[1])*(pp3->at(i3).vel[0]) );

	}

	*/

	std::cout << "Final angular momentum: " << dcms->angM[0] << " " << dcms->angM[1] << " " << dcms->angM[2] << std::endl; 

}




/* Correct particle positions and velocities for angular momentum vector */

void angular_momentum_corrections(struct mass_center *dcms, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5,struct header_h5 *hp){

	float lonNodeCont, diskIncCont;
	float Lxy, LL;
	double pi = std::atan(1.0)*4;

	Lxy = std::sqrt(std::pow(dcms->angM[0],2)+std::pow(dcms->angM[1],2));	
	LL = std::sqrt(std::pow(dcms->angM[0],2)+std::pow(dcms->angM[1],2)+std::pow(dcms->angM[2],2));
	diskIncCont = std::acos((dcms->angM[2])/LL);
	if (dcms->angM[0] >= 0)
	{
		lonNodeCont = std::acos(-(dcms->angM[1])/Lxy);
	}
	else
	{
		lonNodeCont = 2.*pi - std::acos(-(dcms->angM[1])/Lxy);
	}

	std::cout << "Angular momentum angles: " << lonNodeCont << " and " << diskIncCont << " " << std::endl;

	// Call for coordinate rotations

	std::cout << "Rotating coordinates 1 " << std::endl; 
	coordinates_rotate(&lonNodeCont,&diskIncCont,pp1,&(hp->nbody1));		
	std::cout << "Rotating coordinates 2 " << std::endl; 
	coordinates_rotate(&lonNodeCont,&diskIncCont,pp2,&(hp->nbody2));		
	std::cout << "Rotating coordinates 3 " << std::endl; 
	coordinates_rotate(&lonNodeCont,&diskIncCont,pp3,&(hp->nbody3));		
	std::cout << "Rotating coordinates 4 " << std::endl; 
	coordinates_rotate(&lonNodeCont,&diskIncCont,pp4,&(hp->nbody4));
	if(hp->wantBaryons == 1){
		std::cout << "Rotating coordinates 5 " << std::endl; 		
		coordinates_rotate(&lonNodeCont,&diskIncCont,pp5,&(hp->nbody5));		
	}

}


/* Rotation of coordinates for angular momentum inclination of disk */

void coordinates_rotate(float *lonNodeCont, float *diskIncCont, std::vector<struct h5_particle> *pp, int *nPart){

	int iP;
	float xP, yP, zP, vxP, vyP, vzP; //Provvisory coordinates after first correction
	float xS, yS, zS, vxS, vyS, vzS;//Provvisory coordinates after second correction
	float xT, yT, zT, vxT, vyT, vzT; //Provvisory coordinates after third correction

	for(iP=0;iP<(*nPart);iP++){

	        //First rotation angle lonNodeCont around Z axis
	        xP = pp->at(iP).pos[0]*std::cos(-(*lonNodeCont)) - pp->at(iP).pos[1]*std::sin(-(*lonNodeCont));	
	        yP = pp->at(iP).pos[0]*std::sin(-(*lonNodeCont)) + pp->at(iP).pos[1]*std::cos(-(*lonNodeCont));
        	zP = pp->at(iP).pos[2];
        	vxP = pp->at(iP).vel[0]*std::cos(-(*lonNodeCont)) - pp->at(iP).vel[1]*std::sin(-(*lonNodeCont));
        	vyP = pp->at(iP).vel[0]*std::sin(-(*lonNodeCont)) + pp->at(iP).vel[1]*std::cos(-(*lonNodeCont));
        	vzP = pp->at(iP).vel[2];
        	//Second rotation angle diskIncCont around X axis
        	xS = xP;
        	yS = yP*std::cos(-(*diskIncCont)) - zP*std::sin(-(*diskIncCont));
        	zS = yP*std::sin(-(*diskIncCont)) + zP*std::cos(-(*diskIncCont));
        	vxS = vxP;
        	vyS = vyP*std::cos(-(*diskIncCont)) - vzP*std::sin(-(*diskIncCont));
        	vzS = vyP*std::sin(-(*diskIncCont)) + vzP*std::cos(-(*diskIncCont));
        	//Third rotation angle -lonNodeCont around Z axis
        	xT = xS*std::cos(*lonNodeCont) - yS*std::sin(*lonNodeCont);
        	yT = xS*std::sin(*lonNodeCont) + yS*std::cos(*lonNodeCont);
        	zT = zS;
        	vxT = vxS*std::cos(*lonNodeCont) - vyS*std::sin(*lonNodeCont);
        	vyT = vxS*std::sin(*lonNodeCont) + vyS*std::cos(*lonNodeCont);
        	vzT = vzS;
        	pp->at(iP).pos[0] = xT;
        	pp->at(iP).pos[1] = yT;
        	pp->at(iP).pos[2] = zT;
        	pp->at(iP).vel[0] = vxT;
        	pp->at(iP).vel[1] = vyT;
        	pp->at(iP).vel[2] = vzT;   

	}

}

void mass_distribution_flag(std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, int *thisFile){

	/*
	
	Using these arrays:
	bulge_vels_p4[10];
        disk_Vtan_p4[40];
        disk_Vzed_p4[40];
        disk_Zed_p4[40];
        bulge_vels_p5[10];
        disk_Vtan_p5[40];
        disk_Vzed_p5[40];
        disk_Zed_p5[40];

	*/

	float disk_other_mass_p4[40];
	float disk_other_mass_p5[40];
	
	for(int i=0;i<40;i++){
	
		disk_other_mass_p4[i] = 0.;
		disk_other_mass_p5[i] = 0.;		
	
	}	

	int myNbody4 = 0;
	int myNbody5 = 0; 

	std::cout << "Doing zone flagging on particles..." << std::endl;

	for(int aho=0;aho<hp->sat_numb;aho++){

		myNbody4 += hp->myNumbDark[aho];
		myNbody5 += hp->myNumbStar[aho];
		
	}

	float myBar_inDisk = 0.;
	float myBar_inBulge = 0.;
	float myBar_inHalo = 0.;
	float myTot_Bar = 0.;	

	float myDark_inDisk = 0.;
	float myDark_inBulge = 0.;
	float myDark_inHalo = 0.;
	float myTot_Dark = 0.;	

	//hp->nbody4

	/*
	std::cout << "\n NBDOY4 " << hp->nbody4 << "\n\n\n\n\n\n\n\n\n\n\n " << std::endl; 
	std::cout <<
	std::cout << "\n NBDOY5 " << hp->nbody5 << "\n\n\n\n\n\n\n\n\n\n\n " << std::endl;
	*/

	std::cout << "\nParticle size 4: " << pp4->size() << std::endl;
	std::cout << "\nParticle nbody4: " << hp->nbody4 << std::endl;
        std::cout << "\nParticle size 5: " << pp5->size() << std::endl;
        std::cout << "\nParticle nbody5: " << hp->nbody5 << std::endl;


	float potDisk = 0.;
	float potBulge = 0.;
	float potHalo = 0.;
	float kin;
	float cosTheta, sinTheta;
	float totDisk = 0.;
	float totBulge = 0.;
	float totHalo = 0.;	
	float rad, vel, rad3D, zed, vZed, vTan, vRad;
	int iP;
	int j;
	int iBin, iBin_b, iBin_h;
	float debugTotmass;

	std::cout << "	 Doing zone flagging on particle 4" << std::endl;
	
	for(iP = 0; iP < hp->nbody4; iP++){
		debugTotmass += pp4->at(iP).mass;
		if(pp4->at(iP).flag_bound == 0){
		  rad3D = std::sqrt( std::pow(pp4->at(iP).pos[0],2) + std::pow(pp4->at(iP).pos[1],2) + std::pow(pp4->at(iP).pos[2],2));	

			rad = std::sqrt( std::pow(pp4->at(iP).pos[0],2) + std::pow(pp4->at(iP).pos[1],2));
			cosTheta = (pp4->at(iP).pos[0])/rad;	
			sinTheta = (pp4->at(iP).pos[1])/rad;		
			vRad = (pp4->at(iP).vel[0])*cosTheta + (pp4->at(iP).vel[1])*sinTheta;			
			vTan = -1.*(pp4->at(iP).vel[0])*sinTheta + (pp4->at(iP).vel[1])*cosTheta;
			zed = pp4->at(iP).pos[2];
			vZed = pp4->at(iP).vel[2];
			myTot_Dark += pp4->at(iP).mass;
			//if ((rad >= 20.)){
			if ((rad >= 20.) or (std::abs(zed) >= 2.)){

				pp4->at(iP).flag_zone = 3; //halo particle				
				hp->halo_Mass_p4_tot += pp4->at(iP).mass;
				myDark_inHalo += pp4->at(iP).mass;
				//std::cout << "I found one!!!!!! wow \n" << std::endl;	

			}else{

				if ((rad3D >= 1.0)){

					pp4->at(iP).flag_zone = 1; //disk particle
					myDark_inDisk += pp4->at(iP).mass;	
		
					/*

					for(iBin = 0;iBin<40;iBin++){
						if((rad >= (iBin*0.5)) and (rad < ((iBin+1)*0.5))){
						  //if( ( std::abs(zed-(hp->disk_Zed_p1[iBin])) < 3.*(hp->disk_SigZed_p1[iBin]) ) and ( std::abs(vZed-(hp->disk_Vzed_p1[iBin])) < 3.*(hp->disk_SigVzed_p1[iBin]) ) ){
							//if( abs(vZed-(hp->disk_Vzed_p1[iBin])) < 3.*(hp->disk_SigVzed_p1[iBin]) ){
								pp4->at(iP).flag_zone = 1; //disk particle
								hp->disk_Mass_p4[iBin] += pp4->at(iP).mass;
								disk_other_mass_p4[iBin]+= pp4->at(iP).mass;
								hp->disk_Zed_p4[iBin] += pp4->at(iP).mass*(pp4->at(iP).pos[2]);
								hp->disk_Vzed_p4[iBin] += pp4->at(iP).mass*(pp4->at(iP).vel[2]);
								hp->disk_Vtan_p4[iBin]	 += pp4->at(iP).mass*vTan;
								if(iBin>0){

									hp->disk_Zed_p4[iBin-1] += pp4->at(iP).mass*(pp4->at(iP).pos[2]);
									hp->disk_Vzed_p4[iBin-1] += pp4->at(iP).mass*(pp4->at(iP).vel[2]);
									hp->disk_Vtan_p4[iBin-1] += pp4->at(iP).mass*vTan;
									disk_other_mass_p4[iBin-1] += pp4->at(iP).mass;
								}		
								if(iBin<39){

									hp->disk_Zed_p4[iBin+1] += pp4->at(iP).mass*(pp4->at(iP).pos[2]);
									hp->disk_Vzed_p4[iBin+1] += pp4->at(iP).mass*(pp4->at(iP).vel[2]);
									hp->disk_Vtan_p4[iBin+1] += pp4->at(iP).mass*vTan;
									disk_other_mass_p4[iBin+1]+= pp4->at(iP).mass;

								}

								break; //added later on 6 Nov 2017
								
							

							}else{
			
								pp4->at(iP).flag_zone = 3; //halo particle
								hp->halo_Mass_p4_tot += pp4->at(iP).mass;
								break;
										
							}

								

					        } // end if((rad >= (iBin*0.5)) and (rad < ((iBin+1)*0.5))){  	 

					} //for iBin...

					*/

				}else{

					pp4->at(iP).flag_zone = 2; //bulge particle
					myDark_inBulge += pp4->at(iP).mass;

					/*

					vel = std::sqrt( std::pow(pp4->at(iP).vel[0],2) + std::pow(pp4->at(iP).vel[1],2) + std::pow(pp4->at(iP).vel[2],2) );
					for(iBin_b = 0;iBin_b<10;iBin_b++){

						if((rad3D >= (iBin_b*0.1)) and (rad3D < ((iBin_b+1)*0.1))){
						
							if(vel<(3.*hp->bulge_Sigma_p2[iBin_b])){			
	
								pp4->at(iP).flag_zone = 2; //bulge particle
								hp->bulge_Mass_p4[iBin_b] += pp4->at(iP).mass;	
								//pp4->at(iP).flag_zone = 2; //bulge particle				
								break;		
							
							}else{

												
								pp4->at(iP).flag_zone = 3; //halo particle
								hp->halo_Mass_p4_tot += pp4->at(iP).mass;
								break;						
							}

							 	
						}
				
					}

					*/

				} //end else on rad3d
									
			} //end else on cylinder/on radius xy

		} // end if on p4 not bound

	}//end loop on p4 particles

  
		//if(pp4->at(iP).flag_bound == 0){ //so if the particle is out of satellite

		//	kin = 1./2.*(pp4->at(iP).mass)*( pow((pp4->at(iP).vel[0]),2) + pow((pp4->at(iP).vel[1]),2) + pow((pp4->at(iP).vel[2]),2) );			

			/*for(j=0;j< hp->nbody1; j++){

				radius = sqrt( pow((pp4->at(iP).pos[0] - pp1->at(j).pos[0]),2) + pow((pp4->at(iP).pos[1] - pp1->at(j).pos[1]),2) + pow((pp4->at(iP).pos[2] - pp1->at(j).pos[2]),2) );  
				potDisk += -1.*(pp4->at(iP).mass * pp1->at(j).mass )/(radius);
			
			}*/	

			/*for(j=0;j< hp->nbody2; j++){

				radius = sqrt( pow((pp4->at(iP).pos[0] - pp2->at(j).pos[0]),2) + pow((pp4->at(iP).pos[1] - pp2->at(j).pos[1]),2) + pow((pp4->at(iP).pos[2] - pp2->at(j).pos[2]),2) );  
				potBulge += -1.*(pp4->at(iP).mass * pp2->at(j).mass )/(radius);
			
			}*/	

			/*for(j=0;j< hp->nbody3; j++){

				radius = sqrt( pow((pp4->at(iP).pos[0] - pp3->at(j).pos[0]),2) + pow((pp4->at(iP).pos[1] - pp3->at(j).pos[1]),2) + pow((pp4->at(iP).pos[2] - pp3->at(j).pos[2]),2) );  
				potHalo += -1.*(pp4->at(iP).mass * pp3->at(j).mass )/(radius);
			
			}*/	
		
			/*totDisk = potDisk + kin;
			totBulge = potBulge + kin;
			totHalo = potHalo + kin;

			if((totDisk >= 0) and (totHalo >= 0) and (totBulge >= 0)){
				pp4->at(iP).flag_zone = 0; //The particle is lost
			 			
		
			}else{	

				if((totDisk < totBulge) and (totDisk < totHalo)){
					pp4->at(iP).flag_zone = 1;  // it's a disk particle
					hp->sat_massDisk += pp4->at(iP).mass;
					hp->sat_massDisk_p4 += pp4->at(iP).mass;
				}else if((totBulge < totHalo) and (totBulge < totDisk)){			
					pp4->at(iP).flag_zone = 2;  //it's a bulge particle
					hp->sat_massBulge += pp4->at(iP).mass;
					hp->sat_massBulge_p4 += pp4->at(iP).mass;
				}else if((totHalo < totDisk) and (totHalo < totBulge)){	
					pp4->at(iP).flag_zone = 3;  //it's a halo particle		
					hp->sat_massHalo += pp4->at(iP).mass;	
					hp->sat_massHalo_p4 += pp4->at(iP).mass;
				}
	
			} //end if on particle-out-of-satellite condition	

		}else{

			pp4->at(iP).flag_zone = 4; // it's still a satellite particle. This is equivalent to flag bound == 0
		
		}*/

	if (hp->wantBaryons == 1){
		std::cout << "	 Doing zone flagging on particle 5" << std::endl;
		for(iP = 0; iP < hp->nbody5; iP++){
			if(pp5->at(iP).flag_bound == 0){
				rad3D = std::sqrt( std::pow(pp5->at(iP).pos[0],2) + std::pow(pp5->at(iP).pos[1],2) + std::pow(pp5->at(iP).pos[2],2));	
				rad = sqrt( std::pow(pp5->at(iP).pos[0],2) + std::pow(pp5->at(iP).pos[1],2));
				cosTheta = (pp5->at(iP).pos[0])/rad;	
				sinTheta = (pp5->at(iP).pos[1])/rad;		
				vRad = (pp5->at(iP).vel[0])*cosTheta + (pp5->at(iP).vel[1])*sinTheta;			
				vTan = -1.*(pp5->at(iP).vel[0])*sinTheta + (pp5->at(iP).vel[1])*cosTheta;
				zed = pp5->at(iP).pos[2];
				vZed = pp5->at(iP).vel[2];
				myTot_Bar += pp5->at(iP).mass;	
				//if ((rad >= 20.)){
				if ((rad >= 20.) or (std::abs(zed) >= 2.)){

					pp5->at(iP).flag_zone = 3; //halo particle
					hp->halo_Mass_p5_tot += pp5->at(iP).mass;
					myBar_inHalo += pp5->at(iP).mass;
				}else{

					if ((rad3D >= 1.0)){

						pp5->at(iP).flag_zone = 1; //disk particle
						myBar_inDisk += pp5->at(iP).mass;	

						/*

						for(iBin = 0;iBin<40;iBin++){
							if((rad >= (iBin*0.5)) and (rad < ((iBin+1.)*0.5))){

								//if( ( std::abs(zed-(hp->disk_Zed_p1[iBin])) < 3.*(hp->disk_SigZed_p1[iBin]) ) and ( std::abs(vZed-(hp->disk_Vzed_p1[iBin])) < 3.*(hp->disk_SigVzed_p1[iBin]) ) ){
								//if( abs(vZed-(hp->disk_Vzed_p1[iBin])) < 3.*(hp->disk_SigVzed_p1[iBin]) ){
							
									
									pp5->at(iP).flag_zone = 1; //disk particle
									
									

									hp->disk_Mass_p5[iBin] += pp5->at(iP).mass;
									disk_other_mass_p5[iBin] += pp5->at(iP).mass;
									hp->disk_Zed_p5[iBin] += pp5->at(iP).mass*(pp5->at(iP).pos[2]);
									hp->disk_Vzed_p5[iBin] += pp5->at(iP).mass*(pp5->at(iP).vel[2]);
									hp->disk_Vtan_p5[iBin]	 += pp5->at(iP).mass*vTan;
									if(iBin>0){

										disk_other_mass_p5[iBin-1] += pp5->at(iP).mass;
										hp->disk_Zed_p5[iBin-1] += pp5->at(iP).mass*(pp5->at(iP).pos[2]);
										hp->disk_Vzed_p5[iBin-1] += pp5->at(iP).mass*(pp5->at(iP).vel[2]);
										hp->disk_Vtan_p5[iBin-1] += pp5->at(iP).mass*vTan;

									}		
									if(iBin<39){

										disk_other_mass_p5[iBin+1] += pp5->at(iP).mass;
										hp->disk_Zed_p5[iBin+1] += pp5->at(iP).mass*(pp5->at(iP).pos[2]);
										hp->disk_Vzed_p5[iBin+1] += pp5->at(iP).mass*(pp5->at(iP).vel[2]);
										hp->disk_Vtan_p5[iBin+1] += pp5->at(iP).mass*vTan;

									}

										
		
									break;	
								
								}else{
			
									pp5->at(iP).flag_zone = 3; //halo particle
									hp->halo_Mass_p5_tot += pp5->at(iP).mass;
									break;

								}

							}

						}

						*/
					
					}else{

						pp5->at(iP).flag_zone = 2; //bulge particle
						myBar_inBulge += pp5->at(iP).mass;
						/*

						vel = std::sqrt( std::pow(pp5->at(iP).vel[0],2) + std::pow(pp5->at(iP).vel[1],2) + std::pow(pp5->at(iP).vel[2],2));
						for(iBin_b = 0;iBin_b<10;iBin_b++){

							if((rad3D >= (iBin_b*0.1)) and (rad3D < ((iBin_b+1.)*0.1))){
													
								
								if(vel<(3.*hp->bulge_Sigma_p2[iBin_b])){				
									pp5->at(iP).flag_zone = 2; //bulge particle				
									hp->bulge_Mass_p5[iBin_b] += pp5->at(iP).mass;

									//pp5->at(iP).flag_zone = 2; //bulge particle				
									//hp->bulge_Mass_p5[iBin_b] += pp5->at(iP).mass;
									break;

								

								}else{

									pp5->at(iP).flag_zone = 3; //halo particle
									hp->halo_Mass_p5_tot += pp5->at(iP).mass;
									break;
								}

								 	

							}
				
						}

						*/

					}
									
				}

			}  //End if p5 is bound or not

		} //End for loop on p5 particles


			/*if(iP%1000 == 0){

				std::cout << iP << std::endl;

			}
	
			if(pp5->at(iP).flag_bound == 0){ //so if the particle is out of satellite

				kin = 1./2.*(pp5->at(iP).mass)*( pow((pp5->at(iP).vel[0]),2) + pow((pp5->at(iP).vel[1]),2) + pow((pp5->at(iP).vel[2]),2) );			

				for(j=0;j< hp->nbody1; j++){

					radius = sqrt( pow((pp5->at(iP).pos[0] - pp1->at(j).pos[0]),2) + pow((pp5->at(iP).pos[1] - pp1->at(j).pos[1]),2) + pow((pp5->at(iP).pos[2] - pp1->at(j).pos[2]),2) );  
					potDisk += -1.*(pp5->at(iP).mass * pp1->at(j).mass )/(radius);
			
				}	

				for(j=0;j< hp->nbody2; j++){

					radius = sqrt( pow((pp5->at(iP).pos[0] - pp2->at(j).pos[0]),2) + pow((pp5->at(iP).pos[1] - pp2->at(j).pos[1]),2) + pow((pp5->at(iP).pos[2] - pp2->at(j).pos[2]),2) );  
					potBulge += -1.*(pp5->at(iP).mass * pp2->at(j).mass )/(radius);
			
				}	

				for(j=0;j< hp->nbody3; j++){

					radius = sqrt( pow((pp5->at(iP).pos[0] - pp3->at(j).pos[0]),2) + pow((pp5->at(iP).pos[1] - pp3->at(j).pos[1]),2) + pow((pp5->at(iP).pos[2] - pp3->at(j).pos[2]),2) );  
					potHalo += -1.*(pp5->at(iP).mass * pp3->at(j).mass )/(radius);
			
				}	
		
				totDisk = potDisk + kin;
				totBulge = potBulge + kin;
				totHalo = potHalo + kin;

				if((totDisk >= 0) and (totHalo >= 0) and (totBulge >= 0)){
					pp5->at(iP).flag_zone = 0; //The particle is lost
			 					
				}else{	

					if((totDisk < totBulge) and (totDisk < totHalo)){
						pp5->at(iP).flag_zone = 1;  // it's a disk particle	
						hp->sat_massDisk += pp5->at(iP).mass;
						hp->sat_massDisk_p5 += pp5->at(iP).mass;
					}else if((totBulge < totHalo) and (totBulge < totDisk)){	
						pp5->at(iP).flag_zone = 2;  //it's a bulge particle
						hp->sat_massBulge += pp5->at(iP).mass;
						hp->sat_massBulge_p5 += pp5->at(iP).mass;	
					}else if((totHalo < totDisk) and (totHalo < totBulge)){	
						pp5->at(iP).flag_zone = 3;  //it's a halo particle
						hp->sat_massHalo += pp5->at(iP).mass;
						hp->sat_massHalo_p5 += pp5->at(iP).mass;
					}
	
				} //end if on particle-out-of-satellite condition	

			}else{

				pp5->at(iP).flag_zone = 4; // it's still a satellite particle
		
			}*/
		
	} //end if on wantBaryons

	std::string fileStrN; 
	if (*thisFile < 10){
		fileStrN = "00"+myString(*thisFile);
	}else if((*thisFile >= 10) and (*thisFile < 100)){
		 fileStrN = "0"+myString(*thisFile);	
	}else{
		fileStrN = myString(*thisFile);
	}

	std::ofstream myFractionsTextFile((hp->myHDF5Root+"_snapshot_"+fileStrN+"_fractions_mass.txt").c_str());

	myFractionsTextFile <<  std::scientific << std::setprecision(6) << myDark_inDisk << "\t" << myDark_inBulge << "\t" << myDark_inHalo << "\t" << myTot_Dark << std::endl;
	myFractionsTextFile <<  std::scientific << std::setprecision(6) << myBar_inDisk << "\t" << myBar_inBulge << "\t" << myBar_inHalo << "\t" << myTot_Bar << std::endl;

	myFractionsTextFile.close();

}


//mass distribution on disk

void mass_distribution_disk(std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp){

	std::cout << " Calculating disk matter properties " << std::endl;

	int i1, i4, i5, iBin, zBin;
	float OldvRad, OldvTan, OldvZed, Oldzed;		
	float vRad, vTan, vZed, zed;
	float rad, theta;
	float cosTheta, sinTheta;

	float provvZed[40];

	float disk_other_mass_p1[40];
	float disk_other_mass_p4[40];
	float disk_other_mass_p5[40];

	for(int i=0;i<40;i++){

		disk_other_mass_p1[i] = 0.;
		disk_other_mass_p4[i] = 0.;
		disk_other_mass_p5[i] = 0.;
		provvZed[i] = 0.;

	}

	//Disk 

	for(i1=0;i1< hp->nbody1;i1++){
	
		//if(abs(pp1->at(i1).pos[2])<=2.0){ //1.96
			rad = sqrt( std::pow(pp1->at(i1).pos[0],2) + std::pow(pp1->at(i1).pos[1],2));
			cosTheta = (pp1->at(i1).pos[0])/rad;	
			sinTheta = (pp1->at(i1).pos[1])/rad;		
			vRad = (pp1->at(i1).vel[0])*cosTheta + (pp1->at(i1).vel[1])*sinTheta;
			vTan = -1.*(pp1->at(i1).vel[0])*sinTheta + (pp1->at(i1).vel[1])*cosTheta; 
			vZed = pp1->at(i1).vel[2];
			zed = pp1->at(i1).pos[2];
			for(iBin = 0;iBin<40;iBin++){

				if( (rad >= (iBin*0.5)) and (rad < ((iBin+1.)*0.5) ) ){
					hp->disk_Vtan_p1[iBin] += vTan*(pp1->at(i1).mass);
					hp->disk_Vzed_p1[iBin] += vZed*(pp1->at(i1).mass);
					provvZed[iBin] += zed*(pp1->at(i1).mass);
					hp->disk_Mass_p1[iBin] += pp1->at(i1).mass;
					disk_other_mass_p1[iBin] += pp1->at(i1).mass;
					
					if(iBin>0){
						hp->disk_Vtan_p1[iBin-1] += vTan*(pp1->at(i1).mass);
						hp->disk_Vzed_p1[iBin-1] += vZed*(pp1->at(i1).mass);
						provvZed[iBin-1] += zed*(pp1->at(i1).mass);
						disk_other_mass_p1[iBin-1] += pp1->at(i1).mass;
					}
					if(iBin<39){
						hp->disk_Vtan_p1[iBin+1] += vTan*(pp1->at(i1).mass);
						hp->disk_Vzed_p1[iBin+1] += vZed*(pp1->at(i1).mass);
						provvZed[iBin+1] += zed*(pp1->at(i1).mass);
						disk_other_mass_p1[iBin+1] += pp1->at(i1).mass;
					}		

					break;			
				}
			}
		//}
	}

	for(iBin = 0;iBin<40;iBin++){

		provvZed[iBin] /= disk_other_mass_p1[iBin];
		hp->disk_Mass_p1[iBin] = 0.;
		disk_other_mass_p1[iBin] = 0.;

	}

	for(i1=0;i1< hp->nbody1;i1++){

		//if(abs(pp1->at(i1).pos[2])<=2.0){ //1.96
			rad = sqrt( std::pow(pp1->at(i1).pos[0],2) + std::pow(pp1->at(i1).pos[1],2));
			cosTheta = (pp1->at(i1).pos[0])/rad;	
			sinTheta = (pp1->at(i1).pos[1])/rad;		
			vRad = (pp1->at(i1).vel[0])*cosTheta + (pp1->at(i1).vel[1])*sinTheta;		
			vTan = -1.*(pp1->at(i1).vel[0])*sinTheta + (pp1->at(i1).vel[1])*cosTheta; 
			vZed = pp1->at(i1).vel[2];
			zed = pp1->at(i1).pos[2];
			for(iBin = 0;iBin<40;iBin++){

				if( (rad >= (iBin*0.5)) and (rad < ((iBin+1.)*0.5) ) ){

				if (std::abs(zed - provvZed[iBin]) < 2.0){

					hp->disk_Vtan_p1[iBin] += vTan*(pp1->at(i1).mass);
					hp->disk_Vzed_p1[iBin] += vZed*(pp1->at(i1).mass);
					hp->disk_Zed_p1[iBin] += zed*(pp1->at(i1).mass);
					hp->disk_Mass_p1[iBin] += pp1->at(i1).mass;
					disk_other_mass_p1[iBin] += pp1->at(i1).mass;
					
					if(iBin>0){
						hp->disk_Vtan_p1[iBin-1] += vTan*(pp1->at(i1).mass);
						hp->disk_Vzed_p1[iBin-1] += vZed*(pp1->at(i1).mass);
						hp->disk_Zed_p1[iBin-1] += zed*(pp1->at(i1).mass);
						disk_other_mass_p1[iBin-1] += pp1->at(i1).mass;
					}
					if(iBin<39){
						hp->disk_Vtan_p1[iBin+1] += vTan*(pp1->at(i1).mass);
						hp->disk_Vzed_p1[iBin+1] += vZed*(pp1->at(i1).mass);
						hp->disk_Zed_p1[iBin+1] += zed*(pp1->at(i1).mass);
						disk_other_mass_p1[iBin+1] += pp1->at(i1).mass;
					}		

				}
				
				break;			

				}
			}
		//}
	}

	

	//Baryons 

	for(i4=0;i4< hp->nbody4;i4++){
	
	  //if((pp4->at(i4).flag_zone == 1) and (abs(pp4->at(i4).pos[2])<=1.96) and (pp4->at(i4).flag_bound == 0) ){
	  if((pp4->at(i4).flag_zone == 1) and (pp4->at(i4).flag_bound == 0) ){

			rad = sqrt( pow(pp4->at(i4).pos[0],2) + pow(pp4->at(i4).pos[1],2));
			if(pp4->at(i4).pos[0] != 0.){

				theta = std::atan((pp4->at(i4).pos[1])/(pp4->at(i4).pos[0]));			
			}else{
				theta = 0.;			
			}	
			vRad = (pp4->at(i4).vel[0])*cos(theta) + (pp4->at(i4).vel[1])*sin(theta);
			vZed = (pp4->at(i4).vel[2]);			
			vTan = -1.*(pp4->at(i4).vel[0])*sin(theta) + (pp4->at(i4).vel[1])*cos(theta); 
			zed = pp4->at(i4).pos[2];
			for(iBin = 0;iBin<40;iBin++){
	
				if((rad >= iBin*0.5) and (rad < (iBin+1.)*0.5)){

					hp->disk_Vtan_p4[iBin] += vTan*(pp4->at(i4).mass);
					hp->disk_Mass_p4[iBin] += pp4->at(i4).mass;
					hp->disk_Vzed_p4[iBin] += vZed*(pp4->at(i4).mass);
					hp->disk_Zed_p4[iBin] += zed*(pp4->at(i4).mass);	
					break;			
				}
			}
		}
	}

	//Dark 

	if(hp->wantBaryons == 1){

		for(i5=0;i5< hp->nbody5;i5++){

		  //if((pp5->at(i5).flag_zone == 1) and (abs(pp5->at(i5).pos[2])<=1.96) and (pp5->at(i5).flag_bound == 0) ){
		  if((pp5->at(i5).flag_zone == 1) and (pp5->at(i5).flag_bound == 0) ){
 				rad = sqrt( pow(pp5->at(i5).pos[0],2) + pow(pp5->at(i5).pos[1],2));
				if(pp5->at(i5).pos[0] != 0.){
					theta = std::atan((pp5->at(i5).pos[1])/(pp5->at(i5).pos[0]));	
				}else{
					theta = 0.;			
				}	

				vRad = (pp5->at(i5).vel[0])*cos(theta) + (pp5->at(i5).vel[1])*sin(theta);
				vTan = -1.*(pp5->at(i5).vel[0])*sin(theta) + (pp5->at(i5).vel[1])*cos(theta);
				vZed = (pp5->at(i5).vel[2]); 
				zed = pp5->at(i5).pos[2];	
				for(iBin = 0;iBin<40;iBin++){
					if((rad >= iBin*0.5) and (rad < (iBin+1.)*0.5)){
						hp->disk_Vtan_p5[iBin] += vTan*(pp5->at(i5).mass);
						hp->disk_Vzed_p5[iBin] += vZed*(pp5->at(i5).mass);	
						hp->disk_Mass_p5[iBin] += pp5->at(i5).mass;
						hp->disk_Zed_p5[iBin] += zed*(pp5->at(i5).mass);
						break;			
					}
				}	
			}
		}
	}

	//Final divisions		
	
	for(iBin = 0; iBin<40;iBin++){

		//Final calculation disk vTan	

		if(disk_other_mass_p1[iBin] != 0){

			hp->disk_Vtan_p1[iBin] /= disk_other_mass_p1[iBin];
			hp->disk_Vzed_p1[iBin] /= disk_other_mass_p1[iBin];
			hp->disk_Zed_p1[iBin]  /= disk_other_mass_p1[iBin];
			
		}else{				
			hp->disk_Vtan_p1[iBin] = 0.;
			hp->disk_Vzed_p1[iBin] = 0.;
			hp->disk_Zed_p1[iBin]  = 0.;
		}

	}
	

	for(i1=0;i1< hp->nbody1;i1++){
		
		//if(abs(pp1->at(i1).pos[2])<=2.0){ //1.96
			rad = std::sqrt( std::pow(pp1->at(i1).pos[0],2) + std::pow(pp1->at(i1).pos[1],2));
			cosTheta = (pp1->at(i1).pos[0])/rad;	
			sinTheta = (pp1->at(i1).pos[1])/rad;		
			vRad = (pp1->at(i1).vel[0])*cosTheta + (pp1->at(i1).vel[1])*sinTheta;	
			vTan = -1.*(pp1->at(i1).vel[0])*sinTheta + (pp1->at(i1).vel[1])*cosTheta;  
			vZed = pp1->at(i1).vel[2];
			zed = pp1->at(i1).pos[2];
			for(iBin = 0;iBin<40;iBin++){
		
				if( (rad >= (iBin*0.5)) and (rad < ((iBin+1.)*0.5) ) ){
					
					float vDiffer = std::pow((vZed-hp->disk_Vzed_p1[iBin]),2);
					float zDiffer = std::pow((zed-hp->disk_Zed_p1[iBin]),2);
					hp->disk_SigVzed_p1[iBin] += vDiffer*(pp1->at(i1).mass);
					hp->disk_SigZed_p1[iBin] += zDiffer*(pp1->at(i1).mass);
					
					if(zDiffer < 4.){

					if(iBin>0){
						
						hp->disk_SigVzed_p1[iBin-1] += vDiffer*(pp1->at(i1).mass);
						hp->disk_SigZed_p1[iBin-1] += zDiffer*(pp1->at(i1).mass);
					}
					if(iBin<39){
						hp->disk_SigVzed_p1[iBin+1] += vDiffer*(pp1->at(i1).mass);
						hp->disk_SigZed_p1[iBin+1] += zDiffer*(pp1->at(i1).mass);
					}		
		
					}	

					break;			
				}
			}
		//}
	}

	for(iBin = 0; iBin<40;iBin++){	

		if(disk_other_mass_p1[iBin] != 0){

			hp->disk_SigVzed_p1[iBin] = std::sqrt((hp->disk_SigVzed_p1[iBin])/(disk_other_mass_p1[iBin]));
			hp->disk_SigZed_p1[iBin] = std::sqrt((hp->disk_SigZed_p1[iBin])/(disk_other_mass_p1[iBin]));
									
		}else{				
			hp->disk_SigVzed_p1[iBin] = 0.;
			hp->disk_SigZed_p1[iBin] = 0.;
		}

	}


	for(iBin = 0; iBin<40;iBin++){

       	//Final calculation baryon vTan	
	
		if(hp->disk_Mass_p4[iBin] != 0){
                    hp->disk_Vtan_p4[iBin] /= hp->disk_Mass_p4[iBin];	    
   		    hp->disk_Zed_p4[iBin] /= hp->disk_Mass_p4[iBin];
                    hp->disk_Vzed_p4[iBin] /= hp->disk_Mass_p4[iBin];
		}else{				
			hp->disk_Vtan_p4[iBin] = 0.;
			hp->disk_Vzed_p4[iBin] = 0.;
			hp->disk_Zed_p4[iBin] = 0.;
		}

		if(hp->wantBaryons == 1){	

			//Final calculation dark vTan	
	
			if(hp->disk_Mass_p5[iBin] != 0){

			  hp->disk_Vtan_p5[iBin] /= hp->disk_Mass_p5[iBin];	
                          hp->disk_Zed_p5[iBin] /= hp->disk_Mass_p5[iBin];
			  hp->disk_Vzed_p5[iBin] /= hp->disk_Mass_p5[iBin];		

			}else{
				hp->disk_Vtan_p5[iBin] = 0.;
                                hp->disk_Vzed_p5[iBin] = 0.;
                                hp->disk_Zed_p5[iBin] = 0.;
			}
		}
	}

	
	//Now we calculate how much mass lies in vertical for all disk extension, directly, then you separate in radii
	//Total

	/*
	//Disk
	
	for(i1=0;i1< hp->nbody1;i1++){		

		for(zBin=0;zBin<7;zBin++){

			if( ( (pp1->at(i1).pos[2]) >= zBin*0.280) and ( (pp1->at(i1).pos[2]) < (zBin+1)*0.280) ){
				hp->disk_vertical_total_p1[zBin] += pp1->at(i1).mass;
				break;	
			}else if( (pp1->at(i1).pos[2]) >= 2.0){
			
				hp->disk_vertical_total_p1[7] += pp1->at(i1).mass;
			}		

	}

	//Baryonic

	for(i4=0;i4< hp->nbody4;i4++){		

		if((pp4->at(i4).flag_zone == 1) and (pp4->at(i4).flag_bound == 0)) {

			for(zBin=0;zBin<7;zBin++){

				if( ( (pp4->at(i4).pos[2]) >= zBin*0.280) and ( (pp4->at(i4).pos[2]) < (zBin+1)*0.280) ){
		
					hp->disk_vertical_total_p4[zBin] += pp4->at(i4).mass;
					break;	
				}else if( (pp4->at(i4).pos[2]) >= 2.0){
			
					hp->disk_vertical_total_p4[7] += pp4->at(i4).mass;
				}
			}

		}			

	}

	//Dark

	if(hp->wantBaryons == 1){

		for(i5=0;i5< hp->nbody5;i5++){		

			if((pp5->at(i5).flag_zone == 1) and (pp5->at(i5).flag_bound == 0)){

				for(zBin=0;zBin<7;zBin++){

					if( ( (pp5->at(i5).pos[2]) >= zBin*0.280) and ( (pp5->at(i5).pos[2]) < (zBin+1)*0.280) ){
		
						hp->disk_vertical_total_p5[zBin] += pp5->at(i5).mass;
						break;	
					}else if( (pp5->at(i5).pos[2]) >= 2.0){
				
						hp->disk_vertical_total_p5[7] += pp5->at(i5).mass;
					}
				}
			}			
		}	
	}


	/*       Now we separate in radii	*/
 
	//Disk

	/*
	
	for(i1=0;i1< hp->nbody1;i1++){		

		for(iBin=0;iBin<20;iBin++){

			rad = sqrt( pow(pp1->at(i1).pos[0],2) + pow(pp1->at(i1).pos[1],2));	
			if((rad >= iBin*0.5) and (rad < (iBin+1.)*0.5)){
			
				if( abs(pp1->at(i1).pos[2]) < 0.280 ){
					hp->disk_vertical_zone_p1[0][iBin] += pp1->at(i1).mass;
					break;	
				}else if( ( abs(pp1->at(i1).pos[2]) >= 0.280 ) and ( abs(pp1->at(i1).pos[2]) < 1.0 ) ){
					hp->disk_vertical_zone_p1[1][iBin] += pp1->at(i1).mass;
					break;	
				}else if( ( abs(pp1->at(i1).pos[2]) >= 1.000 ) and ( abs(pp1->at(i1).pos[2]) < 1.96 ) ){
					hp->disk_vertical_zone_p1[2][iBin] += pp1->at(i1).mass;					
					break;
				}else{
					hp->disk_vertical_zone_p1[3][iBin] += pp1->at(i1).mass;					
					break;					
				}

			}

		}

	}

		

	//Baryonic

	for(i4=0;i4< hp->nbody4;i4++){		

		if((pp4->at(i4).flag_zone == 1) and (pp4->at(i4).flag_bound == 0)){

			rad = sqrt( pow(pp4->at(i4).pos[0],2) + pow(pp4->at(i4).pos[1],2));

			for(iBin=0;iBin<20;iBin++){

				if((rad >= iBin*0.5) and (rad < (iBin+1.)*0.5)){
			
					if( abs(pp4->at(i4).pos[2]) < 0.280 ){
						hp->disk_vertical_zone_p4[0][iBin] += pp4->at(i4).mass;
						break;	
					}else if( ( abs(pp4->at(i4).pos[2]) >= 0.280 ) and ( abs(pp4->at(i4).pos[2]) < 1.0 ) ){
						hp->disk_vertical_zone_p4[1][iBin] += pp4->at(i4).mass;
						break;	
					}else if( ( abs(pp4->at(i4).pos[2]) >= 1.000 ) and ( abs(pp4->at(i4).pos[2]) < 1.96 ) ){
						hp->disk_vertical_zone_p4[2][iBin] += pp4->at(i4).mass;					
						break;
					}else{
						hp->disk_vertical_zone_p4[3][iBin] += pp4->at(i4).mass;					
						break;					
					}
	
				}
	
			}
	
		}

	}

	//Dark

	if(hp->wantBaryons == 1){

		for(i5=0;i5< hp->nbody5;i5++){	

			if((pp5->at(i5).flag_zone == 1)and (pp5->at(i5).flag_bound == 0)){	

				rad = sqrt( pow(pp5->at(i5).pos[0],2) + pow(pp5->at(i5).pos[1],2));
	
				for(iBin=0;iBin<20;iBin++){

					if((rad >= iBin*0.5) and (rad < (iBin+1.)*0.5)){
			
						if( abs(pp5->at(i5).pos[2]) < 0.280 ){
							hp->disk_vertical_zone_p5[0][iBin] += pp5->at(i5).mass;
							break;	
						}else if( ( abs(pp5->at(i5).pos[2]) >= 0.280 ) and ( abs(pp5->at(i5).pos[2]) < 1.0 ) ){
							hp->disk_vertical_zone_p5[1][iBin] += pp5->at(i5).mass;
							break;	
						}else if( ( abs(pp5->at(i5).pos[2]) >= 1.000 ) and ( abs(pp5->at(i5).pos[2]) < 1.96 ) ){
							hp->disk_vertical_zone_p5[2][iBin] += pp5->at(i5).mass;					
							break;
						}else{
							hp->disk_vertical_zone_p5[3][iBin] += pp5->at(i5).mass;					
							break;					
						}
					}
				}
			}		
		}
	}
	
	*/
}		



void mass_distribution_bulge(std::vector<struct h5_particle> *pp2, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp){

	std::cout << " Calculating bulge matter properties " << std::endl;

	int i2, i4, i5, iBin, vBin;	
	float rad, radxy;
	float theta;
	float phiAng;
	float vSQ;
	float vPhi, vTheta, vRad;
	float sinTheta, sinPhiAng, cosTheta, cosPhiAng;

	float bulge_other_mass_p2[10];
	float bulge_other_mass_p4[10];
	float bulge_other_mass_p5[10];

	for(int i=0;i<10;i++){

		bulge_other_mass_p2[i] = 0.;
		bulge_other_mass_p4[i] = 0.;
		bulge_other_mass_p5[i] = 0.;

	}

	//Bulge 

	for(i2=0;i2< hp->nbody2;i2++){
	
		rad = std::sqrt( pow(pp2->at(i2).pos[0],2) + std::pow(pp2->at(i2).pos[1],2) + std::pow(pp2->at(i2).pos[2],2));
		radxy = std::sqrt( pow(pp2->at(i2).pos[0],2) + std::pow(pp2->at(i2).pos[1],2));
		cosTheta = (pp2->at(i2).pos[2])/rad;	
		sinTheta = radxy/rad;
		cosPhiAng = (pp2->at(i2).pos[0])/radxy;			
		sinPhiAng = (pp2->at(i2).pos[1])/radxy;			
		//vRad = (pp2->at(i2).vel[0])*sin(theta)*cos(phiAng) + (pp2->at(i2).vel[1])*sin(theta)*sin(phiAng) + (pp2->at(i2).vel[2])*cos(theta);
		vRad = (pp2->at(i2).vel[0])*sinTheta*cosPhiAng + (pp2->at(i2).vel[1])*sinTheta*sinPhiAng + (pp2->at(i2).vel[2])*cosTheta;		
		//vTheta = (pp2->at(i2).vel[0])*cos(theta)*sin(phiAng) + (pp2->at(i2).vel[1])*cos(theta)*sin(phiAng) - sin(theta)*(pp2->at(i2).vel[2]); 
		vTheta = (pp2->at(i2).vel[0])*cosTheta*sinPhiAng + (pp2->at(i2).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp2->at(i2).vel[2]); 		
		//vPhi = -1.*(pp2->at(i2).vel[0])*sin(phiAng) + (pp2->at(i2).vel[1])*cos(phiAng);	
		vPhi = -1.*(pp2->at(i2).vel[0])*sinPhiAng + (pp2->at(i2).vel[1])*cosPhiAng;	
		vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);
		for(iBin = 0;iBin<10;iBin++){
			if((rad >= iBin*0.1) and (rad < (iBin+1.)*0.1)){
				hp->bulge_Vrad_p2[iBin] += vRad*(pp2->at(i2).mass);		
				hp->bulge_Vtheta_p2[iBin] += vTheta*(pp2->at(i2).mass);
				hp->bulge_Vphi_p2[iBin] += vPhi*(pp2->at(i2).mass);	
				hp->bulge_Mass_p2[iBin] += pp2->at(i2).mass;
				bulge_other_mass_p2[iBin] += pp2->at(i2).mass;
				if(iBin>0){

					hp->bulge_Vrad_p2[iBin-1] += vRad*(pp2->at(i2).mass);		
					hp->bulge_Vtheta_p2[iBin-1] += vTheta*(pp2->at(i2).mass);
					hp->bulge_Vphi_p2[iBin-1] += vPhi*(pp2->at(i2).mass);	
					bulge_other_mass_p2[iBin-1] += pp2->at(i2).mass;
				}

				if(iBin<9){

					hp->bulge_Vrad_p2[iBin+1] += vRad*(pp2->at(i2).mass);		
					hp->bulge_Vtheta_p2[iBin+1] += vTheta*(pp2->at(i2).mass);
					hp->bulge_Vphi_p2[iBin+1] += vPhi*(pp2->at(i2).mass);	
					bulge_other_mass_p2[iBin+1] += pp2->at(i2).mass;
				}


				//for(vBin=0;vBin<40;vBin++){
				//	if( (vSQ >= (pow(vBin*5.08,2))) and (vSQ <(pow((vBin+1)*5.08,2)) ) ){
				//		hp->halo_vels_p2[vBin] += pp2->at(i2).mass; 
				//		break;						
				//	}
				//}	
				break;			
			}
		
		}

	}

	//Baryons 


	//Final divisions		
	
	for(iBin = 0; iBin<10;iBin++){

		//Final calculation disk vRad	
		if(bulge_other_mass_p2[iBin] != 0){
			hp->bulge_Vrad_p2[iBin] /= bulge_other_mass_p2[iBin];
		}else{				
			hp->bulge_Vrad_p2[iBin] = 0.;
		}

		//Final calculation disk vTheta	
		if(bulge_other_mass_p2[iBin] != 0){
			hp->bulge_Vtheta_p2[iBin] /= bulge_other_mass_p2[iBin];	
		}else{				
			hp->bulge_Vtheta_p2[iBin] = 0.;
		}

		//Final calculation disk vPhi	
		if(bulge_other_mass_p2[iBin] != 0){
			hp->bulge_Vphi_p2[iBin] /= bulge_other_mass_p2[iBin];	
		}else{				
			hp->bulge_Vphi_p2[iBin] = 0.;
		}

	}

	//Bulge 

	for(i2=0;i2< hp->nbody2;i2++){
	
		rad = std::sqrt( pow(pp2->at(i2).pos[0],2) + std::pow(pp2->at(i2).pos[1],2) + pow(pp2->at(i2).pos[2],2));
		radxy = std::sqrt( pow(pp2->at(i2).pos[0],2) + std::pow(pp2->at(i2).pos[1],2));
		cosTheta = (pp2->at(i2).pos[2])/rad;	
		sinTheta = radxy/rad;
		cosPhiAng = (pp2->at(i2).pos[0])/radxy;			
		sinPhiAng = (pp2->at(i2).pos[1])/radxy;			
		//vRad = (pp2->at(i2).vel[0])*sin(theta)*cos(phiAng) + (pp2->at(i2).vel[1])*sin(theta)*sin(phiAng) + (pp2->at(i2).vel[2])*cos(theta);
		vRad = (pp2->at(i2).vel[0])*sinTheta*cosPhiAng + (pp2->at(i2).vel[1])*sinTheta*sinPhiAng + (pp2->at(i2).vel[2])*cosTheta;		
		//vTheta = (pp2->at(i2).vel[0])*cos(theta)*sin(phiAng) + (pp2->at(i2).vel[1])*cos(theta)*sin(phiAng) - sin(theta)*(pp2->at(i2).vel[2]); 
		vTheta = (pp2->at(i2).vel[0])*cosTheta*sinPhiAng + (pp2->at(i2).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp2->at(i2).vel[2]); 		
		//vPhi = -1.*(pp2->at(i2).vel[0])*sin(phiAng) + (pp2->at(i2).vel[1])*cos(phiAng);	
		vPhi = -1.*(pp2->at(i2).vel[0])*sinPhiAng + (pp2->at(i2).vel[1])*cosPhiAng;	
		for(iBin = 0;iBin<10;iBin++){
			if((rad >= iBin*0.1) and (rad < (iBin+1.)*0.1)){
				hp->bulge_Sigrad_p2[iBin] += std::pow((vRad - hp->bulge_Vrad_p2[iBin]),2)*(pp2->at(i2).mass);		
				hp->bulge_Sigtheta_p2[iBin] += std::pow((vTheta - hp->bulge_Vtheta_p2[iBin]),2)*(pp2->at(i2).mass);
				hp->bulge_Sigphi_p2[iBin] += std::pow((vPhi - hp->bulge_Vphi_p2[iBin]),2)*(pp2->at(i2).mass);
				
				if(iBin>0){
				hp->bulge_Sigrad_p2[iBin-1] += std::pow((vRad - hp->bulge_Vrad_p2[iBin]),2)*(pp2->at(i2).mass);		
				hp->bulge_Sigtheta_p2[iBin-1] += std::pow((vTheta - hp->bulge_Vtheta_p2[iBin]),2)*(pp2->at(i2).mass);
				hp->bulge_Sigphi_p2[iBin-1] += std::pow((vPhi - hp->bulge_Vphi_p2[iBin]),2)*(pp2->at(i2).mass);		
				}
				if(iBin<9){
				hp->bulge_Sigrad_p2[iBin+1] += std::pow((vRad - hp->bulge_Vrad_p2[iBin]),2)*(pp2->at(i2).mass);		
				hp->bulge_Sigtheta_p2[iBin+1] += std::pow((vTheta - hp->bulge_Vtheta_p2[iBin]),2)*(pp2->at(i2).mass);
				hp->bulge_Sigphi_p2[iBin+1] += std::pow((vPhi - hp->bulge_Vphi_p2[iBin]),2)*(pp2->at(i2).mass);		
				}
	

				break;			
			}
		}

	}

	//Final divisions		
	
	for(iBin = 0; iBin<10;iBin++){

		//Final calculation bulge vRad	
		if(bulge_other_mass_p2[iBin] != 0){
			hp->bulge_Sigrad_p2[iBin] = sqrt(hp->bulge_Sigrad_p2[iBin]/bulge_other_mass_p2[iBin]);	
		}else{				
			hp->bulge_Sigrad_p2[iBin] = 0.;
		}

		//Final calculation bulge vTheta	
		if(bulge_other_mass_p2[iBin] != 0){
			hp->bulge_Sigtheta_p2[iBin] = sqrt(hp->bulge_Sigtheta_p2[iBin]/bulge_other_mass_p2[iBin]);
		}else{				
			hp->bulge_Sigtheta_p2[iBin] = 0.;
		}

		//Final calculation bulge vPhi	
		if(bulge_other_mass_p2[iBin] != 0){
			hp->bulge_Sigphi_p2[iBin] = sqrt(hp->bulge_Sigphi_p2[iBin]/bulge_other_mass_p2[iBin]);
		}else{				
			hp->bulge_Sigphi_p2[iBin] = 0.;
		}

		hp->bulge_Sigma_p2[iBin] = std::sqrt( (hp->bulge_Sigrad_p2[iBin])*(hp->bulge_Sigrad_p2[iBin]) + (hp->bulge_Sigtheta_p2[iBin])*(hp->bulge_Sigtheta_p2[iBin]) + (hp->bulge_Sigphi_p2[iBin])*(hp->bulge_Sigphi_p2[iBin]) );

	}

	//particle 4

	for(i4=0;i4< hp->nbody4;i4++){
	
		if((pp4->at(i4).flag_zone == 2) and (pp4->at(i4).flag_bound == 0) ){

		rad = std::sqrt( pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2) + std::pow(pp4->at(i4).pos[2],2));
		radxy = std::sqrt( pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2));
		cosTheta = (pp4->at(i4).pos[2])/rad;	
		sinTheta = radxy/rad;
		cosPhiAng = (pp4->at(i4).pos[0])/radxy;			
		sinPhiAng = (pp4->at(i4).pos[1])/radxy;			
		//vRad = (pp2->at(i2).vel[0])*sin(theta)*cos(phiAng) + (pp2->at(i2).vel[1])*sin(theta)*sin(phiAng) + (pp2->at(i2).vel[2])*cos(theta);
		vRad = (pp4->at(i4).vel[0])*sinTheta*cosPhiAng + (pp4->at(i4).vel[1])*sinTheta*sinPhiAng + (pp4->at(i4).vel[2])*cosTheta;		
		//vTheta = (pp2->at(i2).vel[0])*cos(theta)*sin(phiAng) + (pp2->at(i2).vel[1])*cos(theta)*sin(phiAng) - sin(theta)*(pp2->at(i2).vel[2]); 
		vTheta = (pp4->at(i4).vel[0])*cosTheta*sinPhiAng + (pp4->at(i4).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp4->at(i4).vel[2]); 		
		//vPhi = -1.*(pp2->at(i2).vel[0])*sin(phiAng) + (pp2->at(i2).vel[1])*cos(phiAng);	
		vPhi = -1.*(pp4->at(i4).vel[0])*sinPhiAng + (pp4->at(i4).vel[1])*cosPhiAng;	
		vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);
		for(iBin = 0;iBin<10;iBin++){
			if((rad >= iBin*0.1) and (rad < (iBin+1.)*0.1)){
				hp->bulge_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
				hp->bulge_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
				hp->bulge_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
				hp->bulge_Mass_p4[iBin] += pp4->at(i4).mass;
				bulge_other_mass_p4[iBin] += pp4->at(i4).mass;
				if(iBin>0){

					hp->bulge_Vrad_p4[iBin-1] += vRad*(pp4->at(i4).mass);		
					hp->bulge_Vtheta_p4[iBin-1] += vTheta*(pp4->at(i4).mass);
					hp->bulge_Vphi_p4[iBin-1] += vPhi*(pp4->at(i4).mass);	
					bulge_other_mass_p4[iBin-1] += pp4->at(i4).mass;
				}

				if(iBin<9){

					hp->bulge_Vrad_p4[iBin+1] += vRad*(pp4->at(i4).mass);		
					hp->bulge_Vtheta_p4[iBin+1] += vTheta*(pp4->at(i4).mass);
					hp->bulge_Vphi_p4[iBin+1] += vPhi*(pp4->at(i4).mass);	
					bulge_other_mass_p4[iBin+1] += pp4->at(i4).mass;
				}


				//for(vBin=0;vBin<40;vBin++){
				//	if( (vSQ >= (pow(vBin*5.08,2))) and (vSQ <(pow((vBin+1)*5.08,2)) ) ){
				//		hp->halo_vels_p2[vBin] += pp2->at(i2).mass; 
				//		break;						
				//	}
				//}	
				break;			
			}
		
		}

		}

	}


	//particle 5

	for(i5=0;i5< hp->nbody5;i5++){
	
		if((pp5->at(i5).flag_zone == 2) and (pp5->at(i5).flag_bound == 0) ){

		rad = std::sqrt( pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2) + std::pow(pp5->at(i5).pos[2],2));
		radxy = std::sqrt( pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2));
		cosTheta = (pp5->at(i5).pos[2])/rad;	
		sinTheta = radxy/rad;
		cosPhiAng = (pp5->at(i5).pos[0])/radxy;			
		sinPhiAng = (pp5->at(i5).pos[1])/radxy;			
		//vRad = (pp2->at(i2).vel[0])*sin(theta)*cos(phiAng) + (pp2->at(i2).vel[1])*sin(theta)*sin(phiAng) + (pp2->at(i2).vel[2])*cos(theta);
		vRad = (pp5->at(i5).vel[0])*sinTheta*cosPhiAng + (pp5->at(i5).vel[1])*sinTheta*sinPhiAng + (pp5->at(i5).vel[2])*cosTheta;		
		//vTheta = (pp2->at(i2).vel[0])*cos(theta)*sin(phiAng) + (pp2->at(i2).vel[1])*cos(theta)*sin(phiAng) - sin(theta)*(pp2->at(i2).vel[2]); 
		vTheta = (pp5->at(i5).vel[0])*cosTheta*sinPhiAng + (pp5->at(i5).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp5->at(i5).vel[2]); 		
		//vPhi = -1.*(pp2->at(i2).vel[0])*sin(phiAng) + (pp2->at(i2).vel[1])*cos(phiAng);	
		vPhi = -1.*(pp5->at(i5).vel[0])*sinPhiAng + (pp5->at(i5).vel[1])*cosPhiAng;	
		vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);
		for(iBin = 0;iBin<10;iBin++){
			if((rad >= iBin*0.1) and (rad < (iBin+1.)*0.1)){
				hp->bulge_Vrad_p5[iBin] += vRad*(pp5->at(i5).mass);		
				hp->bulge_Vtheta_p5[iBin] += vTheta*(pp5->at(i5).mass);
				hp->bulge_Vphi_p5[iBin] += vPhi*(pp5->at(i5).mass);	
				hp->bulge_Mass_p5[iBin] += pp5->at(i5).mass;
				bulge_other_mass_p5[iBin] += pp5->at(i5).mass;
				if(iBin>0){

					hp->bulge_Vrad_p5[iBin-1] += vRad*(pp5->at(i5).mass);		
					hp->bulge_Vtheta_p5[iBin-1] += vTheta*(pp5->at(i5).mass);
					hp->bulge_Vphi_p5[iBin-1] += vPhi*(pp5->at(i5).mass);	
					bulge_other_mass_p5[iBin-1] += pp5->at(i5).mass;
				}

				if(iBin<9){

					hp->bulge_Vrad_p5[iBin+1] += vRad*(pp5->at(i5).mass);		
					hp->bulge_Vtheta_p5[iBin+1] += vTheta*(pp5->at(i5).mass);
					hp->bulge_Vphi_p5[iBin+1] += vPhi*(pp5->at(i5).mass);	
					bulge_other_mass_p5[iBin+1] += pp5->at(i5).mass;
				}


				//for(vBin=0;vBin<40;vBin++){
				//	if( (vSQ >= (pow(vBin*5.08,2))) and (vSQ <(pow((vBin+1)*5.08,2)) ) ){
				//		hp->halo_vels_p2[vBin] += pp2->at(i2).mass; 
				//		break;						
				//	}
				//}	
				break;			
			}
		
		}

		}

	}

	for(iBin = 0;iBin<10;iBin++){

		//Final calculation baryon vRad	
		if(hp->disk_Mass_p4[iBin] != 0){
			hp->bulge_Vrad_p4[iBin] = sqrt(hp->bulge_Vrad_p4[iBin]/hp->bulge_Mass_p4[iBin]);		
		}else{				
			hp->bulge_Vrad_p4[iBin] = 0.;
		}

		if(hp->wantBaryons == 1){

			//Final calculation dark vRad	
			if(hp->disk_Mass_p5[iBin] != 0){
				hp->bulge_Vrad_p5[iBin] = sqrt(hp->bulge_Vrad_p5[iBin]/hp->bulge_Mass_p5[iBin]);
			}else{
				hp->bulge_Vrad_p5[iBin] = 0.;
			}
		}

		
		//Final calculation baryon vTheta	
		if(hp->bulge_Mass_p4[iBin] != 0){
			hp->bulge_Vtheta_p4[iBin] = sqrt(hp->bulge_Vtheta_p4[iBin]/hp->bulge_Mass_p4[iBin]);	
		}else{				
			hp->bulge_Vtheta_p4[iBin] = 0.;
		}

		if(hp->wantBaryons == 1){

			//Final calculation dark vTheta	
			if(hp->bulge_Mass_p5[iBin] != 0){
				hp->bulge_Vtheta_p5[iBin] = sqrt(hp->bulge_Vtheta_p5[iBin]/hp->bulge_Mass_p5[iBin]);
			}else{
				hp->bulge_Vtheta_p5[iBin] = 0.;
			}
		}

		
		//Final calculation baryon vPhi	
		if(hp->bulge_Mass_p4[iBin] != 0){
			hp->bulge_Vphi_p4[iBin] = sqrt(hp->bulge_Vphi_p4[iBin]/hp->bulge_Mass_p4[iBin]);		
		}else{				
			hp->bulge_Vphi_p4[iBin] = 0.;
		}
	
		if(hp->wantBaryons == 1){

			//Final calculation dark vPhi	
			if(hp->bulge_Mass_p5[iBin] != 0){
				hp->bulge_Vphi_p5[iBin] = sqrt(hp->bulge_Vphi_p5[iBin]/hp->bulge_Mass_p5[iBin]);	
			}else{
				hp->bulge_Vphi_p5[iBin] = 0.;
			}
		}
	}

}

//Mass distribution in space and velocity in halo

	
void mass_distribution_halo(std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, struct header_h5 *hp){

	float G_const = 1.; //4.30e-6; //const of gravitation 
	float R200 =  242.8; //assume c.a. the r_200 for the halo of the milky way; r_scale= 25 kpc
	float c_par = hp->MW_concentration;
	float r_scale = R200/c_par;
	float M_tot = hp->tot_halo_mass;
	float M_hern;
	float mu, mfac, myPhi, En, Lang_x, Lang_y, Lang_z, LL;
	float a_halo = r_scale*std::sqrt( 2.*( log(1.+c_par)-c_par/(1.+c_par) ) );
	int i3, i4, i5, iBin,vBin, iBin_hiRes, iBin_hiRes2, iBin_hiRes0;	
	float rad, radxy;
	float theta, cosTheta, sinTheta, cosPhiAng, sinPhiAng;
	float phiAng;
	float vRad, vTheta, vPhi, vSQ;

	hp->mass_out_of_242_p4 = 0.;
	hp->mass_out_of_242_p5 = 0.;

	//Halo particles 

	for(i3=0;i3< hp->nbody3;i3++){

		pp3->at(i3).flag_bound = 1; //we flag all halo particles as bound for simplicity
	
		rad = std::sqrt( std::pow(pp3->at(i3).pos[0],2) + std::pow(pp3->at(i3).pos[1],2) + std::pow(pp3->at(i3).pos[2],2));
		radxy = std::sqrt( std::pow(pp3->at(i3).pos[0],2) + std::pow(pp3->at(i3).pos[1],2));
		cosTheta = (pp3->at(i3).pos[2])/rad;	
		sinTheta = radxy/rad;
		cosPhiAng = (pp3->at(i3).pos[0])/radxy;			
		sinPhiAng = (pp3->at(i3).pos[1])/radxy;			
		vRad = (pp3->at(i3).vel[0])*sinTheta*cosPhiAng + (pp3->at(i3).vel[1])*sinTheta*sinPhiAng + (pp3->at(i3).vel[2])*cosTheta;		
		vTheta = (pp3->at(i3).vel[0])*cosTheta*sinPhiAng + (pp3->at(i3).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp3->at(i3).vel[2]); 		
		vPhi = -1.*(pp3->at(i3).vel[0])*sinPhiAng + (pp3->at(i3).vel[1])*cosPhiAng;	
		vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);
		for(iBin = 0;iBin<40;iBin++){
			if((rad >= iBin*6.07) and (rad < (iBin+1.)*6.07)){
					M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
					mfac =  M_hern*(pp3->at(i3).mass);
					mu =  mfac/(M_hern+(pp3->at(i3).mass));		
					En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
					Lang_x = ( (pp3->at(i3).pos[1])*(pp3->at(i3).vel[2]) - (pp3->at(i3).pos[2])*(pp3->at(i3).vel[1]) )*mu;
					Lang_y = ( (pp3->at(i3).pos[2])*(pp3->at(i3).vel[0]) - (pp3->at(i3).pos[0])*(pp3->at(i3).vel[2]) )*mu;
					Lang_z = ( (pp3->at(i3).pos[0])*(pp3->at(i3).vel[1]) - (pp3->at(i3).pos[1])*(pp3->at(i3).vel[0]) )*mu;	
					LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
					hp->halo_Vrad_p3[iBin] += vRad*(pp3->at(i3).mass);		
					hp->halo_Vtheta_p3[iBin] += vTheta*(pp3->at(i3).mass);
					hp->halo_Vphi_p3[iBin] += vPhi*(pp3->at(i3).mass);	
					hp->halo_ecc_p3[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp3->at(i3).mass);
					hp->smj_ax_p3[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p3[iBin])*(hp->halo_ecc_p3[iBin]) )*G_const*M_hern )*(pp3->at(i3).mass); 
					hp->halo_Mass_p3[iBin] += pp3->at(i3).mass;					
					break;						
			}
		}

	        for(iBin_hiRes = 0;iBin_hiRes<50;iBin_hiRes++){

			if((rad >= iBin_hiRes*0.5) and (rad < (iBin_hiRes+1.)*0.5)){
	
				/*

				M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
				mfac =  M_hern*(pp4->at(i4).mass);
				mu =  mfac/(M_hern+(pp4->at(i4).mass));		
				En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
				Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
				Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
				Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
				LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
				hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
				hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
				hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
				hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
				hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
				//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
				//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
				//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;
				
				*/
	
				hp->halo_Mass_hiRes_p3[iBin_hiRes] += pp3->at(i3).mass; 
				break;			

		       	}

		}

		for(iBin_hiRes2 = 0;iBin_hiRes2<250;iBin_hiRes2++){

			if((rad >= iBin_hiRes2*0.1) and (rad < (iBin_hiRes2+1.)*0.1)){
	
				/*
				M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
				mfac =  M_hern*(pp4->at(i4).mass);
				mu =  mfac/(M_hern+(pp4->at(i4).mass));		
				En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
				Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
				Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
				Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
				LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
				hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
				hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
				hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
				hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
				hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
				//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
				//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
				//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;
				*/

				hp->halo_Mass_hiRes2_p3[iBin_hiRes2] += pp3->at(i3).mass; 
				break;			

			}

		}

		for(iBin_hiRes0 = 0;iBin_hiRes0<25;iBin_hiRes0++){

			if((rad >= iBin_hiRes0*1.0) and (rad < (iBin_hiRes0+1.)*1.0)){
	
				/*
				M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
				mfac =  M_hern*(pp4->at(i4).mass);
				mu =  mfac/(M_hern+(pp4->at(i4).mass));		
				En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
				Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
				Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
				Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
				LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
				hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
				hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
				hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
				hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
				hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
				//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
				//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
				//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;
				*/

				hp->halo_Mass_hiRes0_p3[iBin_hiRes0] += pp3->at(i3).mass; 
				break;			

			}

		}


	}

	//Baryons 

	for (int nSat=0;nSat<hp->sat_numb;nSat++){

		int nStartDark, nStopDark, nStartStar, nStopStar;
	
		if(nSat == 0){

			nStartDark = 0;
			nStopDark = hp->myNumbDark[0];	
			nStartStar = 0;	
			nStopStar = hp->myNumbStar[0];
			
		}else{

			nStartDark += hp->myNumbDark[nSat-1];
			nStopDark += hp->myNumbDark[nSat];
                        nStartStar += hp->myNumbStar[nSat-1];	
                        nStopStar += hp->myNumbStar[nSat];

		}

		for(i4=nStartDark;i4<nStopDark;i4++){

		  //if((pp4->at(i4).flag_zone == 3)and (pp4->at(i4).flag_bound == 0)){
		  	if(pp4->at(i4).flag_bound == 0){

				rad = std::sqrt( std::pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2) + std::pow(pp4->at(i4).pos[2],2));
				radxy = std::sqrt( std::pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2));
				cosTheta = (pp4->at(i4).pos[2])/rad;	
				sinTheta = radxy/rad;
				cosPhiAng = (pp4->at(i4).pos[0])/radxy;			
				sinPhiAng = (pp4->at(i4).pos[1])/radxy;			
				vRad = (pp4->at(i4).vel[0])*sinTheta*cosPhiAng + (pp4->at(i4).vel[1])*sinTheta*sinPhiAng + (pp4->at(i4).vel[2])*cosTheta;
				vTheta = (pp4->at(i4).vel[0])*cosTheta*sinPhiAng + (pp4->at(i4).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp4->at(i4).vel[2]); 
				vPhi = -1.*(pp4->at(i4).vel[0])*sinPhiAng + (pp4->at(i4).vel[1])*cosPhiAng;	
				vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);		
				for(iBin = 0;iBin<40;iBin++){

					if((rad >= iBin*6.07) and (rad < (iBin+1.)*6.07)){
	
						M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
						mfac =  M_hern*(pp4->at(i4).mass);
						mu =  mfac/(M_hern+(pp4->at(i4).mass));		
						En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
						Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
						Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
						Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
						LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
						hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
						hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
						hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
						hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
						hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
						//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
						//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
						//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;	
						hp->halo_Mass_p4[iBin+40*nSat] += pp4->at(i4).mass; 
 						break;			

					}

				}

				if (rad > 242.8){

					hp->mass_out_of_242_p4 += pp4->at(i4).mass;
						
				}

				for(iBin_hiRes = 0;iBin_hiRes<50;iBin_hiRes++){

					if((rad >= iBin_hiRes*0.5) and (rad < (iBin_hiRes+1.)*0.5)){
	
						/*
						M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
						mfac =  M_hern*(pp4->at(i4).mass);
						mu =  mfac/(M_hern+(pp4->at(i4).mass));		
						En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
						Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
						Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
						Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
						LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
						hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
						hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
						hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
						hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
						hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
						//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
						//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
						//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;

						*/
	
						hp->halo_Mass_hiRes_p4[iBin_hiRes+50*nSat] += pp4->at(i4).mass; 
 						break;			

					}

				}


				for(iBin_hiRes2 = 0;iBin_hiRes2<250;iBin_hiRes2++){

					if((rad >= iBin_hiRes2*0.1) and (rad < (iBin_hiRes2+1.)*0.1)){
	
						/*
						M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
						mfac =  M_hern*(pp4->at(i4).mass);
						mu =  mfac/(M_hern+(pp4->at(i4).mass));		
						En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
						Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
						Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
						Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
						LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
						hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
						hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
						hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
						hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
						hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
						//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
						//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
						//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;

						*/
	
						hp->halo_Mass_hiRes2_p4[iBin_hiRes2+250*nSat] += pp4->at(i4).mass; 
 						break;			

					}

				}

				for(iBin_hiRes0 = 0;iBin_hiRes0<25;iBin_hiRes0++){
	
					if((rad >= iBin_hiRes0*1.0) and (rad < (iBin_hiRes0+1.)*1.0)){
		
						/*
						M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
						mfac =  M_hern*(pp4->at(i4).mass);
						mu =  mfac/(M_hern+(pp4->at(i4).mass));		
						En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
						Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
						Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
						Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
						LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
						hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
						hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
						hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
						hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
						hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
						//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
						//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
						//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;
						*/
		
						hp->halo_Mass_hiRes0_p4[iBin_hiRes0+25*nSat] += pp4->at(i4).mass; 
						break;			
		
					}
	
				}
		
			}

		}	 

		//Dark 

		if(hp->wantBaryons == 1){

			for(i5=nStartStar;i5< nStopStar;i5++){

				//if((pp5->at(i5).flag_zone == 3)and (pp5->at(i5).flag_bound == 0)){
				if(pp5->at(i5).flag_bound == 0){
					rad = std::sqrt( std::pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2) + std::pow(pp5->at(i5).pos[2],2));
					radxy = std::sqrt( std::pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2));
					cosTheta = (pp5->at(i5).pos[2])/rad;	
					sinTheta = radxy/rad;
					cosPhiAng = (pp5->at(i5).pos[0])/radxy;			
					sinPhiAng = (pp5->at(i5).pos[1])/radxy;			
					vRad = (pp5->at(i5).vel[0])*sinTheta*cosPhiAng + (pp5->at(i5).vel[1])*sinTheta*sinPhiAng + (pp5->at(i5).vel[2])*cosTheta;
					vTheta = (pp5->at(i5).vel[0])*cosTheta*sinPhiAng + (pp5->at(i5).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp5->at(i5).vel[2]); 				
					vPhi = -1.*(pp5->at(i5).vel[0])*sinPhiAng + (pp5->at(i5).vel[1])*cosPhiAng;	
					vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);	
					for(iBin = 0;iBin<40;iBin++){
						if((rad >= iBin*6.07) and (rad < (iBin+1.)*6.07)){
							M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
							mfac =  M_hern*(pp5->at(i5).mass);
							mu =  mfac/(M_hern+(pp5->at(i5).mass));		
							En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
							Lang_x = ( (pp5->at(i5).pos[1])*(pp5->at(i5).vel[2]) - (pp5->at(i5).pos[2])*(pp5->at(i5).vel[1]) )*mu;
							Lang_y = ( (pp5->at(i5).pos[2])*(pp5->at(i5).vel[0]) - (pp5->at(i5).pos[0])*(pp5->at(i5).vel[2]) )*mu;
							Lang_z = ( (pp5->at(i5).pos[0])*(pp5->at(i5).vel[1]) - (pp5->at(i5).pos[1])*(pp5->at(i5).vel[0]) )*mu;	
							LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
							hp->halo_Vrad_p5[iBin] += vRad*(pp5->at(i5).mass);		
							hp->halo_Vtheta_p5[iBin] += vTheta*(pp5->at(i5).mass);
							hp->halo_Vphi_p5[iBin] += vPhi*(pp5->at(i5).mass);	
							hp->halo_ecc_p5[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp5->at(i5).mass);
							hp->smj_ax_p5[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p5[iBin])*(hp->halo_ecc_p5[iBin]) )*G_const*M_hern )*(pp5->at(i5).mass); 
							hp->halo_Mass_p5[iBin+40*nSat] += pp5->at(i5).mass;
							break;			

						}

					}

					if (rad > 242.8){

						hp->mass_out_of_242_p5 += pp5->at(i5).mass;
					
					}
	
					for(iBin_hiRes = 0;iBin_hiRes<50;iBin_hiRes++){

						if((rad >= iBin_hiRes*0.5) and (rad < (iBin_hiRes+1.)*0.5)){
	
							/*
							M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
							mfac =  M_hern*(pp4->at(i4).mass);
							mu =  mfac/(M_hern+(pp4->at(i4).mass));		
							En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
							Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
							Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
							Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
							LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
							hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
							hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
							hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
							hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
							hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
							//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
							//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
							//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;
	
							*/
	
							hp->halo_Mass_hiRes_p5[iBin_hiRes+50*nSat] += pp5->at(i5).mass; 
	 						break;			

						}

					}


					for(iBin_hiRes2 = 0;iBin_hiRes2<250;iBin_hiRes2++){

						if((rad >= iBin_hiRes2*0.1) and (rad < (iBin_hiRes2+1.)*0.1)){
	
							/*
							M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
							mfac =  M_hern*(pp4->at(i4).mass);
							mu =  mfac/(M_hern+(pp4->at(i4).mass));		
							En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
							Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
							Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
							Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
							LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
							hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
							hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
							hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
							hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
							hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
							//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
							//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
							//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;
	
							*/
		
							hp->halo_Mass_hiRes2_p5[iBin_hiRes2+250*nSat] += pp5->at(i5).mass; 
	 						break;			
	
						}
	
					}

					for(iBin_hiRes0 = 0;iBin_hiRes0<25;iBin_hiRes0++){
	
						if((rad >= iBin_hiRes0*1.0) and (rad < (iBin_hiRes0+1.)*1.0)){
		
							/*
							M_hern = M_tot*(rad*rad)/((rad+a_halo)*(rad+a_halo));
							mfac =  M_hern*(pp4->at(i4).mass);
							mu =  mfac/(M_hern+(pp4->at(i4).mass));		
							En = -1.*G_const*mfac/(rad+a_halo) + 1./2.*mu*vSQ; 
							Lang_x = ( (pp4->at(i4).pos[1])*(pp4->at(i4).vel[2]) - (pp4->at(i4).pos[2])*(pp4->at(i4).vel[1]) )*mu;
							Lang_y = ( (pp4->at(i4).pos[2])*(pp4->at(i4).vel[0]) - (pp4->at(i4).pos[0])*(pp4->at(i4).vel[2]) )*mu;
							Lang_z = ( (pp4->at(i4).pos[0])*(pp4->at(i4).vel[1]) - (pp4->at(i4).pos[1])*(pp4->at(i4).vel[0]) )*mu;	
							LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
							hp->halo_Vrad_p4[iBin] += vRad*(pp4->at(i4).mass);		
							hp->halo_Vtheta_p4[iBin] += vTheta*(pp4->at(i4).mass);
							hp->halo_Vphi_p4[iBin] += vPhi*(pp4->at(i4).mass);	
							hp->halo_ecc_p4[iBin] += std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) )*(pp4->at(i4).mass);
							hp->smj_ax_p4[iBin] += std::pow((LL/mu),2)*1./( (1.-(hp->halo_ecc_p4[iBin])*(hp->halo_ecc_p4[iBin]) )*G_const*M_hern )*(pp4->at(i4).mass); 
							//hp->halo_Lx_p4[iBin] += Lang_x*(pp4->at(i4).mass)/mu;		
							//hp->halo_Ly_p4[iBin] += Lang_y*(pp4->at(i4).mass)/mu;
							//hp->halo_Lz_p4[iBin] += Lang_z*(pp4->at(i4).mass)/mu;
							*/
		
							hp->halo_Mass_hiRes0_p5[iBin_hiRes0+25*nSat] += pp5->at(i5).mass; 
							break;			
		
						}
	
					}
		
				}

			}

		}

	}

	//Final divisions		
	
	for(iBin = 0; iBin<40;iBin++){

		//Final calculations	
		if(hp->halo_Mass_p3[iBin] != 0){
			hp->halo_Vrad_p3[iBin] /=   hp->halo_Mass_p3[iBin];
			hp->halo_Vtheta_p3[iBin] /= hp->halo_Mass_p3[iBin];
			hp->halo_Vphi_p3[iBin] /= hp->halo_Mass_p3[iBin];
			hp->halo_ecc_p3[iBin] /= hp->halo_Mass_p3[iBin];
			hp->smj_ax_p3[iBin] /= hp->halo_Mass_p3[iBin];
										
		}else{				
			hp->halo_Vrad_p3[iBin] = 0.;
			hp->halo_Vtheta_p3[iBin] = 0.;
			hp->halo_Vphi_p3[iBin] = 0.;
			hp->halo_ecc_p3[iBin] = 0.;
			hp->smj_ax_p3[iBin] = 0.;
		}

		//part 4	
		if(hp->halo_Mass_p4[iBin] != 0){
			hp->halo_Vrad_p4[iBin] /= hp->halo_Mass_p4[iBin];
			hp->halo_Vtheta_p4[iBin] /= hp->halo_Mass_p4[iBin];
			hp->halo_Vphi_p4[iBin] /= hp->halo_Mass_p4[iBin];
			hp->halo_ecc_p4[iBin] /= hp->halo_Mass_p4[iBin];
			hp->smj_ax_p4[iBin] /= hp->halo_Mass_p4[iBin];
										
		}else{				
			hp->halo_Vrad_p4[iBin] = 0.;
			hp->halo_Vtheta_p4[iBin] = 0.;
			hp->halo_Vphi_p4[iBin] = 0.;
			hp->halo_ecc_p4[iBin] = 0.;
			hp->smj_ax_p4[iBin] = 0.;
		}
		
		if(hp->wantBaryons == 1){

			//part 5	
			if(hp->halo_Mass_p5[iBin] != 0){
				hp->halo_Vrad_p5[iBin] /= hp->halo_Mass_p5[iBin];
				hp->halo_Vtheta_p5[iBin] /= hp->halo_Mass_p5[iBin];
				hp->halo_Vphi_p5[iBin] /= hp->halo_Mass_p5[iBin];
				hp->halo_ecc_p5[iBin] /= hp->halo_Mass_p5[iBin];
				hp->smj_ax_p5[iBin] /= hp->halo_Mass_p5[iBin];
										
			}else{				

				hp->halo_Vrad_p5[iBin] = 0.;
				hp->halo_Vtheta_p5[iBin] = 0.;
				hp->halo_Vphi_p5[iBin] = 0.;
				hp->halo_ecc_p5[iBin] = 0.;
				hp->smj_ax_p5[iBin] = 0.;

			}

		}

	}

	for(i3=0;i3< hp->nbody3;i3++){


		rad = std::sqrt( std::pow(pp3->at(i3).pos[0],2) + std::pow(pp3->at(i3).pos[1],2) + std::pow(pp3->at(i3).pos[2],2));
		radxy = std::sqrt( std::pow(pp3->at(i3).pos[0],2) + std::pow(pp3->at(i3).pos[1],2));
		cosTheta = (pp3->at(i3).pos[2])/rad;	
		sinTheta = radxy/rad;
		cosPhiAng = (pp3->at(i3).pos[0])/radxy;			
		sinPhiAng = (pp3->at(i3).pos[1])/radxy;			
		vRad = (pp3->at(i3).vel[0])*sinTheta*cosPhiAng + (pp3->at(i3).vel[1])*sinTheta*sinPhiAng + (pp3->at(i3).vel[2])*cosTheta;		
		vTheta = (pp3->at(i3).vel[0])*cosTheta*sinPhiAng + (pp3->at(i3).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp3->at(i3).vel[2]); 		
		vPhi = -1.*(pp3->at(i3).vel[0])*sinPhiAng + (pp3->at(i3).vel[1])*cosPhiAng;	
		vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);	
		for(iBin = 0;iBin<40;iBin++){
			if((rad >= iBin*6.07) and (rad < (iBin+1.)*6.07)){
				hp->halo_SigVrad_p3[iBin] += std::pow((vRad - hp->halo_Vrad_p3[iBin]),2)*(pp3->at(i3).mass);		
				hp->halo_SigVtheta_p3[iBin] += std::pow((vTheta - hp->halo_Vtheta_p3[iBin]),2)*(pp3->at(i3).mass);
				hp->halo_SigVphi_p3[iBin] += std::pow((vPhi - hp->halo_Vphi_p3[iBin]),2)*(pp3->at(i3).mass);	
				break;			
			}
		}

	}


	//Baryons 

	for(i4=0;i4< hp->nbody4;i4++){

		if((pp4->at(i4).flag_zone == 3)and (pp4->at(i4).flag_bound == 0)){
	
			rad = std::sqrt( std::pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2) + std::pow(pp4->at(i4).pos[2],2));
			radxy = std::sqrt( std::pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2));
			cosTheta = (pp4->at(i4).pos[2])/rad;	
			sinTheta = radxy/rad;
			cosPhiAng = (pp4->at(i4).pos[0])/radxy;			
			sinPhiAng = (pp4->at(i4).pos[1])/radxy;			
			vRad = (pp4->at(i4).vel[0])*sinTheta*cosPhiAng + (pp4->at(i4).vel[1])*sinTheta*sinPhiAng + (pp4->at(i4).vel[2])*cosTheta;	
			vTheta = (pp4->at(i4).vel[0])*cosTheta*sinPhiAng + (pp4->at(i4).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp4->at(i4).vel[2]); 
			vPhi = -1.*(pp4->at(i4).vel[0])*sinPhiAng + (pp4->at(i4).vel[1])*cosPhiAng;	
			vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);	
			for(iBin = 0;iBin<40;iBin++){

				if((rad >= iBin*6.07) and (rad < (iBin+1.)*6.07)){
					hp->halo_SigVrad_p4[iBin] += std::pow((vRad - hp->halo_Vrad_p4[iBin]),2)*(pp4->at(i4).mass);		
					hp->halo_SigVtheta_p4[iBin] += std::pow((vTheta - hp->halo_Vtheta_p4[iBin]),2)*(pp4->at(i4).mass);
					hp->halo_SigVphi_p4[iBin] += std::pow((vPhi - hp->halo_Vphi_p4[iBin]),2)*(pp4->at(i4).mass);	
					break;			
				}
			}
		}

	}

	//Dark 

	if(hp->wantBaryons == 1){

		for(i5=0;i5< hp->nbody5;i5++){

			if(( pp5->at(i5).flag_zone == 3)and (pp5->at(i5).flag_bound == 0)){
				rad = std::sqrt( std::pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2) + std::pow(pp5->at(i5).pos[2],2));
				radxy = std::sqrt( std::pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2));
				cosTheta = (pp5->at(i5).pos[2])/rad;	
				sinTheta = radxy/rad;
				cosPhiAng = (pp5->at(i5).pos[0])/radxy;			
				sinPhiAng = (pp5->at(i5).pos[1])/radxy;			
				vRad = (pp5->at(i5).vel[0])*sinTheta*cosPhiAng + (pp5->at(i5).vel[1])*sinTheta*sinPhiAng + (pp5->at(i5).vel[2])*cosTheta;
				vTheta = (pp5->at(i5).vel[0])*cosTheta*sinPhiAng + (pp5->at(i5).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp5->at(i5).vel[2]); 				
				vPhi = -1.*(pp5->at(i5).vel[0])*sinPhiAng + (pp5->at(i5).vel[1])*cosPhiAng;	
				vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);	
				for(iBin = 0;iBin<40;iBin++){

					if((rad >= iBin*6.07) and (rad < (iBin+1.)*6.07)){
						hp->halo_SigVrad_p5[iBin] += std::pow((vRad - hp->halo_Vrad_p5[iBin]),2)*(pp5->at(i5).mass);		
						hp->halo_SigVtheta_p5[iBin] += std::pow((vTheta - hp->halo_Vtheta_p5[iBin]),2)*(pp5->at(i5).mass);
						hp->halo_SigVphi_p5[iBin] += std::pow((vPhi - hp->halo_Vphi_p5[iBin]),2)*(pp5->at(i5).mass);	
						break;			
					}
				}
			}
		}
	}


	//Final divisions		
	
	for(iBin = 0; iBin<40;iBin++){

		//Final calculation halo vRad	
		if(hp->halo_Mass_p3[iBin] != 0){
			hp->halo_SigVrad_p3[iBin] /= hp->halo_Mass_p3[iBin];
			hp->halo_SigVtheta_p3[iBin] /= hp->halo_Mass_p3[iBin];
			hp->halo_SigVphi_p3[iBin] /= hp->halo_Mass_p3[iBin];			
												
		}else{	
				
			hp->halo_SigVrad_p3[iBin] = 0.;
			hp->halo_SigVtheta_p3[iBin] = 0.;
			hp->halo_SigVphi_p3[iBin] = 0.;			

		}
	
		//part 4	
		if(hp->halo_Mass_p4[iBin] != 0){

			hp->halo_SigVrad_p4[iBin] /= hp->halo_Mass_p4[iBin];
			hp->halo_SigVtheta_p4[iBin] /= hp->halo_Mass_p4[iBin];	
			hp->halo_SigVphi_p4[iBin] /= hp->halo_Mass_p4[iBin];
				
		}else{	
			
			hp->halo_SigVrad_p4[iBin] = 0.;
			hp->halo_SigVtheta_p4[iBin] = 0.;
			hp->halo_SigVphi_p4[iBin] = 0.;
	
		}

		if(hp->wantBaryons == 1){

			//Final calculation dark vRad	
			if(hp->halo_Mass_p5[iBin] != 0){

				hp->halo_SigVrad_p5[iBin] /= hp->halo_Mass_p5[iBin];
				hp->halo_SigVtheta_p5[iBin] /= hp->halo_Mass_p5[iBin];
				hp->halo_SigVphi_p5[iBin] /= hp->halo_Mass_p5[iBin];
									
			}else{				

				hp->halo_SigVrad_p5[iBin] = 0.;
				hp->halo_SigVtheta_p5[iBin] = 0.;
				hp->halo_SigVphi_p5[iBin] = 0.;

			}
	
		}

	}

	float myMass = 0.;
	for(int i = 0; i< hp->nbody4;i++){
		if(pp4->at(i).flag_bound == 1){
			myMass += pp4->at(i).mass;
		}
	}
	
	if(hp->wantBaryons == 1){
		for(int i = 0; i< hp->nbody5;i++){
			if(pp5->at(i).flag_bound == 1){
				myMass += pp5->at(i).mass;
			}
		}	
	}


	//float radDC = std::sqrt( std::pow((dcms_Sat->pos[0]-dcms->pos[0]),2) + std::pow((dcms_Sat->pos[1]-dcms->pos[1]),2) + std::pow((dcms_Sat->pos[2]-dcms->pos[2]),2) );
	//float vsqDC = std::sqrt( std::pow((dcms_Sat->vel[0]-dcms->vel[0]),2) + std::pow((dcms_Sat->vel[1]-dcms->vel[1]),2) + std::pow((dcms_Sat->vel[2]-dcms->vel[2]),2) );	
	//M_hern = M_tot*(radDC*radDC)/((radDC+a_halo)*(radDC+a_halo));
	//mfac =  M_hern*myMass;
	//mu =  mfac/(M_hern+myMass);		
	//En = -1.*G_const*mfac/(radDC+a_halo) + 1./2.*mu*vsqDC*vsqDC; 
	//Lang_x = ( (dcms_Sat->pos[1]-dcms->pos[1])*(dcms_Sat->vel[2]-dcms->vel[2]) - (dcms_Sat->pos[2]-dcms->pos[2])*(dcms_Sat->vel[1]-dcms->vel[1]) )*mu;
	//Lang_y = ( (dcms_Sat->pos[2]-dcms->pos[2])*(dcms_Sat->vel[0]-dcms->vel[0]) - (dcms_Sat->pos[0]-dcms->pos[0])*(dcms_Sat->vel[2]-dcms->vel[2]) )*mu;
	//Lang_z = ( (dcms_Sat->pos[0]-dcms->pos[0])*(dcms_Sat->vel[1]-dcms->vel[1]) - (dcms_Sat->pos[1]-dcms->pos[1])*(dcms_Sat->vel[0]-dcms->vel[0]) )*mu;
	//LL = std::sqrt( (Lang_x)*(Lang_x) + (Lang_y)*(Lang_y) + (Lang_z)*(Lang_z) );
	//hp->sat_ecc[] = std::sqrt( 1. + 2.*En*LL*LL/( mu*( G_const*G_const*mfac*mfac ) ) );
	//hp->sat_ax  = std::pow((LL/mu),2)*1./( (1.-(hp->sat_ecc)*(hp->sat_ecc) )*G_const*M_hern ); 

}

void calc_tensor_inertia(std::vector<struct h5_particle> *pp3,std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, std::string filenumb){

	/* Write files for Peter */

	std::string newfilenumb;

	if(atoi((filenumb).c_str()) < 10){
		newfilenumb = "00"+filenumb;
	}else if(atoi((filenumb).c_str()) >= 10){
		newfilenumb = "0"+filenumb;
	}
	

	std::ofstream myFileDebrisStar((hp->myHDF5Root+"_"+newfilenumb+"_starsDebris.txt").c_str());
	std::ofstream myFileDebrisDark((hp->myHDF5Root+"_"+newfilenumb+"_darkDebris.txt").c_str());	

	myFileDebrisDark << 0 << std::endl;
	myFileDebrisDark << 0 << std::endl;
	myFileDebrisDark << std::scientific << std::setprecision(6) << 0 << std::endl;

	myFileDebrisStar << 0 << std::endl;
	myFileDebrisStar << 0 << std::endl;
	myFileDebrisStar << std::scientific << std::setprecision(6) << 0 << std::endl;


	for (int ii=0;ii<50;ii++){

		for(int jjj=0;jjj<3;jjj++){

			for(int kkk=0;kkk<3;kkk++){

				hp->geomTensorDark_Unbound[jjj][kkk][ii] = 0.;
				hp->geomTensorStar_Unbound[jjj][kkk][ii] = 0.;
				hp->geomTensorHalo[jjj][kkk][ii] = 0.;
	
			}

		}

	}


	float rad;

	for (int iP =0;iP<hp->nbody4;iP++){

		rad = std::sqrt( std::pow(pp4->at(iP).pos[0],2) + std::pow(pp4->at(iP).pos[1],2) + std::pow(pp4->at(iP).pos[2],2));

		for(int iR=0;iR<50;iR++){

			if ( ( iR*0.5 < rad )	and ( ( iR + 1. )*0.5 >= rad ) ){

				hp->inertiaTensorDark[0][0][iR] += pp4->at(iP).mass*( std::pow(pp4->at(iP).pos[1],2) + std::pow(pp4->at(iP).pos[2],2) );
				hp->inertiaTensorDark[1][1][iR] += pp4->at(iP).mass*( std::pow(pp4->at(iP).pos[2],2) + std::pow(pp4->at(iP).pos[0],2) );
				hp->inertiaTensorDark[2][2][iR] += pp4->at(iP).mass*( std::pow(pp4->at(iP).pos[0],2) + std::pow(pp4->at(iP).pos[1],2) );
				hp->inertiaTensorDark[0][1][iR] -= pp4->at(iP).mass*( pp4->at(iP).pos[0]*pp4->at(iP).pos[1] );
				hp->inertiaTensorDark[0][2][iR] -= pp4->at(iP).mass*( pp4->at(iP).pos[0]*pp4->at(iP).pos[2] );
				hp->inertiaTensorDark[1][2][iR] -= pp4->at(iP).mass*( pp4->at(iP).pos[1]*pp4->at(iP).pos[2] );

				if(pp4->at(iP).flag_bound == 0){

					hp->inertiaTensorDark_Unbound[0][0][iR] += pp4->at(iP).mass*( std::pow(pp4->at(iP).pos[1],2) + std::pow(pp4->at(iP).pos[2],2) );
					hp->inertiaTensorDark_Unbound[1][1][iR] += pp4->at(iP).mass*( std::pow(pp4->at(iP).pos[2],2) + std::pow(pp4->at(iP).pos[0],2) );
					hp->inertiaTensorDark_Unbound[2][2][iR] += pp4->at(iP).mass*( std::pow(pp4->at(iP).pos[0],2) + std::pow(pp4->at(iP).pos[1],2) );
					hp->inertiaTensorDark_Unbound[0][1][iR] -= pp4->at(iP).mass*( pp4->at(iP).pos[0]*pp4->at(iP).pos[1] );
					hp->inertiaTensorDark_Unbound[0][2][iR] -= pp4->at(iP).mass*( pp4->at(iP).pos[0]*pp4->at(iP).pos[2] );
					hp->inertiaTensorDark_Unbound[1][2][iR] -= pp4->at(iP).mass*( pp4->at(iP).pos[1]*pp4->at(iP).pos[2] );

					hp->geomTensorDark_Unbound[0][0][iR] += pp4->at(iP).mass*( std::pow(pp4->at(iP).pos[0],2) );
					hp->geomTensorDark_Unbound[1][1][iR] += pp4->at(iP).mass*( std::pow(pp4->at(iP).pos[1],2) );
					hp->geomTensorDark_Unbound[2][2][iR] += pp4->at(iP).mass*( std::pow(pp4->at(iP).pos[2],2) );
					hp->geomTensorDark_Unbound[0][1][iR] += pp4->at(iP).mass*( pp4->at(iP).pos[0]*pp4->at(iP).pos[1] );
					hp->geomTensorDark_Unbound[0][2][iR] += pp4->at(iP).mass*( pp4->at(iP).pos[0]*pp4->at(iP).pos[2] );
					hp->geomTensorDark_Unbound[1][2][iR] += pp4->at(iP).mass*( pp4->at(iP).pos[1]*pp4->at(iP).pos[2] );

					myFileDebrisDark << std::scientific << std::setprecision(6) << pp4->at(iP).idf << "\t" << pp4->at(iP).mass << "\t"  << pp4->at(iP).pos[0] << "\t"  << pp4->at(iP).pos[1] << "\t" << pp4->at(iP).pos[2] << "\t" << pp4->at(iP).vel[0] << "\t" << pp4->at(iP).vel[1] << "\t" << pp4->at(iP).vel[2] << std::endl;


				}

				break;

			}

		}

	}

	myFileDebrisDark.close();

	for (int iP =0;iP<hp->nbody5;iP++){

		rad = std::sqrt( std::pow(pp5->at(iP).pos[0],2) + std::pow(pp5->at(iP).pos[1],2) + std::pow(pp5->at(iP).pos[2],2));

		for(int iR=0;iR<50;iR++){

			if ( ( iR*0.5 < rad )	and ( ( iR + 1. )*0.5 >= rad ) ){

				hp->inertiaTensorStar[0][0][iR] += pp5->at(iP).mass*( std::pow(pp5->at(iP).pos[1],2) + std::pow(pp5->at(iP).pos[2],2) );
				hp->inertiaTensorStar[1][1][iR] += pp5->at(iP).mass*( std::pow(pp5->at(iP).pos[2],2) + std::pow(pp5->at(iP).pos[0],2) );
				hp->inertiaTensorStar[2][2][iR] += pp5->at(iP).mass*( std::pow(pp5->at(iP).pos[0],2) + std::pow(pp5->at(iP).pos[1],2) );
				hp->inertiaTensorStar[0][1][iR] -= pp5->at(iP).mass*( pp5->at(iP).pos[0]*pp5->at(iP).pos[1] );
				hp->inertiaTensorStar[0][2][iR] -= pp5->at(iP).mass*( pp5->at(iP).pos[0]*pp5->at(iP).pos[2] );
				hp->inertiaTensorStar[1][2][iR] -= pp5->at(iP).mass*( pp5->at(iP).pos[1]*pp5->at(iP).pos[2] );

				if(pp5->at(iP).flag_bound == 0){

					hp->inertiaTensorStar_Unbound[0][0][iR] += pp5->at(iP).mass*( std::pow(pp5->at(iP).pos[1],2) + std::pow(pp5->at(iP).pos[2],2) );
					hp->inertiaTensorStar_Unbound[1][1][iR] += pp5->at(iP).mass*( std::pow(pp5->at(iP).pos[2],2) + std::pow(pp5->at(iP).pos[0],2) );
					hp->inertiaTensorStar_Unbound[2][2][iR] += pp5->at(iP).mass*( std::pow(pp5->at(iP).pos[0],2) + std::pow(pp5->at(iP).pos[1],2) );
					hp->inertiaTensorStar_Unbound[0][1][iR] -= pp5->at(iP).mass*( pp5->at(iP).pos[0]*pp5->at(iP).pos[1] );
					hp->inertiaTensorStar_Unbound[0][2][iR] -= pp5->at(iP).mass*( pp5->at(iP).pos[0]*pp5->at(iP).pos[2] );
					hp->inertiaTensorStar_Unbound[1][2][iR] -= pp5->at(iP).mass*( pp5->at(iP).pos[1]*pp5->at(iP).pos[2] );

					hp->geomTensorStar_Unbound[0][0][iR] += pp5->at(iP).mass*( std::pow(pp5->at(iP).pos[0],2) );
					hp->geomTensorStar_Unbound[1][1][iR] += pp5->at(iP).mass*( std::pow(pp5->at(iP).pos[1],2) );
					hp->geomTensorStar_Unbound[2][2][iR] += pp5->at(iP).mass*( std::pow(pp5->at(iP).pos[2],2) );
					hp->geomTensorStar_Unbound[0][1][iR] += pp5->at(iP).mass*( pp5->at(iP).pos[0]*pp5->at(iP).pos[1] );
					hp->geomTensorStar_Unbound[0][2][iR] += pp5->at(iP).mass*( pp5->at(iP).pos[0]*pp5->at(iP).pos[2] );
					hp->geomTensorStar_Unbound[1][2][iR] += pp5->at(iP).mass*( pp5->at(iP).pos[1]*pp5->at(iP).pos[2] );

					myFileDebrisStar << std::scientific << std::setprecision(6) << pp5->at(iP).idf << "\t" << pp5->at(iP).mass << "\t"  << pp5->at(iP).pos[0] << "\t"  << pp5->at(iP).pos[1] << "\t" << pp5->at(iP).pos[2] << "\t" << pp5->at(iP).vel[0] << "\t" << pp5->at(iP).vel[1] << "\t" << pp5->at(iP).vel[2] << std::endl;

				}

				break;

			}

		}

	}

	myFileDebrisStar.close();	

	for (int iP =0;iP<hp->nbody3;iP++){

		rad = std::sqrt( std::pow(pp3->at(iP).pos[0],2) + std::pow(pp3->at(iP).pos[1],2) + std::pow(pp3->at(iP).pos[2],2));

		for(int iR=0;iR<50;iR++){

			if ( ( iR*0.5 < rad )	and ( ( iR + 1. )*0.5 >= rad ) ){

				hp->inertiaTensorHalo[0][0][iR] += pp3->at(iP).mass*( std::pow(pp3->at(iP).pos[1],2) + std::pow(pp3->at(iP).pos[2],2) );
				hp->inertiaTensorHalo[1][1][iR] += pp3->at(iP).mass*( std::pow(pp3->at(iP).pos[2],2) + std::pow(pp3->at(iP).pos[0],2) );
				hp->inertiaTensorHalo[2][2][iR] += pp3->at(iP).mass*( std::pow(pp3->at(iP).pos[0],2) + std::pow(pp3->at(iP).pos[1],2) );
				hp->inertiaTensorHalo[0][1][iR] -= pp3->at(iP).mass*( pp3->at(iP).pos[0]*pp3->at(iP).pos[1] );
				hp->inertiaTensorHalo[0][2][iR] -= pp3->at(iP).mass*( pp3->at(iP).pos[0]*pp3->at(iP).pos[2] );
				hp->inertiaTensorHalo[1][2][iR] -= pp3->at(iP).mass*( pp3->at(iP).pos[1]*pp3->at(iP).pos[2] );

				hp->geomTensorHalo[0][0][iR] += pp3->at(iP).mass*( std::pow(pp3->at(iP).pos[0],2) );
				hp->geomTensorHalo[1][1][iR] += pp3->at(iP).mass*( std::pow(pp3->at(iP).pos[1],2) );
				hp->geomTensorHalo[2][2][iR] += pp3->at(iP).mass*( std::pow(pp3->at(iP).pos[2],2) );
				hp->geomTensorHalo[0][1][iR] += pp3->at(iP).mass*( pp3->at(iP).pos[0]*pp3->at(iP).pos[1] );
				hp->geomTensorHalo[0][2][iR] += pp3->at(iP).mass*( pp3->at(iP).pos[0]*pp3->at(iP).pos[2] );
				hp->geomTensorHalo[1][2][iR] += pp3->at(iP).mass*( pp3->at(iP).pos[1]*pp3->at(iP).pos[2] );


				if(pp3->at(iP).flag_bound == 0){

					hp->inertiaTensorHalo_Unbound[0][0][iR] += pp3->at(iP).mass*( std::pow(pp3->at(iP).pos[1],2) + std::pow(pp3->at(iP).pos[2],2) );
					hp->inertiaTensorHalo_Unbound[1][1][iR] += pp3->at(iP).mass*( std::pow(pp3->at(iP).pos[2],2) + std::pow(pp3->at(iP).pos[0],2) );
					hp->inertiaTensorHalo_Unbound[2][2][iR] += pp3->at(iP).mass*( std::pow(pp3->at(iP).pos[0],2) + std::pow(pp3->at(iP).pos[1],2) );
					hp->inertiaTensorHalo_Unbound[0][1][iR] += pp3->at(iP).mass*( pp3->at(iP).pos[0]*pp3->at(iP).pos[1] );
					hp->inertiaTensorHalo_Unbound[0][2][iR] += pp3->at(iP).mass*( pp3->at(iP).pos[0]*pp3->at(iP).pos[2] );
					hp->inertiaTensorHalo_Unbound[1][2][iR] += pp3->at(iP).mass*( pp3->at(iP).pos[1]*pp3->at(iP).pos[2] );

				}

			}

		}

	}



	/*------------------------------------------------------------------------------------------------------------*/

	for (int iR = 0; iR<50;iR++){

		hp->inertiaTensorDark[1][0][iR] = hp->inertiaTensorDark[0][1][iR];
		hp->inertiaTensorDark[2][0][iR] = hp->inertiaTensorDark[0][2][iR];
		hp->inertiaTensorDark[2][1][iR] = hp->inertiaTensorDark[1][2][iR];
		hp->inertiaTensorDark_Unbound[1][0][iR] = hp->inertiaTensorDark_Unbound[0][1][iR];
		hp->inertiaTensorDark_Unbound[2][0][iR] = hp->inertiaTensorDark_Unbound[0][2][iR];
		hp->inertiaTensorDark_Unbound[2][1][iR] = hp->inertiaTensorDark_Unbound[1][2][iR];
		hp->inertiaTensorStar[1][0][iR] = hp->inertiaTensorStar[0][1][iR];
        	hp->inertiaTensorStar[2][0][iR] = hp->inertiaTensorStar[0][2][iR];
        	hp->inertiaTensorStar[2][1][iR] = hp->inertiaTensorStar[1][2][iR];
		hp->inertiaTensorStar_Unbound[1][0][iR] = hp->inertiaTensorStar_Unbound[0][1][iR];
        	hp->inertiaTensorStar_Unbound[2][0][iR] = hp->inertiaTensorStar_Unbound[0][2][iR];
        	hp->inertiaTensorStar_Unbound[2][1][iR] = hp->inertiaTensorStar_Unbound[1][2][iR];
		hp->inertiaTensorHalo[1][0][iR] = hp->inertiaTensorHalo[0][1][iR];
        	hp->inertiaTensorHalo[2][0][iR] = hp->inertiaTensorHalo[0][2][iR];
        	hp->inertiaTensorHalo[2][1][iR] = hp->inertiaTensorHalo[1][2][iR];
		hp->inertiaTensorHalo_Unbound[1][0][iR] = hp->inertiaTensorHalo_Unbound[0][1][iR];
        	hp->inertiaTensorHalo_Unbound[2][0][iR] = hp->inertiaTensorHalo_Unbound[0][2][iR];
        	hp->inertiaTensorHalo_Unbound[2][1][iR] = hp->inertiaTensorHalo_Unbound[1][2][iR];

		hp->geomTensorDark[1][0][iR] = hp->geomTensorDark[0][1][iR];
		hp->geomTensorDark[2][0][iR] = hp->geomTensorDark[0][2][iR];
		hp->geomTensorDark[2][1][iR] = hp->geomTensorDark[1][2][iR];
		hp->geomTensorDark_Unbound[1][0][iR] = hp->geomTensorDark_Unbound[0][1][iR];
		hp->geomTensorDark_Unbound[2][0][iR] = hp->geomTensorDark_Unbound[0][2][iR];
		hp->geomTensorDark_Unbound[2][1][iR] = hp->geomTensorDark_Unbound[1][2][iR];
		hp->geomTensorStar[1][0][iR] = hp->geomTensorStar[0][1][iR];
        	hp->geomTensorStar[2][0][iR] = hp->geomTensorStar[0][2][iR];
        	hp->geomTensorStar[2][1][iR] = hp->geomTensorStar[1][2][iR];
		hp->geomTensorStar_Unbound[1][0][iR] = hp->geomTensorStar_Unbound[0][1][iR];
        	hp->geomTensorStar_Unbound[2][0][iR] = hp->geomTensorStar_Unbound[0][2][iR];
        	hp->geomTensorStar_Unbound[2][1][iR] = hp->geomTensorStar_Unbound[1][2][iR];
		hp->geomTensorHalo[1][0][iR] = hp->geomTensorHalo[0][1][iR];
        	hp->geomTensorHalo[2][0][iR] = hp->geomTensorHalo[0][2][iR];
        	hp->geomTensorHalo[2][1][iR] = hp->geomTensorHalo[1][2][iR];
	}


/********************************************************************/

	//int initPart4, initPart5;   //old
        //int finPart4, finPart5;     //old

	/*

	for (int ind_sat=0;ind_sat<hp->sat_numb;ind_sat++){

		std::cout << "Calculating tidal mass for satellite " << ind_sat+1 << " within tidal radius " <<  hp->sat_tidRadius[ind_sat] << std::endl; 

		int i4, i5;

		if(ind_sat == 0){

			initPart4 = 0;
			finPart4 = hp->myNumbDark[0];
			initPart5 = 0;
			finPart5 = hp->myNumbStar[0];

		}else{

			initPart4 += hp->myNumbDark[ind_sat-1];
			finPart4 += hp->myNumbDark[ind_sat];
                      	initPart5 += hp->myNumbStar[ind_sat-1];
                        finPart5 += hp->myNumbStar[ind_sat];


		}

		for(i4=initPart4;i4<finPart4;i4++){

			if(std::sqrt(std::pow((pp4->at(i4).pos[0] - dcms_Sat->at(ind_sat).pos[0] + dcms->pos[0] ),2) + std::pow((pp4->at(i4).pos[1] - dcms_Sat->at(ind_sat).pos[1]  + dcms->pos[1] ),2) + std::pow((pp4->at(i4).pos[2] - dcms_Sat->at(ind_sat).pos[2]  + dcms->pos[2] ),2) ) > hp->sat_tidRadius[ind_sat]){

				myFileDebrisDark << std::scientific << std::setprecision(6) << pp4->at(i4).idf << "\t" << pp4->at(i4).mass << "\t"  << pp4->at(i4).pos[0] << "\t"  << pp4->at(i4).pos[1] << "\t" << pp4->at(i4).pos[2] << "\t" << pp4->at(i4).vel[0] << "\t" << pp4->at(i4).vel[1] << "\t" << pp4->at(i4).vel[2] << std::endl;

			}

		}
		
		std::cout << "for satellite " << ind_sat+1 << "there are " << hp->sat_tidalMass[ind_sat] << " unit masses inside and " << hp->sat_outTidMass[ind_sat] << " outside" << std::endl; 


		if(hp->wantBaryons == 1){	

			for(i5=initPart5;i5<finPart5;i5++){

				if(std::sqrt( std::pow((pp5->at(i5).pos[0] - dcms_Sat->at(ind_sat).pos[0]  + dcms->pos[0] ),2) + std::pow((pp5->at(i5).pos[1] - dcms_Sat->at(ind_sat).pos[1]  + dcms->pos[1] ),2) + std::pow((pp5->at(i5).pos[2] - dcms_Sat->at(ind_sat).pos[2]  + dcms->pos[2] ),2) ) > hp->sat_tidRadius[ind_sat]  ){

					myFileDebrisStar << std::scientific << std::setprecision(6) << pp5->at(i5).idf << "\t" << pp5->at(i5).mass << "\t"  << pp5->at(i5).pos[0] << "\t"  << pp5->at(i5).pos[1] << "\t" << pp5->at(i5).pos[2] << "\t" << pp5->at(i5).vel[0] << "\t" << pp5->at(i5).vel[1] << "\t" << pp5->at(i5).vel[2] << std::endl;

				}

			}

		}
	
	}

	*/

/***********************************************************************************************************/

	/* this is again mass radial profile to check with peters data*/

        float myMassProfileDark[50];
        float myMassProfileStar[50];
        float massRad;

	for(int iRad=0;iRad<50;iRad++){

		myMassProfileDark[iRad] = 0.;
		myMassProfileStar[iRad] = 0.;

	}

	for(int iP4=0;iP4<hp->nbody4;iP4++){

		if(pp4->at(iP4).flag_bound == 0){

			massRad = std::sqrt( (pp4->at(iP4).pos[0])*(pp4->at(iP4).pos[0]) + (pp4->at(iP4).pos[1])*(pp4->at(iP4).pos[1]) + (pp4->at(iP4).pos[2])*(pp4->at(iP4).pos[2]) );

			for(int iRad=0;iRad<50;iRad++){

				if( (massRad > iRad*0.5) and (massRad <= (iRad + 1)*0.5) ){

					myMassProfileDark[iRad] += pp4->at(iP4).mass;
					break;			
	
				}

			}
	
		}

	}

	for(int iP5=0;iP5<hp->nbody5;iP5++){

		if(pp5->at(iP5).flag_bound == 0){
	
			massRad = std::sqrt( (pp5->at(iP5).pos[0])*(pp5->at(iP5).pos[0]) + (pp5->at(iP5).pos[1])*(pp5->at(iP5).pos[1]) + (pp5->at(iP5).pos[2])*(pp5->at(iP5).pos[2]) );

			for(int iRad=0;iRad<50;iRad++){

				if( (massRad > iRad*0.5) and (massRad <= (iRad + 1)*0.5) ){

					myMassProfileStar[iRad] += pp5->at(iP5).mass;
					break;			
	
				}

			}

		}

	}


	std::ofstream mProfFile("mass_profiles_040_Inertia.txt");

	for(int iRad=0;iRad<50;iRad++){

		mProfFile << std::scientific << std::setprecision(6) << myMassProfileDark[iRad] << "\t" << myMassProfileStar[iRad] <<std::endl;

	}

	mProfFile.close();

}

void calc_En_Ang_25kpc(std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, std::string myStringHalo, std::string myStringDark, std::string myStringStar){

	float rad, kin, ang_z;
	float en_tot, ang_tot, mass_tot;
	int iEn, iL;

	float myMap_Star[600][400];
	float myMap_Dark[600][400];
	float myMap_Halo[600][400];

	for(int iMapX = 0; iMapX<400;iMapX++){

		for(int iMapY = 0; iMapY<600;iMapY++){

			myMap_Star[iMapY][iMapX] = 0.;
			myMap_Dark[iMapY][iMapX] = 0.;
			myMap_Halo[iMapY][iMapX] = 0.;			
 
		}

	}

	en_tot = 0.;
	ang_tot = 0.;
	mass_tot = 0.;

	std::ofstream myOutEnAngHalo(myStringHalo.c_str());
	std::ofstream myOutEnAngHalo_tot(("total_"+myStringHalo).c_str());
	std::ofstream myOutEnAngHalo_map(("map_"+myStringHalo).c_str());
	std::ofstream myOutEnAngDark(myStringDark.c_str());
	std::ofstream myOutEnAngDark_tot(("total_"+myStringDark).c_str());
	std::ofstream myOutEnAngDark_map(("map_"+myStringDark).c_str());
	std::ofstream myOutEnAngStar(myStringStar.c_str());
	std::ofstream myOutEnAngStar_tot(("total_"+myStringStar).c_str());
	std::ofstream myOutEnAngStar_map(("map_"+myStringStar).c_str());

	for (int iP = 0; iP < hp->nbody3;iP++){
		rad = std::sqrt( std::pow(pp3->at(iP).pos[0],2) + std::pow(pp3->at(iP).pos[1],2) + std::pow(pp3->at(iP).pos[2],2) );	 
		if ( rad < 25. ){
			kin = ( std::pow(pp3->at(iP).vel[0],2) + std::pow(pp3->at(iP).vel[1],2) + std::pow(pp3->at(iP).vel[2],2) );
			ang_z = (pp3->at(iP).pos[0]*pp3->at(iP).vel[1] - pp3->at(iP).pos[1]*pp3->at(iP).vel[0]);
			en_tot += kin+pp3->at(iP).phi/2.;
			ang_tot += ang_z;
			mass_tot += pp3->at(iP).mass;
			myOutEnAngHalo << std::scientific << std::setprecision(6) << kin+pp3->at(iP).phi << "\t" << ang_z << "\t" << pp3->at(iP).mass << std::endl;


				/*
				for(iEn=0;iEn<400;iEn++){
					if( ( -1.*std::pow(10.,(4.2 - iEn*0.04)) < kin+pp3->at(iP).phi ) and ( -1.*std::pow(10.,(4.2 - (iEn+1.)*0.04)) >= kin+pp3->at(iP).phi ) ){
					
							for(iL=0;iL<300;iL++){
								if( ( -1.*std::pow(10.,(10. - iL*0.05)) < ang_z ) and ( -1.*std::pow(10.,(10. - (iL+1.)*0.05) ) >= ang_z ) ){ 
									myMap_Halo[iL][iEn] += pp3->at(iP).mass;
									break;
								}
							}

						
						

					        }else if (ang_z >= 0){	
		
							for(iL=0;iL<300;iL++){
								if( ( std::pow(10.,(-5. + iL*0.05)) < ang_z ) and ( std::pow(10.,(-5. + (iL+1.)*0.05 )) >= ang_z ) ){ 
									myMap_Halo[iL+300][iEn] += pp3->at(iP).mass;
									break;
								}
							}
					        }
							
							

				      		break;	
					}						

				}*/

		}

	}

	myOutEnAngHalo.close();
	
	myOutEnAngHalo_tot << std::scientific << std::setprecision(6) << en_tot	<< "\t" << ang_tot << "\t" << mass_tot << "\n" << "0.000000e+00" << "\t" << "0.000000e+00" << "\t" << "0.000000e+00" << std::endl;  

	myOutEnAngHalo_tot.close();

	//use dark and star particles

	en_tot = 0.;
	ang_tot = 0.;
	mass_tot = 0.;


	for (int iP = 0; iP < hp->nbody4;iP++){
		if (pp4->at(iP).flag_bound == 0){
			rad = std::sqrt( std::pow(pp4->at(iP).pos[0],2) + std::pow(pp4->at(iP).pos[1],2) + std::pow(pp4->at(iP).pos[2],2) );	 
			if ( rad < 25. ){
				kin = ( std::pow(pp4->at(iP).vel[0],2) + std::pow(pp4->at(iP).vel[1],2) + std::pow(pp4->at(iP).vel[2],2) );
				ang_z = (pp4->at(iP).pos[0]*pp4->at(iP).vel[1] - pp4->at(iP).pos[1]*pp4->at(iP).vel[0]);
				en_tot += kin+pp4->at(iP).phi/2.;
				ang_tot += ang_z;
				mass_tot += pp4->at(iP).mass;
				myOutEnAngDark << std::scientific << std::setprecision(6) << kin+pp4->at(iP).phi << "\t" << ang_z << "\t" << pp4->at(iP).mass << std::endl;
				//if(kin+pp4->at(iP).phi < 0){

					for(iEn=0;iEn<400;iEn++){
						if( ( -1e5 + iEn*625.0  < kin+pp4->at(iP).phi ) and ( -1e5 +(iEn + 1.)*625.0  >= kin+pp4->at(iP).phi ) ){
						
							//if (ang_z < 0){	
		
								for(iL=0;iL<600;iL++){
									if( ( -5000 + iL*16.66 < ang_z ) and ( -5000 + (iL+1.)*16.66 >= ang_z ) ){ 
										myMap_Dark[iL][iEn] += pp4->at(iP).mass;
										break;
									}
								}
											
							/*
					      		}else if (ang_z>= 0){	
		
								for(iL=0;iL<300;iL++){
									if( ( std::pow(10.,(-3.7 + iL*0.025)) < ang_z ) and ( std::pow(10.,(-3.7 + (iL+1.)*0.025)) >= ang_z ) ){ 
										myMap_Dark[iL+300][iEn] += pp4->at(iP).mass;
										break;
									}
								}
					        	}	
							*/
				      			break;	
						}							
					}	     
				//}
		
				/*

				if(kin+pp4->at(iP).phi >= 0){

					for(iEn=0;iEn<200;iEn++){
						if( ( 1.*std::pow(10.,(-3.8 + iEn*0.04)) < kin+pp4->at(iP).phi ) and ( 1.*std::pow(10.,(-3.8 + (iEn+1.)*0.04)) >= kin+pp4->at(iP).phi ) ){
						
							if (ang_z < 0){	
		
								for(iL=0;iL<300;iL++){
									if( ( -1.*std::pow(10.,(3.8 - iL*0.025) < ang_z )) and ( -1.*std::pow(10.,(3.8 - (iL+1.)*0.025)) >= ang_z ) ){ 
										myMap_Dark[iL][iEn+200] += pp4->at(iP).mass;
										break;
									}
								}
						      	}else if (ang_z>= 0){	
		
								for(iL=0;iL<300;iL++){
									if( ( std::pow(10.,(-3.7 + iL*0.025)) < ang_z ) and ( std::pow(10.,(-3.7 + (iL+1.)*0.025)) >= ang_z ) ){ 
										myMap_Dark[iL+300][iEn+200] += pp4->at(iP).mass;
										break;
									}
								}
					        	}	
				      			break;	
						}							
					}	     
				}

				*/

			}
		}
	}
	myOutEnAngDark.close();
	
	myOutEnAngDark_tot << std::scientific << std::setprecision(6) << en_tot	<< "\t" << ang_tot << "\t" << mass_tot << "\n" << "0.000000e+00" << "\t" << "0.000000e+00" << "\t" << "0.000000e+00" << std::endl;  

	myOutEnAngDark_tot.close();

	en_tot = 0.;
	ang_tot = 0.;
	mass_tot = 0;
	
	for (int iP = 0; iP < hp->nbody5;iP++){
		if (pp5->at(iP).flag_bound == 0){
			rad = std::sqrt( std::pow(pp5->at(iP).pos[0],2) + std::pow(pp5->at(iP).pos[1],2) + std::pow(pp5->at(iP).pos[2],2) );	 
			if ( rad < 25. ){
				kin = ( std::pow(pp5->at(iP).vel[0],2) + std::pow(pp5->at(iP).vel[1],2) + std::pow(pp5->at(iP).vel[2],2) );
				ang_z = (pp5->at(iP).pos[0]*pp5->at(iP).vel[1] - pp5->at(iP).pos[1]*pp5->at(iP).vel[0]);
				en_tot += kin+pp5->at(iP).phi/2.;
				ang_tot += ang_z;
				mass_tot += pp5->at(iP).mass;
				myOutEnAngStar << std::scientific << std::setprecision(6) << kin+pp5->at(iP).phi << "\t" << ang_z << "\t" << pp5->at(iP).mass <<  std::endl;

				//if(kin+pp5->at(iP).phi < 0){

					for(iEn=0;iEn<400;iEn++){
						if( ( -1e5 + iEn*625.0  < kin+pp5->at(iP).phi ) and ( -1e5 +(iEn + 1.)*625.0  >= kin+pp5->at(iP).phi ) ){
						
							//if (ang_z < 0){	
		
								for(iL=0;iL<600;iL++){
									if( ( -5000 + iL*16.66 < ang_z ) and ( -5000 + (iL+1.)*16.66 >= ang_z ) ){ 
										myMap_Star[iL][iEn] += pp5->at(iP).mass;
										break;
									}
								}
											
							/*
					      		}else if (ang_z>= 0){	
		
								for(iL=0;iL<300;iL++){
									if( ( std::pow(10.,(-3.7 + iL*0.025)) < ang_z ) and ( std::pow(10.,(-3.7 + (iL+1.)*0.025)) >= ang_z ) ){ 
										myMap_Dark[iL+300][iEn] += pp4->at(iP).mass;
										break;
									}
								}
					        	}	
							*/
				      			break;	
						}							
					}
			}
		}
	}

	myOutEnAngStar.close();

	myOutEnAngStar_tot << std::scientific << std::setprecision(6) << en_tot	<< "\t" << ang_tot << "\t" << mass_tot << "\n" << "0.000000e00" << "\t" << "0.000000e00" << "\t" << "0.000000e00" << std::endl;  

	myOutEnAngStar_tot.close();

	myOutEnAngDark_map << std::scientific << std::setprecision(6);
	myOutEnAngStar_map << std::scientific << std::setprecision(6);
	myOutEnAngHalo_map << std::scientific << std::setprecision(6);

	for(iL=0;iL<600;iL++){
		for(iEn=0;iEn<400;iEn++){	

			myOutEnAngStar_map << myMap_Star[iL][iEn] << "\t" ;
			myOutEnAngDark_map << myMap_Dark[iL][iEn] << "\t" ;
			myOutEnAngHalo_map << myMap_Halo[iL][iEn] << "\t" ;
		}
		
		myOutEnAngStar_map << std::endl;
		myOutEnAngDark_map << std::endl;
		myOutEnAngHalo_map << std::endl;
	
	}	

	myOutEnAngStar_map.close();
	myOutEnAngDark_map.close();
	myOutEnAngHalo_map.close();

}



void calc_Toomre_Space(std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, std::string myStringHalo, std::string myStringDark, std::string myStringStar){

	int i5, i4, i3, iPhi, iRz;
	float rad, cosTheta, sinTheta, cosPhiAng, sinPhiAng, vRad, vTheta, vPhi, vSQ, radxy;

	float myMap_Star[400][600];
	float myMap_Dark[400][600];
	float myMap_Halo[400][600];

	for(int iMapX = 0; iMapX<600;iMapX++){

		for(int iMapY = 0; iMapY<400;iMapY++){

			myMap_Star[iMapY][iMapX] = 0.;
			myMap_Dark[iMapY][iMapX] = 0.;
			myMap_Halo[iMapY][iMapX] = 0.; 

		}

	}

	std::ofstream myOutToomreDark(("toomre_"+myStringDark).c_str());
	std::ofstream myOutToomreStar(("toomre_"+myStringStar).c_str());
	std::ofstream myOutToomreHalo(("toomre_"+myStringHalo).c_str());

	std::ofstream myOutToomreStar_map(("map_v_"+myStringStar).c_str());
	std::ofstream myOutToomreDark_map(("map_v_"+myStringDark).c_str());	
	std::ofstream myOutToomreHalo_map(("map_v_"+myStringHalo).c_str());



	//Halo

	for(i3=0;i3< hp->nbody3;i3++){

		rad = std::sqrt( std::pow(pp3->at(i3).pos[0],2) + std::pow(pp3->at(i3).pos[1],2) + std::pow(pp3->at(i3).pos[2],2));
		if(rad <= 25.){
			radxy = std::sqrt( std::pow(pp3->at(i3).pos[0],2) + std::pow(pp3->at(i3).pos[1],2));
 			cosTheta = (pp3->at(i3).pos[2])/rad;	
			sinTheta = radxy/rad;
			cosPhiAng = (pp3->at(i3).pos[0])/radxy;			
			sinPhiAng = (pp3->at(i3).pos[1])/radxy;			
			vRad = (pp3->at(i3).vel[0])*sinTheta*cosPhiAng + (pp3->at(i3).vel[1])*sinTheta*sinPhiAng + (pp3->at(i3).vel[2])*cosTheta;
			vTheta = (pp3->at(i3).vel[0])*cosTheta*sinPhiAng + (pp3->at(i3).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp3->at(i3).vel[2]); 
			vPhi = -1.*(pp3->at(i3).vel[0])*sinPhiAng + (pp3->at(i3).vel[1])*cosPhiAng;	
			vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);	

			myOutToomreHalo << std::scientific << std::setprecision(6) << vRad << "\t" << vPhi << "\t" << pp3->at(i3).vel[2] << "\t" << pp3->at(i3).mass <<  std::endl;

			for(iRz=0;iRz<400;iRz++){
				if( ( std::pow(10.,(0. + iRz*0.0075)) < std::sqrt( std::pow(vRad,2) + std::pow(pp3->at(i3).vel[2],2) ) ) and ( std::pow(10.,(0. + (iRz+1.)*0.0075)) >= std::sqrt( std::pow(vRad,2) + std::pow(pp3->at(i3).vel[2],2) ) ) ){
						
					if (vPhi < 0){	
		
						for(iPhi=0;iPhi<300;iPhi++){
							if( ( -1.*std::pow(10.,(3. - iPhi*0.03)) < vPhi ) and ( -1.*std::pow(10.,(3. - (iPhi+1.)*0.03) ) >= vPhi ) ){ 
								myMap_Halo[iRz][iPhi] += pp3->at(i3).mass;
								break;
							}
						}
				        }else if (vPhi >= 0){	
	
						for(iPhi=0;iPhi<300;iPhi++){
							if( ( std::pow(10.,(-6. + iPhi*0.03)) < vPhi ) and ( std::pow(10.,(-6. + (iPhi+1.)*0.03 )) >= vPhi ) ){ 
								myMap_Halo[iRz][iPhi+300] += pp3->at(i3).mass;
								break;
							}
						}
				        }
			      		break;	
				}						
			}

		}
			
	}


	//Dark 

	for(i4=0;i4< hp->nbody4;i4++){

		if(pp4->at(i4).flag_bound == 0){

			rad = std::sqrt( std::pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2) + std::pow(pp4->at(i4).pos[2],2));

			if(rad <= 25.){

				radxy = std::sqrt( std::pow(pp4->at(i4).pos[0],2) + std::pow(pp4->at(i4).pos[1],2));
				cosTheta = (pp4->at(i4).pos[2])/rad;	
				sinTheta = radxy/rad;
				cosPhiAng = (pp4->at(i4).pos[0])/radxy;			
				sinPhiAng = (pp4->at(i4).pos[1])/radxy;			
				vRad = (pp4->at(i4).vel[0])*sinTheta*cosPhiAng + (pp4->at(i4).vel[1])*sinTheta*sinPhiAng + (pp4->at(i4).vel[2])*cosTheta;
				vTheta = (pp4->at(i4).vel[0])*cosTheta*sinPhiAng + (pp4->at(i4).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp4->at(i4).vel[2]); 
				vPhi = -1.*(pp4->at(i4).vel[0])*sinPhiAng + (pp4->at(i4).vel[1])*cosPhiAng;	
				vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);	

				myOutToomreDark << std::scientific << std::setprecision(6) << vRad << "\t" << vPhi << "\t" << pp4->at(i4).vel[2] << "\t" << pp4->at(i4).mass <<  std::endl;

				for(iRz=0;iRz<400;iRz++){
					if( ( 0. + iRz*1.5 < std::sqrt( std::pow(vRad,2) + std::pow(pp4->at(i4).vel[2],2) ) ) and ( 0. + (iRz+1.)*1.5  >= std::sqrt( std::pow(vRad,2) + std::pow(pp4->at(i4).vel[2],2) ) ) ){
						
						//if (vPhi < 0){	
		
						for(iPhi=0;iPhi<600;iPhi++){
							if( ( (-600. + iPhi*2.) < vPhi ) and ( -600. + (iPhi+1.)*2. >= vPhi ) ){ 
								myMap_Dark[iRz][iPhi] += pp4->at(i4).mass;
								break;
							}
						}
						
					        /*
						}else if (vPhi >= 0){	
		
							for(iPhi=0;iPhi<300;iPhi++){
								if( ( std::pow(10.,(-6. + iPhi*0.03)) < vPhi ) and ( std::pow(10.,(-6. + (iPhi+1.)*0.03 )) >= vPhi ) ){ 
									myMap_Dark[iRz][iPhi+300] += pp4->at(i4).mass;
									break;
								}
							}
					        }
						*/

				      		break;	

					}						
				}
			}		
		}
	}

	//Baryons 

	for(i5=0;i5< hp->nbody5;i5++){

		if(pp5->at(i5).flag_bound == 0){

			rad = std::sqrt( std::pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2) + std::pow(pp5->at(i5).pos[2],2));

			if(rad <= 25.){

				radxy = std::sqrt( std::pow(pp5->at(i5).pos[0],2) + std::pow(pp5->at(i5).pos[1],2));
				cosTheta = (pp5->at(i5).pos[2])/rad;	
				sinTheta = radxy/rad;
				cosPhiAng = (pp5->at(i5).pos[0])/radxy;			
				sinPhiAng = (pp5->at(i5).pos[1])/radxy;			
				vRad = (pp5->at(i5).vel[0])*sinTheta*cosPhiAng + (pp5->at(i5).vel[1])*sinTheta*sinPhiAng + (pp5->at(i5).vel[2])*cosTheta;	
				vTheta = (pp5->at(i5).vel[0])*cosTheta*sinPhiAng + (pp5->at(i5).vel[1])*cosTheta*sinPhiAng - sinTheta*(pp5->at(i5).vel[2]); 
				vPhi = -1.*(pp5->at(i5).vel[0])*sinPhiAng + (pp5->at(i5).vel[1])*cosPhiAng;	
				vSQ = std::pow(vRad,2)+std::pow(vTheta,2)+std::pow(vPhi,2);	

				myOutToomreStar << std::scientific << std::setprecision(6) << vRad << "\t" << vPhi << "\t" << pp5->at(i5).vel[2] << "\t" << pp5->at(i5).mass <<  std::endl;

				for(iRz=0;iRz<400;iRz++){
					if( ( 0. + iRz*1.5 < std::sqrt( std::pow(vRad,2) + std::pow(pp5->at(i5).vel[2],2) ) ) and ( 0. + (iRz+1.)*1.5  >= std::sqrt( std::pow(vRad,2) + std::pow(pp5->at(i5).vel[2],2) ) ) ){


						/*	

						if (vPhi < 0){	
		
						*/	

						
						for(iPhi=0;iPhi<600;iPhi++){
							if( ( -600. + iPhi*2. < vPhi ) and ( -600. + (iPhi+1.)*2. >= vPhi ) ){ 
								myMap_Star[iRz][iPhi] += pp5->at(i5).mass;
								break;
							}
						}

					       /* } else if (vPhi >= 0){	
		
							for(iPhi=0;iPhi<300;iPhi++){
								if( ( std::pow(10.,(-6. + iPhi*0.03)) < vPhi ) and ( std::pow(10.,(-6. + (iPhi+1.)*0.03 )) >= vPhi ) ){ 
									myMap_Star[iRz][iPhi] += pp5->at(i5).mass;
									break;
								}
							}
					        }*/
				      		break;	
					}						
				}
			}
		}
	}	

	myOutToomreDark.close();
	myOutToomreStar.close();
	myOutToomreHalo.close();




	myOutToomreStar_map << std::scientific << std::setprecision(6);
	myOutToomreDark_map << std::scientific << std::setprecision(6);
	myOutToomreHalo_map << std::scientific << std::setprecision(6);

	for(iRz=0;iRz<400;iRz++){
		for(iPhi=0;iPhi<600;iPhi++){	

			myOutToomreStar_map << myMap_Star[iRz][iPhi] << "\t" ;
			myOutToomreDark_map << myMap_Dark[iRz][iPhi] << "\t" ;
			myOutToomreHalo_map << myMap_Halo[iRz][iPhi] << "\t" ;

		}
		
		myOutToomreStar_map << std::endl;
		myOutToomreDark_map << std::endl;
		myOutToomreHalo_map << std::endl;
	
	}	

	myOutToomreStar_map.close();
	myOutToomreDark_map.close();
	myOutToomreHalo_map.close();

}


