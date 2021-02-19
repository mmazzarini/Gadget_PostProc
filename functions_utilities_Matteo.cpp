#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "library_MMazzarini_structs.h"
#include "library_MMazzarini_utilities.h"


std::string myString(int iFile){

		std::string fileNumber;
		std::ostringstream numberStream;   
		numberStream << iFile;      
		fileNumber = numberStream.str();
		return fileNumber;

}

/*		Reset values of variables		*/

void reset_structures(struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, std::vector<struct h5_particle> *pp0, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3,std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, struct array_hd *ap){

	hp->mass_out_of_242_p4 = 0.;
	hp->mass_out_of_242_p5 = 0.;


	for(int s=0;s<3;s++){

		for(int r=0;r<3;r++){

			for(int iR=0;iR<50;iR++){	
			
				hp->inertiaTensorDark[s][r][iR] = 0.;
				hp->inertiaTensorDark_Unbound[s][r][iR] = 0.;
				hp->inertiaTensorStar[s][r][iR] = 0.;
	                        hp->inertiaTensorStar_Unbound[s][r][iR] = 0.;
				hp->inertiaTensorHalo[s][r][iR] = 0.;
	                        hp->inertiaTensorHalo_Unbound[s][r][iR] = 0.;

			}

		}

	}

	//Density centers
	for (int n=0; n<3;n++){

		dcms->pos[n] = 0.;
		dcms->vel[n] = 0.;
		dcms->angM[n] = 0.;

		for (int j=0;j<hp->sat_numb;j++){

			dcms_Sat->at(j).pos[n] = 0.;
			dcms_Sat->at(j).vel[n] = 0.;
			dcms_Sat->at(j).angM[n] = 0.;

		}

		/*

		dcms_Baryon->pos[n] = 0.;
		dcms_Baryon->vel[n] = 0.;
		dcms_Baryon->angM[n] = 0.;
		dcms_Dark->pos[n] = 0.;
		dcms_Dark->vel[n] = 0.;
		dcms_Dark->angM[n] = 0.;
	
		*/	

	}

	for (int j=0;j<hp->sat_numb;j++){

		dcms_Sat->at(j).mass = 0.;

	}

	dcms->mass = 0.;
	//dcms_Baryon->mass = 0.;
	//dcms_Dark->mass = 0.;

	//particle reset

	pp0->resize(0);
	pp1->resize(0);
	pp2->resize(0);
	pp3->resize(0);
	pp4->resize(0);
	
	if(hp->wantBaryons == 1){
	
		pp5->resize(0);

	}

	//particle header reset	

	hp->myHDF5File = "";	
	hp->myHDF5Post = "";
	hp->myHDF5Maps = "";
	hp->nbodies = 0;
	hp->nbody0 = 0;
	hp->nbody1 = 0;
	hp->nbody2 = 0;
	hp->nbody3 = 0;
	hp->nbody4 = 0;
	hp->nbody5 = 0;
	hp->time = 0.;

	std::cout << "sat numb is: " << hp->sat_numb << std::endl;

	hp->tot_halo_mass = 0.;

	for (int j=0;j<hp->sat_numb;j++){

		hp->nbodySat.at(j) = 0;
		hp->sat_massBound.at(j) = 0.; //Describes the fracion of the initial satellite mass that is now still bound
		hp->sat_massBound_p4.at(j) = 0.;
		hp->sat_massBound_p5.at(j) = 0.;
	  	hp->sat_totMass.at(j) = 0.; //Describes the initial satellite mass	
	  	hp->sat_totMass_p4.at(j) = 0.;
	  	hp->sat_totMass_p5.at(j) = 0.;
	  	hp->sat_tidRadius.at(j) = 0.; //Tidal radius of satellite
  		hp->MW_integratedMass.at(j) = 0.; //MW mass within dens. center of MW and of satellite	   
  		hp->sat_tidalMass.at(j) = 0.;
		hp->sat_tidalMass_p4.at(j) = 0.;	
  		hp->sat_tidalMass_p5.at(j) = 0.;
  
        	hp->sat_outTidMass.at(j) = 0.;
		hp->sat_outTidMass_p4.at(j) = 0.;
		hp->sat_outTidMass_p5.at(j) = 0.;

	}	

		/*
		hp->sat_massDisk = 0.;
		hp->sat_massDisk_p4 = 0.;	
	  	hp->sat_massDisk_p5 = 0.;
	  	hp->sat_massBulge = 0.;
		hp->sat_massBulge_p4 = 0.;	
	  	hp->sat_massBulge_p5 = 0.;
  		hp->sat_massHalo = 0.;
		hp->sat_massHalo_p4 = 0.;	
  		hp->sat_massHalo_p5 = 0.;
		*/


		hp->halo_Mass_p4_tot = 0.;
		hp->halo_Mass_p5_tot = 0.;


	//array header reset		

	ap->xMin_rel = 0.;	
	ap->yMin_rel = 0.;	
	ap->zMin_rel = 0.;
	ap->xMax_rel = 0.;
	ap->yMax_rel = 0.;
	ap->zMax_rel = 0.;		
	for (int i = 0; i<6;i++){	

		ap->nX[i] = 0;	
		ap->nY[i] = 0;
		ap->nZ[i] = 0; 
		ap->res[i] = 0.;
	
	}	

	for (int i=0;i<40;i++){

		hp->disk_Vtan_p1[i] = 0.;
                hp->disk_Zed_p1[i] = 0.;
                hp->disk_Vzed_p1[i] = 0.;
                hp->disk_SigVzed_p1[i] = 0.;
                hp->disk_SigZed_p1[i] = 0.;
		hp->disk_Mass_p1[i] = 0.;
		hp->disk_Vtan_p4[i] = 0.;
 		hp->disk_Vzed_p4[i] = 0.;
 		hp->disk_Zed_p4[i]  = 0.;
		hp->disk_Mass_p4[i] = 0.;
		hp->disk_Vtan_p5[i] = 0.;
		hp->disk_Vzed_p5[i] = 0.;
		hp->disk_Zed_p5[i]  = 0.;
		hp->disk_Mass_p5[i] = 0.;

		hp->halo_ecc_p3[i] = 0.;
		hp->halo_ecc_p4[i] = 0.;
		hp->halo_ecc_p5[i] = 0.;	
		hp->smj_ax_p3[i] = 0.;
		hp->smj_ax_p4[i] = 0.;
		hp->smj_ax_p5[i] = 0.;
		hp->halo_Mass_p3[i] = 0.;
		for (int jj=0; jj<hp->sat_numb; jj++){

			hp->halo_Mass_p4[i+jj*40] = 0.;
			hp->halo_Mass_p5[i+jj*40] = 0.;

		}
		hp->halo_Vrad_p3[i] = 0.;
		hp->halo_Vtheta_p3[i] = 0.;
		hp->halo_Vphi_p3[i] = 0.;
		hp->halo_Vrad_p4[i] = 0.;
		hp->halo_Vtheta_p4[i] = 0.;
		hp->halo_Vphi_p4[i] = 0.;
		hp->halo_Vrad_p5[i] = 0.;
		hp->halo_Vtheta_p5[i] = 0.;
		hp->halo_Vphi_p5[i] = 0.;
		hp->halo_SigVrad_p3[i] = 0.;
		hp->halo_SigVtheta_p3[i] = 0.;
		hp->halo_SigVphi_p3[i] = 0.;
		hp->halo_SigVrad_p4[i] = 0.;
		hp->halo_SigVtheta_p4[i] = 0.;
		hp->halo_SigVphi_p4[i] = 0.;
		hp->halo_SigVrad_p5[i] = 0.;
		hp->halo_SigVtheta_p5[i] = 0.;
		hp->halo_SigVphi_p5[i] = 0.;		


		/*		
		hp->halo_Vrad_p3[i] = 0.;
		hp->halo_Vtheta_p3[i] = 0.;
		hp->halo_Vphi_p3[i] = 0.;
		hp->halo_Mass_p3[i] = 0.;
		hp->halo_Sigrad_p3[i] = 0.;
		hp->halo_Sigtheta_p3[i] = 0.;
		hp->halo_Sigphi_p3[i] = 0.;								
		hp->bulge_Vrad_p4[i] = 0.;
		hp->bulge_Vtheta_p4[i] = 0.;
		hp->bulge_Vphi_p4[i] = 0.;
		hp->bulge_Sigrad_p4[i] = 0.;
		hp->bulge_Sigtheta_p4[i] = 0.;
		hp->bulge_Sigphi_p4[i] = 0.;
		hp->disk_Vtan_p4[i] = 0.;
		*/		

			//hp->disk_Vtan_p5[i] = 0.;

			
			/*			
			hp->bulge_Vrad_p5[i] = 0.;
			hp->bulge_Vtheta_p5[i] = 0.;
			hp->bulge_Vphi_p5[i] = 0.;
			
			hp->bulge_Sigrad_p5[i] = 0.;
			hp->bulge_Sigtheta_p5[i] = 0.;
			hp->bulge_Sigphi_p5[i] = 0.;
			hp->halo_Vrad_p5[i] = 0.;
			hp->halo_Vtheta_p5[i] = 0.;
			hp->halo_Vphi_p5[i] = 0.;
			hp->halo_Mass_p5[i] = 0.;
			hp->halo_Sigrad_p5[i] = 0.;
			hp->halo_Sigtheta_p5[i] = 0.;
			hp->halo_Sigphi_p5[i] = 0.;
			*/

	
		/*
		for(int j=0;j<4;j++){	

			hp->disk_vertical_zone_p1[j][i] = 0.;
			hp->disk_vertical_zone_p4[j][i] = 0.;
			if(hp->wantBaryons == 1){	
				hp->disk_vertical_zone_p5[j][i] = 0.;		
			}

		}

		*/

		

	}

	for (int i=0;i<50;i++){

		hp->halo_Mass_hiRes_p3[i] = 0.;

		for (int jj=0; jj<hp->sat_numb; jj++){

			hp->halo_Mass_hiRes_p4[i+jj*50] = 0.;
			hp->halo_Mass_hiRes_p5[i+jj*50] = 0.;

		}	
	
	}

	for (int i=0;i<250;i++){

		hp->halo_Mass_hiRes2_p3[i] = 0.;

		for (int jj=0; jj<hp->sat_numb; jj++){

			hp->halo_Mass_hiRes2_p4[i+jj*250] = 0.;
			hp->halo_Mass_hiRes2_p5[i+jj*250] = 0.;

		}	
	
	}

	for (int i=0;i<25;i++){

		hp->halo_Mass_hiRes0_p3[i] = 0.;

		for (int jj=0; jj<hp->sat_numb; jj++){

			hp->halo_Mass_hiRes0_p4[i+jj*25] = 0.;
			hp->halo_Mass_hiRes0_p5[i+jj*25] = 0.;

		}	
	
	}


	for(int i=0;i<10;i++){

		hp->bulge_Vrad_p2[i] = 0.;
		hp->bulge_Vtheta_p2[i] = 0.;
		hp->bulge_Vphi_p2[i] = 0.;
		hp->bulge_Sigrad_p2[i] = 0.;
		hp->bulge_Sigtheta_p2[i] = 0.;
		hp->bulge_Sigphi_p2[i] = 0.;
		hp->bulge_Sigma_p2[i] = 0.;
		hp->bulge_Mass_p2[i] = 0.;	
		hp->bulge_Mass_p4[i] = 0.;	
		hp->bulge_Mass_p5[i] = 0.;	

	}	

	/*

	for(int i=0;i<8;i++){

		hp->disk_vertical_total_p1[i] = 0.;
		hp->disk_vertical_total_p4[i] = 0.;
		if(hp->wantBaryons == 1){
			hp->disk_vertical_total_p5[i] = 0.;
		}	
	}

	for(int i=0;i<20;i++){

		hp->bulge_vels_p2[i] = 0.;
		hp->halo_vels_p3[i] = 0.;
		hp->bulge_vels_p4[i] = 0.;
		hp->halo_vels_p4[i] = 0.;
		if(hp->wantBaryons == 1){
			hp->halo_vels_p5[i] = 0.;	
			hp->bulge_vels_p5[i] = 0.;
		}
	}

	*/	

}



bool acompare(struct h5_particle ptL, struct h5_particle ptR){

	return (ptL.idf < ptR.idf);

}

void read_particle_numbers(struct header_h5 *hp){

	std::cout <<"Debug read number particles" << std::endl;

	std::string myNumberName = hp->AquariusClass;
	std::ifstream myNumberTable(("which_satellite_numbers_Aq_"+myNumberName+".txt").c_str());
	std::string myNumberLine;

	int iLine = 0;

	while (std::getline(myNumberTable,myNumberLine)){

		std::istringstream myNumberStream(myNumberLine);
		int trash;

		myNumberStream >> hp->myNumbDark[iLine];
		myNumberStream >> hp->myNumbStar[iLine];
		myNumberStream >> trash;
	
		++iLine;

	}

	std::cout << "	Qui ci siamo arivati ))) " << std::endl;


}

