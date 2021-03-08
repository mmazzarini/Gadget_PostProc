#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <hdf5.h>
#include "H5Cpp.h"
#include "library_MMazzarini_structs.h"
#include "library_MMazzarini_utilities.h"
#include "library_MMazzarini_density.h"
#include "library_MMazzarini_arrays.h"
#include "library_MMazzarini_readwriteHDF5.h"
#include "library_MMazzarini_mass_distribution.h"
#include "library_MMazzarini_density.h"

using namespace H5;

int main(int argc, char* argv[]){

	std::string myFileRoot;
	int first;
	int last;
	int step;

	struct header_h5 hp;
	std::vector<struct h5_particle> pp0;
	std::vector<struct h5_particle> pp1;
	std::vector<struct h5_particle> pp2;
	std::vector<struct h5_particle> pp3;
	std::vector<struct h5_particle> pp4;
	std::vector<struct h5_particle> pp5;

	mass_center dcms;
	std::vector<struct mass_center> dcms_Sat;
 	mass_center dcms_Baryon;
	mass_center dcms_Dark;
	
        struct array_hd ap; //Vector containing array-header objects...they give info about visual maps for particles of the 6 types

	std::cout << "Starting interactive kernel...\n\n" << std::endl;

	//Launch kernel to ask for input parameters
	kernel_general(&hp, &ap);

	//Set structures to zero

	/*
	*
	* Start kernel and the code main operations...
	*
	*/
		
	hp.nbodySat.resize(hp.sat_numb);  //tell you how many particles per type
	hp.softenings.resize(hp.sat_numb+6);	
	hp.tot_halo_mass = 0.;
	hp.sat_massBound.resize(hp.sat_numb);	
	hp.sat_massBound_p4.resize(hp.sat_numb);
	hp.sat_massBound_p5.resize(hp.sat_numb);
	hp.sat_totMass.resize(hp.sat_numb);
	hp.sat_totMass_p4.resize(hp.sat_numb);
	hp.sat_totMass_p5.resize(hp.sat_numb);
	hp.sat_tidRadius.resize(hp.sat_numb);
	hp.MW_integratedMass.resize(hp.sat_numb);
	hp.sat_tidalMass.resize(hp.sat_numb);
	hp.sat_tidalMass_p4.resize(hp.sat_numb);
	hp.sat_tidalMass_p5.resize(hp.sat_numb);
	hp.sat_ecc.resize(hp.sat_numb);
	hp.sat_ax.resize(hp.sat_numb);
	hp.sat_outTidMass.resize(hp.sat_numb);
	hp.sat_outTidMass_p4.resize(hp.sat_numb);
	hp.sat_outTidMass_p5.resize(hp.sat_numb);
	//Resize dcms
	dcms_Sat.resize(hp.sat_numb);
	hp.halo_Mass_p4.resize(hp.sat_numb*40);
	hp.halo_Mass_p5.resize(hp.sat_numb*40);
	hp.halo_Mass_hiRes0_p4.resize(hp.sat_numb*25);
	hp.halo_Mass_hiRes0_p5.resize(hp.sat_numb*25);
	hp.halo_Mass_hiRes_p4.resize(hp.sat_numb*50);
	hp.halo_Mass_hiRes_p5.resize(hp.sat_numb*50);
	hp.halo_Mass_hiRes2_p4.resize(hp.sat_numb*250);
	hp.halo_Mass_hiRes2_p5.resize(hp.sat_numb*250);
	hp.myNumbDark.resize(hp.sat_numb);
	hp.myNumbStar.resize(hp.sat_numb); 

	for(int isf=0;isf<hp.sat_numb+6;isf++){
		hp.softenings[isf] = ap.res[isf];
	}

	std::cout << "Setting all structures to 0" << std::endl;

	reset_structures(&dcms, &dcms_Sat, &pp0,&pp1, &pp2, &pp3, &pp4, &pp5, &hp, &ap);

	std::cout << "softenings are: ";

	for(int isf=0; isf<hp.sat_numb+6;isf++){

		std::cout << hp.softenings[isf] << " " ;

	}

	
	std::cout << std::endl;

	std::cout << "Reading particle numbers for satellites" << std::endl;

	
	read_particle_numbers(&hp);


	int iFile;
	for(iFile=hp.first;iFile<hp.last;iFile+=hp.step){ 

		//Convert number in string

		std::string fileNumber = myString(iFile);

	        if(iFile<10){
	
			hp.myHDF5File = hp.myHDF5Root+"_snapshot_00"+fileNumber+".hdf5"; 
			hp.myHDF5Post = hp.myHDF5Root+"_snapshot_00"+fileNumber+"_postProc.hdf5"; 

		}else if(iFile>=10 and iFile <100){

			hp.myHDF5File = hp.myHDF5Root+"_snapshot_0"+fileNumber+".hdf5"; 
			hp.myHDF5Post = hp.myHDF5Root+"_snapshot_0"+fileNumber+"_postProc.hdf5";

		}else{					

			hp.myHDF5File = hp.myHDF5Root+"_snapshot_"+fileNumber+".hdf5"; 
			hp.myHDF5Post = hp.myHDF5Root+"_snapshot_"+fileNumber+"_postProc.hdf5";

		}

		/*
	
		}else if(iFile>=10000){					

			hp.myHDF5Post = hp.myHDF5Root+"_snapshot_"+fileNumber+"_postProc.hdf5";

		}	

		*/
		
		//Read hdf5 file and store all data about particles
		std::cout << "Reading HDF5 file..." << hp.myHDF5File << std::endl;
	
		int myTime = iFile;

		read_HDF5_file(&pp0, &pp1, &pp2, &pp3, &pp4, &pp5, &hp);

		for(int ip=0;ip<hp.nbody0;ip++){
			pp0[ip].eps = hp.softenings[0];
		}
		
		for(int ip=0;ip<hp.nbody1;ip++){
			pp1[ip].eps = hp.softenings[1];
		}
		
		for(int ip=0;ip<hp.nbody2;ip++){
			pp2[ip].eps = hp.softenings[2];
		}
		
		for(int ip=0;ip<hp.nbody3;ip++){
			pp3[ip].eps = hp.softenings[3];
		}
		
		//Particle 4-5

		for (int ip=0; ip<hp.nbody4;ip++){

			pp4[ip].eps = hp.softenings[4];
						
		}

		for(int ip=0; ip<hp.nbody5; ip++){

			pp5[ip].eps = hp.softenings[5];

		}

		std::cout << "debug softening pp4: " << pp4[0].eps << std::endl;
		std::cout << "debug softening pp5: " << pp5[0].eps << std::endl;
	
		std::cout << "Sorting satellite dark particles" <<std::endl;

		std::sort(pp4.begin(),pp4.end(),acompare);

                std::cout << "Sorting satellite star particles"	<<std::endl;

		std::sort(pp5.begin(),pp5.end(),acompare);

		std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n My mass for particle random is " << pp4[10].mass << std::endl;
		std::cout << "My softening for particle random is " << pp4[10].eps << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << std::endl;

		//Calculate Milky Way density center -> This version: pure Nbody, 0-type particles have no content
	
		density_center_calc_MW(&dcms, &pp1, &pp2, &pp3, &hp);   //In box-intertial-frame coordinates
		//Calculate Satellite density center 
		density_center_calc_Sat(&dcms_Sat, &pp4, &pp5, &hp);    //In box-intertial-frame coordinates

		//read_densities(&dcms,&dcms_Sat,&hp,&myTime);

		//Calculate Satellite Baryons density center 
		
		/*	
		if(hp.wantBaryons == 1){

		density_center_calc_SatelliteBaryon(&dcms_Baryon, &pp4, &hp);

			//Calculate Satellite Dark Matter density center 
			density_center_calc_SatelliteDark(&dcms_Dark, &pp5, &hp);
		}	
		
		*/

		std::cout << "	Calculating bound mass, integrated MW mass, tidal radius and tidal mass fraction" << std::endl;

		//Calculate fraction of bound mass
		bound_mass_calculate(&dcms_Sat, &pp4, &pp5, &hp);

		std::cout << "\n\n\n" << std::endl;

		for(int i=0;i<hp.sat_numb;i++){

			std::cout << "My mass for satellite i is " << hp.sat_totMass[i] << std::endl;

		}	

		std::cout << "\n\n\n" << std::endl;
	
		//Calculate MW mass within satellite distance to MW density center
		integrate_MW_mass(&dcms, &dcms_Sat, &pp1, &pp2, &pp3, &hp);	

		//Write header in post_proc file
		tidal_radius_calculate(&dcms, &dcms_Sat, &hp);

		//Calculate tidal mass fraction	
		calculate_TidalMass_Fraction(&dcms_Sat, &pp4, &pp5, &hp,fileNumber);

		std::cout << "Relocating coordinates according to MW density center position..." << std::endl;

		//Relocate positions of particles with respect to density center
	
		density_center_relocate(&dcms,&pp1,&(hp.nbody1));
		density_center_relocate(&dcms,&pp2,&(hp.nbody2));				
		density_center_relocate(&dcms,&pp3,&(hp.nbody3));
		density_center_relocate(&dcms,&pp4,&(hp.nbody4));
		if(hp.wantBaryons == 1){
			density_center_relocate(&dcms,&pp5,&(hp.nbody5));
		}

		//Calculate angular momentum of disk
		angular_momentum_MW_calculate(&dcms,&pp1,&pp2,&pp3,&hp);

		//Correct for angular momentum inclination
		angular_momentum_corrections(&dcms,&pp1,&pp2,&pp3,&pp4,&pp5,&hp);	

		//Do flagging on energy-content-based membership of satellite particles
		std::cout << "Flagging satellite particles for distribution analysis..." << std::endl;

		//Calculate distribution properties of baryonic (4) and dark (4) matter

		mass_distribution_flag(&pp4,&pp5,&hp,&iFile);

		mass_distribution_disk(&pp1, &pp4,&pp5,&hp);
		mass_distribution_bulge(&pp2, &pp4,&pp5, &hp);
		
		//Check distribution of debris into the halo

		
		
		mass_distribution_halo(&pp3, &pp4, &pp5, &dcms, &dcms_Sat, &hp);
	
		std::cout << "\nlet's print the vectors\n" << std::endl; 

		for (int s=0;s<40;s++){

			std::cout << hp.halo_Mass_p4[s] << " ";

		}

		std::cout << std::endl;
	
		//Let's store now all header data in HDF5 file... 

		std::cout << "Calculating tensor of inertia..." << std::endl;

		calc_tensor_inertia(&pp3, &pp4, &pp5, &hp, &dcms, &dcms_Sat,fileNumber);

		std::cout << "Storing data in HDF5 postproc. file " << hp.myHDF5Post << std::endl;	
	 	write_Hdf5_postproc(&pp1,&pp2,&pp3,&pp4,&pp5,&hp,&dcms,&dcms_Sat,&dcms_Baryon,&dcms_Dark);	
		std::cout << "Done with current file..." << std::endl;	

		//if (iFile == 40){

			calc_En_Ang_25kpc(&pp3, &pp4, &pp5, &hp, hp.myHDF5Root+"_snapshot_0"+fileNumber+"_Halo.txt", hp.myHDF5Root+"_snapshot_0"+fileNumber+"_Dark.txt", hp.myHDF5Root+"_snapshot_0"+fileNumber+"_Star.txt");
			calc_Toomre_Space(&pp3, &pp4, &pp5, &hp, hp.myHDF5Root+"_snapshot_0"+fileNumber+"_Halo.txt", hp.myHDF5Root+"_snapshot_0"+fileNumber+"_Dark.txt", hp.myHDF5Root+"_snapshot_0"+fileNumber+"_Star.txt");

			calc_Energy_Lang_DcmsSat(&dcms, &dcms_Sat, &pp4, &pp5, &hp);

		//}

		//Reset structures to zero

		std::cout << "the MW density center is" << dcms.pos[0] << " " << dcms.pos[1] << " " << dcms.pos[2] << " " << dcms.vel[0] << " " << dcms.vel[1] << " " << dcms.vel[2] << std::endl;

		reset_structures(&dcms, &dcms_Sat, &pp0,&pp1, &pp2, &pp3, &pp4, &pp5, &hp, &ap);
		
	}

} //End main
