#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <hdf5.h>
#include "H5Cpp.h"
#include "library_MMazzarini_structs.h"
#include "library_MMazzarini_readwriteHDF5.h"

using namespace H5;

void read_HDF5_file(std::vector<struct h5_particle> *pt0, std::vector<struct h5_particle> *pt1, std::vector<struct h5_particle> *pt2, std::vector<struct h5_particle> *pt3, std::vector<struct h5_particle> *pt4, std::vector<struct h5_particle> *pt5,struct header_h5 *hp){

	//Set all
	int Dim2 = 3;
	herr_t ret;
	
	//Hyperslab elements for dataset

        hsize_t start_dset[2];
	hsize_t count_dset[2];
	hsize_t stride_dset[2];
	hsize_t block_dset[2];
	hsize_t block_vert_dset[2];
	hsize_t block_short_dset[2];
	hsize_t block_short_vert_dset[2];
	start_dset[1]  = 0;
 	stride_dset[0] = 1; stride_dset[1] = 1;
	count_dset[0]  = 1; count_dset[1]  = 1;    
	block_dset[0]  = 10; block_dset[1]  = Dim2;
	block_vert_dset[0] = 10; block_vert_dset[1]  = 1;
	block_short_dset[1] = Dim2;
	block_short_vert_dset[1] = 1;

	//Hyperslab elements for buffer; 

        hsize_t start_buff[2];
	hsize_t count_buff[2];
	hsize_t stride_buff[2];
	hsize_t block_buff[2];
	hsize_t block_vert_buff[2];
	hsize_t block_short_buff[2];
	hsize_t block_short_vert_buff[2];
	start_buff[0]  = 0; start_buff[1]  = 0;
 	stride_buff[0] = 1; stride_buff[1] = 1;
	count_buff[0]  = 1; count_buff[1]  = 1;    
	block_buff[0]  = 10; block_buff[1]  = Dim2;
	block_vert_buff[0] = 10; block_vert_buff[1]  = 1;
	block_short_buff[1] = Dim2;
	block_short_vert_buff[1] = 1;
	
	std::cout << "Reading file: " << hp->myHDF5File << "\n" << std::endl;
 
	//Open file

	std::ifstream fileCheck((hp->myHDF5File).c_str());
	if (!fileCheck){
		std::cout << "[ERROR]: File " << hp->myHDF5File << " doesn't exists. Aborting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	hid_t myFileRead = H5Fopen((hp->myHDF5File).c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
	

	//Open header

        hid_t headerRead = H5Gopen(myFileRead,"/Header",H5P_DEFAULT);		


	//Set attribute spaces

	hsize_t attr_dims_vect = 6;
	hsize_t attr_dims_scal = 1;
	
	hid_t vect_Aspace = H5Screate_simple(1,&attr_dims_vect,NULL);
	hid_t scal_Aspace = H5Screate_simple(1,&attr_dims_scal,NULL);


        //Set arrays for Header Attributes
	

	hid_t att_time = H5Aopen_name(headerRead,"Time");
	H5Aread(att_time,H5T_NATIVE_FLOAT,&hp->time);
	H5Aclose(att_time);	

	hid_t att_sfr = H5Aopen_name(headerRead,"Flag_Sfr");
	H5Aread(att_sfr,H5T_NATIVE_INT,&hp->sfr);
	H5Aclose(att_sfr);
	
	hid_t att_cooling = H5Aopen_name(headerRead,"Flag_Cooling");
	H5Aread(att_cooling,H5T_NATIVE_INT,&hp->cool);
	H5Aclose(att_cooling);

	hid_t att_feed = H5Aopen_name(headerRead,"Flag_Feedback");
	H5Aread(att_feed,H5T_NATIVE_INT,&hp->feed);
	H5Aclose(att_feed);

	hid_t att_redsh = H5Aopen_name(headerRead,"Redshift");
	H5Aread(att_redsh,H5T_NATIVE_FLOAT,&hp->redshift);
	H5Aclose(att_redsh);

	hid_t att_box = H5Aopen_name(headerRead,"BoxSize");
	H5Aread(att_box,H5T_NATIVE_FLOAT,&hp->boxSize);
	H5Aclose(att_box);

	hid_t att_omega0 = H5Aopen_name(headerRead,"Omega0");
	H5Aread(att_omega0,H5T_NATIVE_FLOAT,&hp->omega0);
	H5Aclose(att_omega0);

	hid_t att_omegaL = H5Aopen_name(headerRead,"OmegaLambda");
	H5Aread(att_omegaL,H5T_NATIVE_FLOAT,&hp->omegaL);
	H5Aclose(att_omegaL);

	hid_t att_Hubble = H5Aopen_name(headerRead,"HubbleParam");
	H5Aread(att_Hubble,H5T_NATIVE_FLOAT,&hp->Hubble);
	H5Aclose(att_Hubble);
	
	hid_t att_doublePrec = H5Aopen_name(headerRead,"Flag_DoublePrecision");
	H5Aread(att_doublePrec,H5T_NATIVE_INT,&hp->dPrec);
	H5Aclose(att_doublePrec);

	hid_t att_nFile = H5Aopen_name(headerRead,"NumFilesPerSnapshot");
	H5Aread(att_nFile,H5T_NATIVE_INT,&hp->nFile);
	H5Aclose(att_nFile);

	hid_t att_mass = H5Aopen_name(headerRead,"MassTable");
	H5Aread(att_mass,H5T_NATIVE_FLOAT,hp->mass_table);
	H5Aclose(att_mass);	

	hid_t att_nPart_this = H5Aopen_name(headerRead,"NumPart_ThisFile");
	H5Aread(att_nPart_this,H5T_NATIVE_INT,hp->nPart_here_table);
	H5Aclose(att_nPart_this);

	hid_t att_nPart_tot = H5Aopen_name(headerRead,"NumPart_Total");
	H5Aread(att_nPart_tot,H5T_NATIVE_INT,hp->nPart_tot_table);
	H5Aclose(att_nPart_tot);

	hid_t att_nPart_high = H5Aopen_name(headerRead,"NumPart_Total_HighWord");
	H5Aread(att_nPart_high,H5T_NATIVE_INT,hp->nPart_high_table);
	H5Aclose(att_nPart_high);
	
	H5Gclose(headerRead);

	//Update header

	int i;

	hp->nbodies = 0;
	hp->nbody0 = 0;
	hp->nbody1 = 0;
	hp->nbody2 = 0;
	hp->nbody3 = 0;
	hp->nbody4 = 0;
	hp->nbody5 = 0;

        for(i=0;i<6;i++){

	   hp->nbodies += hp->nPart_here_table[i];
	
	}

	hp->nbody0 = hp->nPart_here_table[0];
	hp->nbody1 = hp->nPart_here_table[1];
	hp->nbody2 = hp->nPart_here_table[2];
	hp->nbody3 = hp->nPart_here_table[3];
	hp->nbody4 = hp->nPart_here_table[4];
	hp->nbody5 = hp->nPart_here_table[5];
	
	if(hp->cool == 1){
           hp ->nsph = hp->nPart_here_table[0];	
	}
	else{
	   hp->nsph = 0;
	}

	std::cout << "Reading header parameters: " << std::endl;
	std::cout << " " << std::endl;
	std::cout << "Reading time        : " << hp->time << std::endl;
	std::cout << "Reading sfr         : " << hp->sfr  << std::endl;
	std::cout << "Reading cool        : " << hp->cool << std::endl;
	std::cout << "Reading feed        : " << hp->feed << std::endl;
	std::cout << "Reading redshift    : " << hp->redshift << std::endl;
	std::cout << "Reading boxSize     : " << hp->boxSize << std::endl;
	std::cout << "Reading omega0      : " << hp->omega0 << std::endl;
	std::cout << "Reading omegaL      : " << hp->omegaL << std::endl;
	std::cout << "Reading Hubble      : " << hp->Hubble << std::endl;
	std::cout << "Reading dPrec       : " << hp->dPrec << std::endl;
	std::cout << "Reading nFile       : " << hp->nFile << std::endl;
	std::cout << "Reading nbody0      : " << hp->nbody0 << std::endl;
	std::cout << "Reading nbody1      : " << hp->nbody1 << std::endl;
	std::cout << "Reading nbody2      : " << hp->nbody2 << std::endl;
	std::cout << "Reading nbody3      : " << hp->nbody3 << std::endl;
	std::cout << "Reading nbody4      : " << hp->nbody4 << std::endl;
	std::cout << "Reading nbody5      : " << hp->nbody5 << std::endl;

	std::cout << "Reading nsph        : " << hp->nsph << std::endl;
	std::cout << "Reading nbodies     : " << hp->nbodies << std::endl;

	std::cout << "Reading mass_table  : [ " << hp->mass_table[0] << ", " << hp->mass_table[1] << ", " << hp->mass_table[2] << ", " << hp->mass_table[3] << ", " << hp->mass_table[4] << ", " << hp->mass_table[5] << " ]" << std::endl;
	std::cout << "Reading nParts_here : [ " << hp->nPart_here_table[0] << ", " << hp->nPart_here_table[1] << ", " << hp->nPart_here_table[2] << ", " << hp->nPart_here_table[3] << ", " << hp->nPart_here_table[4] << ", " << hp->nPart_here_table[5] << " ]" << std::endl;
	std::cout << "Reading nPart_tot   : [ " << hp->nPart_tot_table[0] << ", " << hp->nPart_tot_table[1] << ", " << hp->nPart_tot_table[2] << ", " << hp->nPart_tot_table[3] << ", " << hp->nPart_tot_table[4] << ", " << hp->nPart_tot_table[5] << " ]" << std::endl;
	std::cout << "Reading Part_high   : [ " << hp->nPart_high_table[0] << ", " << hp->nPart_high_table[1] << ", " << hp->nPart_high_table[2] <<  ", " << hp->nPart_high_table[3] << ", " << hp->nPart_high_table[4] << ", " << hp->nPart_high_table[5] << " ]\n" << std::endl;


	pt0->resize(hp->nbody0);
	pt1->resize(hp->nbody1);
	pt2->resize(hp->nbody2);
	pt3->resize(hp->nbody3);
	pt4->resize(hp->nbody4);
	pt5->resize(hp->nbody5);

	//Now read particle data

	int grindex;

	for(grindex=0;grindex<6;grindex++){

           std::string indstr;		
           std::ostringstream iconv;
           iconv << grindex;
           indstr = iconv.str(); 
 
	   std::cout << "Reading data for Particle Type "+indstr+"..." << std::endl;

	  	 
	   if(hp->nPart_here_table[grindex] != 0){
	 
              int myNpart = hp->nPart_here_table[grindex];			

	      std::cout << "    Found number of particles " << myNpart << std::endl;	


	      //Create dataSpaces for Groups

	      hsize_t Part_Dims_1D;
	      hsize_t Part_Dims_2D[2];
	      Part_Dims_1D = myNpart;
	      Part_Dims_2D[0] = myNpart;
	      Part_Dims_2D[1] = Dim2;
		
	      hid_t part_Dspace_1D = H5Screate_simple(1,&Part_Dims_1D,NULL);
	      hid_t part_Dspace_2D = H5Screate_simple(2,Part_Dims_2D,NULL);
	      hid_t part_Buffspace_1D = H5Screate_simple(1,&Part_Dims_1D,NULL);
	      hid_t part_Buffspace_2D = H5Screate_simple(2,Part_Dims_2D,NULL);

              int fracPart = ((int)myNpart/10);
	      int restPart = ((int)myNpart%10);	

              float myCoord[10][3];
              float myVels[10][3];
              float myAcc[10][3];
	      float myMass[10];
	      float myPot[10];
	      int myIDs[10];

		

              float myCoord_short[restPart][3];
	      float myVels_short[restPart][3];
              float myAcc_short[restPart][3];
	      float myMass_short[restPart];
	      float myPot_short[restPart];
	      int myIDs_short[restPart]; 	

	      //Open particle group	

	      hid_t myPartRead = H5Gopen(myFileRead,("/PartType"+indstr).c_str(),H5P_DEFAULT);	

	      //Open datasets			
					      		
	      herr_t stat_coord;	
	      herr_t stat_vels;	
	      herr_t stat_acc;	
	      herr_t stat_mass;	
	      herr_t stat_pot;	
	      herr_t stat_IDs;	

	      H5L_info_t trash_info;	

	      stat_coord = H5Lget_info(myPartRead, "Coordinates", &trash_info,H5P_DEFAULT);
	      if (stat_coord == 0){
		 std::cout << "  Object Coordinates found..." 	<< std::endl;
	      }else{
		 std::cout << "\n  #WNG: Object Coordinates NOT found...\n" << std::endl;
		 exit(EXIT_FAILURE); 	
	      }

	      stat_vels = H5Lget_info(myPartRead, "Velocities", &trash_info,H5P_DEFAULT);
	      if (stat_vels == 0){
		 std::cout << "  Object Velocities found..." 	<< std::endl;
	      }else{
		 std::cout << "\n  #WNG: Object Velocities NOT found...\n" << std::endl;
		 exit(EXIT_FAILURE);
	      }

	      stat_acc = H5Lget_info(myPartRead, "Acceleration", &trash_info,H5P_DEFAULT);
	      if (stat_acc == 0){
		 std::cout << "  Object Acceleration found..." 	<< std::endl;
	      }else{
		 std::cout << "\n  #WNG: Object Acceleration NOT found...\n" << std::endl;
		 exit(EXIT_FAILURE);
	      }

	      stat_mass = H5Lget_info(myPartRead, "Masses", &trash_info,H5P_DEFAULT);
	      if (stat_mass == 0){
		 std::cout << "  Object Masses found..." 	<< std::endl;
	      }else{
		 std::cout << "\n  #WNG: Object Masses NOT found...\n" << std::endl;
		 exit(EXIT_FAILURE);
	      }

	      stat_pot = H5Lget_info(myPartRead, "Potential", &trash_info,H5P_DEFAULT);
	      if (stat_pot == 0){
		 std::cout << "  Object Potential found..." 	<< std::endl;
	      }else{
		 std::cout << "\n  #WNG: Object Potential NOT found...\n" << std::endl;
		 exit(EXIT_FAILURE);
	      }

	      stat_IDs = H5Lget_info(myPartRead, "ParticleIDs", &trash_info,H5P_DEFAULT);
	      if (stat_IDs == 0){
		 std::cout << "  Object ParticleIDs found..." 	<< std::endl;
	      }else{
		 std::cout << "\n  #WNG: Object ParticleIDs NOT found...\n" << std::endl;
		 exit(EXIT_FAILURE);
	      }
		
	//      if(stat_coord == 0)	 

	      hid_t pCoord = H5Dopen(myPartRead,"Coordinates",H5P_DEFAULT);
	      hid_t pVels = H5Dopen(myPartRead,"Velocities",H5P_DEFAULT);
	      hid_t pAcc = H5Dopen(myPartRead,"Acceleration",H5P_DEFAULT);
	      hid_t pMass = H5Dopen(myPartRead,"Masses",H5P_DEFAULT);
              hid_t pPot = H5Dopen(myPartRead,"Potential",H5P_DEFAULT);
	      hid_t pIDs = H5Dopen(myPartRead,"ParticleIDs",H5P_DEFAULT);

	      int myid;

	      if(fracPart != 0){
	
	         for (myid = 0;myid<fracPart;myid++){	

         	   start_dset[0]  = myid*10;
 

		   if (stat_coord == 0){	
                      ret = H5Sselect_hyperslab(part_Dspace_2D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_dset);
                      ret = H5Sselect_hyperslab(part_Buffspace_2D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_buff);	
	              ret = H5Dread(pCoord, H5T_NATIVE_FLOAT, part_Buffspace_2D, part_Dspace_2D, H5P_DEFAULT, myCoord);
		      if (ret<0){
			  std::cout << "Cannot read coordinates of type "+indstr+" in fracPart" << myid << " out of " << fracPart << std::endl; 	
		          exit(EXIT_FAILURE);
		      }				   	
		   }

		   if (stat_vels == 0){
                      ret = H5Sselect_hyperslab(part_Dspace_2D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_dset);	
	              ret = H5Sselect_hyperslab(part_Buffspace_2D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_buff);	
	              ret = H5Dread(pVels, H5T_NATIVE_FLOAT, part_Buffspace_2D, part_Dspace_2D, H5P_DEFAULT, myVels);
		      if (ret<0){
			  std::cout << "Cannot read velocities of type "+indstr+" in fracPart" << myid << " out of " << fracPart << std::endl; 		
		          exit(EXIT_FAILURE);
		      }				   	
		   }		           

		   if (stat_acc == 0){
                      ret = H5Sselect_hyperslab(part_Dspace_2D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_dset);	
	              ret = H5Sselect_hyperslab(part_Buffspace_2D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_buff);	
	              ret = H5Dread(pAcc, H5T_NATIVE_FLOAT, part_Buffspace_2D, part_Dspace_2D, H5P_DEFAULT, myAcc);
		      if (ret<0){
			  std::cout << "Cannot read accelerations of type "+indstr+" in fracPart" << myid << " out of " << fracPart << std::endl; 		       
		           exit(EXIT_FAILURE);
		      }				   	
		   }	
	
		   if (stat_mass == 0){	     		
 	              ret = H5Sselect_hyperslab(part_Dspace_1D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_vert_dset);	
	              ret = H5Sselect_hyperslab(part_Buffspace_1D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_vert_buff);	
	              ret = H5Dread(pMass, H5T_NATIVE_FLOAT, part_Buffspace_1D, part_Dspace_1D, H5P_DEFAULT, myMass);	
		      if (ret<0){
			  std::cout << "Cannot read masses of type "+indstr+" in fracPart" << myid << " out of " << fracPart << std::endl; 	
	                  exit(EXIT_FAILURE);
		        }				   		      
		   } 		     
			
		   if (stat_pot == 0){
 	              ret = H5Sselect_hyperslab(part_Dspace_1D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_vert_dset);	
	              ret = H5Sselect_hyperslab(part_Buffspace_1D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_vert_buff);	
	              ret = H5Dread(pPot, H5T_NATIVE_FLOAT, part_Buffspace_1D, part_Dspace_1D, H5P_DEFAULT, myPot);		
		      if (ret<0){
			  std::cout << "Cannot read potentials of type "+indstr+" in fracPart" << myid << " out of " << fracPart << std::endl; 		       
		           exit(EXIT_FAILURE);
		      }				   	      
		   }	

		   if (stat_IDs == 0){
 	              ret = H5Sselect_hyperslab(part_Dspace_1D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_vert_dset);	
	              ret = H5Sselect_hyperslab(part_Buffspace_1D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_vert_buff);	
	              ret = H5Dread(pIDs, H5T_NATIVE_INT, part_Buffspace_1D, part_Dspace_1D, H5P_DEFAULT, myIDs);
		      if (ret<0){
			  std::cout << "Cannot read IDs of type "+indstr+" in fracPart" << myid << " out of " << fracPart << std::endl; 
			   exit(EXIT_FAILURE);
		      }				   			      
		   }
		   
		   int start = myid*10;
                   int end   = myid*10+10;	 		  	

		   myFunc_Assign_H5parts(pt0,pt1,pt2,pt3,pt4,pt5,&grindex,&start,&end,myCoord,myVels,myAcc,myMass,myPot,myIDs);
		       	

         	} //Close loop on partFrac    	

             } //Close if on partFrac    		


             


	     if(restPart != 0){

                start_dset[0]  = fracPart*10;
                block_short_dset[0] = restPart;
	        block_short_vert_dset[0] = restPart;
                block_short_buff[0] = restPart;
	        block_short_vert_buff[0] = restPart; 

		
		if (stat_coord == 0){    
		   ret = H5Sselect_hyperslab(part_Dspace_2D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_short_dset);	
                   ret = H5Sselect_hyperslab(part_Buffspace_2D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_short_buff);	
	           ret = H5Dread(pCoord, H5T_NATIVE_FLOAT, part_Buffspace_2D, part_Dspace_2D, H5P_DEFAULT, myCoord_short);
		   if (ret<0){
			  std::cout << "Cannot read coordinates of type "+indstr+" in restPart" << std::endl;
 		          exit(EXIT_FAILURE);
		      }	
		}

		if (stat_vels == 0){
                   ret = H5Sselect_hyperslab(part_Dspace_2D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_short_dset);	
	           ret = H5Sselect_hyperslab(part_Buffspace_2D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_short_buff);	
	           ret = H5Dread(pVels, H5T_NATIVE_FLOAT, part_Buffspace_2D, part_Dspace_2D, H5P_DEFAULT, myVels_short);
		   if (ret<0){
			  std::cout << "Cannot read velocities of type "+indstr+" in restPart" << std::endl;
 		          exit(EXIT_FAILURE);
		      }
		}

		if (stat_acc == 0){
                   ret = H5Sselect_hyperslab(part_Dspace_2D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_short_dset);	
	           ret = H5Sselect_hyperslab(part_Buffspace_2D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_short_buff);	
	           ret = H5Dread(pAcc, H5T_NATIVE_FLOAT, part_Buffspace_2D, part_Dspace_2D, H5P_DEFAULT, myAcc_short);	
		   if (ret<0){
			  std::cout << "Cannot read accelerations of type "+indstr+" in restPart" << std::endl;
 		          exit(EXIT_FAILURE);
		      }	
		}

		if (stat_mass == 0){	     		
 	           ret = H5Sselect_hyperslab(part_Dspace_1D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_short_vert_dset);	
	           ret = H5Sselect_hyperslab(part_Buffspace_1D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_short_vert_buff);	
	           ret = H5Dread(pMass, H5T_NATIVE_FLOAT, part_Buffspace_1D, part_Dspace_1D, H5P_DEFAULT, myMass_short); 	
		   if (ret<0){
			  std::cout << "Cannot read masses of type "+indstr+" in restPart" << std::endl;
 		          exit(EXIT_FAILURE);
		      }		     
		}

		if (stat_pot == 0){
 	           ret = H5Sselect_hyperslab(part_Dspace_1D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_short_vert_dset);	
	           ret = H5Sselect_hyperslab(part_Buffspace_1D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_short_vert_buff);	
	           ret = H5Dread(pPot, H5T_NATIVE_FLOAT, part_Buffspace_1D, part_Dspace_1D, H5P_DEFAULT, myPot_short);		 
		   if (ret<0){
			  std::cout << "Cannot read potentials of type "+indstr+" in restPart" << std::endl;
 		          exit(EXIT_FAILURE);
		      }     		
		}

		if (stat_IDs == 0){
 	           ret = H5Sselect_hyperslab(part_Dspace_1D, H5S_SELECT_SET, start_dset, stride_dset, count_dset, block_short_vert_dset);	
	           ret = H5Sselect_hyperslab(part_Buffspace_1D, H5S_SELECT_SET, start_buff, stride_buff, count_buff, block_short_vert_buff);	
	           ret = H5Dread(pIDs, H5T_NATIVE_INT, part_Buffspace_1D, part_Dspace_1D, H5P_DEFAULT, myIDs_short);	
		   if (ret<0){
			  std::cout << "Cannot read IDs of type "+indstr+" in restPart" << std::endl;
 		          exit(EXIT_FAILURE);
		      }	      
		}

		int start = fracPart*10;
                int end   = fracPart*10+restPart;	

		myFunc_Assign_H5parts(pt0,pt1,pt2,pt3,pt4,pt5,&grindex,&start,&end,myCoord_short,myVels_short,myAcc_short,myMass_short,myPot_short,myIDs_short);

             } //Close if on restpart

	     H5Dclose(pCoord);	
	     H5Dclose(pVels);
             H5Dclose(pAcc);
	     H5Dclose(pMass);
	     H5Dclose(pPot);
             H5Dclose(pIDs);
	
	     H5Gclose(myPartRead);	

           }  //This closes if on check of how many particles are in particle table at index grindex 


	} //this closes loop on grindex = 0 -> 6

	H5Fclose(myFileRead);
		
	std::cout <<  "File " << hp->myHDF5File << " has been read..." << std::endl;

}
	

//Function that assigns arrays to particle structures

void myFunc_Assign_H5parts(std::vector<struct h5_particle> *pt0, std::vector<struct h5_particle> *pt1, std::vector<struct h5_particle> *pt2, std::vector<struct h5_particle> *pt3, std::vector<struct h5_particle> *pt4, std::vector<struct h5_particle> *pt5, int *grindex,int *start, int*end, float myCoord[][3], float myVels[][3], float myAcc[][3],  float myMass[],  float myPot[], int myIDs[]){


		for(int myVal = *start;myVal<*end; myVal++){

		      int *my_xx = &myVal; 	

		      if (*grindex==0){
  		          for(int my_yy = 0;my_yy <3; my_yy++){	
   
			      pt0->at(*my_xx).pos[my_yy] = myCoord[(*my_xx)%10][my_yy];		
			      pt0->at(*my_xx).vel[my_yy] = myVels[(*my_xx)%10][my_yy];				   	      
			      pt0->at(*my_xx).acc[my_yy] = myAcc[(*my_xx)%10][my_yy];				   	      
			  }    	
   			   	
        		  pt0->at(*my_xx).mass = myMass[(*my_xx)%10];		
   		          pt0->at(*my_xx).phi = myPot[(*my_xx)%10];				   	      
		          pt0->at(*my_xx).idf = myIDs[(*my_xx)%10];				   	      
		       
		      }else if (*grindex==1){
  		          for(int my_yy = 0;my_yy <3; my_yy++){	
   
			      pt1->at(*my_xx).pos[my_yy] = myCoord[(*my_xx)%10][my_yy];		
			      pt1->at(*my_xx).vel[my_yy] = myVels[(*my_xx)%10][my_yy];				   	      
			      pt1->at(*my_xx).acc[my_yy] = myAcc[(*my_xx)%10][my_yy];				   	      
			  }    	
   			   	
        		  pt1->at(*my_xx).mass = myMass[(*my_xx)%10];		
   		          pt1->at(*my_xx).phi = myPot[(*my_xx)%10];				   	      
		          pt1->at(*my_xx).idf = myIDs[(*my_xx)%10];				   	      
		       
		      }else if (*grindex==2){
  		          for(int my_yy = 0;my_yy <3; my_yy++){	
   
			      pt2->at(*my_xx).pos[my_yy] = myCoord[(*my_xx)%10][my_yy];		
			      pt2->at(*my_xx).vel[my_yy] = myVels[(*my_xx)%10][my_yy];				   	      
			      pt2->at(*my_xx).acc[my_yy] = myAcc[(*my_xx)%10][my_yy];				   	      
			  }    	
   			   	
        		  pt2->at(*my_xx).mass = myMass[(*my_xx)%10];		
   		          pt2->at(*my_xx).phi = myPot[(*my_xx)%10];				   	      
		          pt2->at(*my_xx).idf = myIDs[(*my_xx)%10];				   	      
		       
		      }else if (*grindex==3){
  		          for(int my_yy = 0;my_yy <3; my_yy++){	
   
			      pt3->at(*my_xx).pos[my_yy] = myCoord[(*my_xx)%10][my_yy];		
			      pt3->at(*my_xx).vel[my_yy] = myVels[(*my_xx)%10][my_yy];				   	      
			      pt3->at(*my_xx).acc[my_yy] = myAcc[(*my_xx)%10][my_yy];				   	      
			  }    	
   			   	
        		  pt3->at(*my_xx).mass = myMass[(*my_xx)%10];		
   		          pt3->at(*my_xx).phi = myPot[(*my_xx)%10];				   	      
		          pt3->at(*my_xx).idf = myIDs[(*my_xx)%10];				   	      
		       
		      }else if (*grindex==4){
  		          for(int my_yy = 0;my_yy <3; my_yy++){	
   
			      pt4->at(*my_xx).pos[my_yy] = myCoord[(*my_xx)%10][my_yy];		
			      pt4->at(*my_xx).vel[my_yy] = myVels[(*my_xx)%10][my_yy];				   	      
			      pt4->at(*my_xx).acc[my_yy] = myAcc[(*my_xx)%10][my_yy];				   	      
			  }    	
   			   	
        		  pt4->at(*my_xx).mass = myMass[(*my_xx)%10];		
   		          pt4->at(*my_xx).phi = myPot[(*my_xx)%10];				   	      
		          pt4->at(*my_xx).idf = myIDs[(*my_xx)%10];				   	      
		       
		      }else if (*grindex==5){
  		          for(int my_yy = 0;my_yy <3; my_yy++){	
   
			      pt5->at(*my_xx).pos[my_yy] = myCoord[(*my_xx)%10][my_yy];		
			      pt5->at(*my_xx).vel[my_yy] = myVels[(*my_xx)%10][my_yy];				   	      
			      pt5->at(*my_xx).acc[my_yy] = myAcc[(*my_xx)%10][my_yy];				   	      
			  }    	
   			   	
        		  pt5->at(*my_xx).mass = myMass[(*my_xx)%10];		
   		          pt5->at(*my_xx).phi = myPot[(*my_xx)%10];				   	      
		          pt5->at(*my_xx).idf = myIDs[(*my_xx)%10];				   	      
		       
		       }

		 }   

	}

































