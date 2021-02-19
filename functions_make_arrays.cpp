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
#include "library_MMazzarini_arrays.h"
#include "library_MMazzarini_utilities.h"

using namespace H5;

//Make position arrays for particles

void kernel_general(struct header_h5 *hp, struct array_hd *ap){  

	std::string inRoot;
	std::string inClass;
	std::string inFirst;
	std::string inLast;
	std::string inStep;
	std::string inBaryons;
	std::string inSat;
	std::string inCC;

	std::cout << "Enter files root: ";	 	       
 	std::getline(std::cin, inRoot);
	hp->myHDF5Root = inRoot;

        std::cout << "Enter Aq type for Table of numbers (e.g. A2, B2...): ";
        std::getline(std::cin, inClass);
        hp->AquariusClass = inClass;

	std::cout << "Will read " << hp->myHDF5Root << std::endl;

	std::cout << "Enter first output file number: ";	 	       
	std::getline(std::cin, inFirst);
	std::stringstream myStreamFirst(inFirst);
	myStreamFirst >> hp->first;
	
	std::cout << "Enter last output file number: ";	 	       
	std::getline(std::cin, inLast);
	std::stringstream myStreamLast(inLast);
	myStreamLast >> hp->last;

	std::cout << "Enter file step number: ";	 	       
	std::getline(std::cin, inStep);
	std::stringstream myStreamStep(inStep);
	myStreamStep >> hp->step;

	std::cout << "Do you want to include baryons (0/1 = no/yes)? ";	 	       
	std::getline(std::cin, inBaryons);
	std::stringstream myStreamBaryons(inBaryons);
	myStreamBaryons >> hp->wantBaryons;

	std::cout << "How many satellites you have? ";	 	       
	std::getline(std::cin, inSat);
	std::stringstream myStreamNumSat(inSat);
	myStreamNumSat >> hp->sat_numb;

	std::cout << "What is the c param of the MW halo? ";
	std::getline(std::cin, inCC);
	std::stringstream myStreamConcentr(inCC);
	myStreamConcentr >> hp->MW_concentration;

	std::cout << "We will work with a number of resolution values equal to hp->sat_numb+6" << std::endl;

	std::cout << "Segnaposto 1" << std::endl;

	array_kernel(ap);

	std::cout << "Segnaposto 2" << std::endl;

	ap->sat_numb = hp->sat_numb;
	
	std::cout << "Segnaposto 3" << std::endl;

	ap->res.resize(6+ap->sat_numb);
	
	std::cout << "Segnaposto 4" << std::endl;

	std::cout << "qui" << std::endl;

	ap->res[1] = 0.05;
	ap->res[2] = 0.05;
	ap->res[3] = 0.05;
	ap->res[4] = 0.01;	
	ap->res[5] = 0.01;	

	for (int iRes=6;iRes<6+hp->sat_numb;iRes++){

		ap->res[iRes] = 0.025;

	}

	/*

	std::cout << "\nOpening file " << ("softenings_table_"+hp->myHDF5Root+".txt").c_str() << std::endl;
	
	std::string mySoftLine;
	std::ifstream mySoftFile;
	mySoftFile.open(("softenings_table_"+hp->myHDF5Root+".txt").c_str());
	int iRes;
	iRes = 6;
	
	while(std::getline(mySoftFile,mySoftLine)){
		std::cout << "qui debug" << iRes << std::endl;
		float trash;
		std::istringstream mySoftStream(mySoftLine.c_str());
		mySoftStream >> trash;
		mySoftStream >> trash;
		mySoftStream >> ap->res[iRes];
		++iRes;

	}
	mySoftFile.close();
	std::cout << "Qui fine lettura softening" << std::endl;

	*/

}		

void array_kernel(struct array_hd *ap){
	
	int xlen, ylen, zlen; 	
	std::string inRes;
	std::string inXarr = "";
	std::string inYarr = "";
	std::string inZarr = "";
	bool xStatus = false;
	bool yStatus = false;
	bool zStatus = false;

	while(xStatus == false){

		std::cout << "Enter map semi-extension(in kpc)	along x (e.g. 100.): ";	 	       
		std::getline(std::cin, inXarr);
		std::stringstream myStreamXarr(inXarr);
		myStreamXarr >> xlen;
		if(xlen < 0.){

			std::cout << "Invalid x value, must be >=0.";	

		}
		else{
			xStatus = true;
		}

	}
	
	while(yStatus == false){

		std::cout << "Enter map semi-extension(in kpc)	along y (e.g. 100.): ";	 	       
		std::getline(std::cin, inYarr);
		std::stringstream myStreamYarr(inYarr);
		myStreamYarr >> ylen;
		if(ylen < 0.){

			std::cout << "Invalid y value, must be >=0.";	

		}
		else{
			yStatus = true;
		}

	}	

	while(zStatus == false){

		std::cout << "Enter map semi-extension(in kpc)	along z (e.g. 100.): ";	 	       
		std::getline(std::cin, inZarr);
		std::stringstream myStreamZarr(inZarr);
		myStreamZarr >> zlen;
		if(zlen < 0.){

			std::cout << "Invalid z value, must be >=0.";	

		}
		else{
			zStatus = true;
		}

	}

	/*	

	for(int i=0;i<3;i++){

		bool softStatus = false;

		while(softStatus == false){

			std::cout << "Enter resolution for particle type " << i << " : ";	 	       
			std::getline(std::cin, inRes);
			std::stringstream myStreamRes(inRes);		
			myStreamRes >> ap->res[i];
	
			if(ap->res[i] < 0.){
	
				std::cout << "Invalid resolution, must be >= 0." << std::endl;

			}else{

				softStatus = true;

			} //End if

		} //End while	

	} //End for loop


	

	std::cout << "Now enter softenings for satellites..." << std::endl;

	for(int i=6;i6;i++){

		bool softStatus = false;

		while(softStatus == false){

			std::cout << "Enter resolution for satellite " << i << " : ";	 	       
			std::getline(std::cin, inRes);
			std::stringstream myStreamRes(inRes);		
			myStreamRes >> ap->res[i];
	
			if(ap->res[i] < 0.){
	
				std::cout << "Invalid resolution, must be >= 0." << std::endl;

			}else{

				softStatus = true;

			} //End if

		} //End while	

	} //End for loop
	
	
	
	*/

      /*************************************/

	/*


	ap->xMin_rel = -xlen;
	ap->xMax_rel = +xlen;
	ap->yMin_rel = -ylen;
	ap->yMax_rel = +ylen;	
	ap->zMin_rel = -zlen;
	ap->zMax_rel = +zlen;

	for(int i = 0;i<6;i++){
		if(ap->res[i] > 0){
			ap->nX[i] = (int)((ap->xMax_rel - ap->xMin_rel)/(ap->res[i]));
			ap->nY[i] = (int)((ap->yMax_rel - ap->yMin_rel)/(ap->res[i]));	
			ap->nZ[i] = (int)((ap->zMax_rel - ap->zMin_rel)/(ap->res[i]));	
		}else{
			ap->nX[i] = 0;
			ap->nY[i] = 0;
			ap->nZ[i] = 0;			
		}
		
	}


	*/
	
}




// The following: to be called for(int i_part=0;i_part<6;i_part++){

void make_xyz_arrays(std::string myHdf5File, struct mass_center *dcms, std::vector<struct h5_particle> *pt, struct header_h5 *hp,struct array_hd *ap,int *pType){

	//	Calculate array for each particle type, using the softening

	int nPart = pt->size();
	int iPart;
	int ix;
	int iy;
	int iz;
	int RANK_ARR = 2;
	herr_t stat;
	hsize_t array_dims[2];
	hid_t arraySpace;	

	int nX = ap->nX[*pType];
	int nY = ap->nY[*pType]; 
	int nZ = ap->nZ[*pType];

	float mXY[nY][nX];
	float mZY[nY][nZ];
	float mXZ[nZ][nX];
	
	for (iPart=0;iPart<nPart;iPart++){

		//xy array

		for(ix=0; ix< nX;ix++){

			if(  ( (ap->xMin_rel + ix*(ap->res[*pType]))< pt->at(iPart).pos[0]) and (pt->at(iPart).pos[0] <= (ap->xMin_rel + (ix+1.)*(ap->res[*pType])) )  ){ 			
	
				for(iy=0;iy< nY;iy++){
				
					if(  ( (ap->yMin_rel + iy*(ap->res[*pType]))< pt->at(iPart).pos[1]) and (pt->at(iPart).pos[1] <= (ap->yMin_rel + (iy+1.)*(ap->res[*pType])) )  ){

						mXY[iy][ix] += pt->at(iPart).mass;
					}
				}   
			}
		}

		//yz array

		for(iy=0; iy< nY;iy++){

			if(  ( (ap->yMin_rel + iy*(ap->res[*pType]))< pt->at(iPart).pos[1]) and (pt->at(iPart).pos[1] <= (ap->yMin_rel + (iy+1.)*(ap->res[*pType])) )  ){ 			
	
				for(iz=0;iz< nZ;iz++){
				
					if(  ( (ap->zMin_rel + iz*(ap->res[*pType]))< pt->at(iPart).pos[2]) and (pt->at(iPart).pos[2] <= (ap->zMin_rel + (iz+1.)*(ap->res[*pType])) )  ){

						mZY[iy][iz] += pt->at(iPart).mass;
					}
				}   
			}
		}

		//zx array
		
		for(iz=0; iz< nZ;iz++){

			if(  ( (ap->zMin_rel + iz*(ap->res[*pType]))< pt->at(iPart).pos[2]) and (pt->at(iPart).pos[2] <= (ap->zMin_rel + (iz+1.)*(ap->res[*pType])) )  ){ 			
	
				for(ix=0;ix< nX;ix++){
				
					if(  ( (ap->xMin_rel + ix*(ap->res[*pType]))< pt->at(iPart).pos[0]) and (pt->at(iPart).pos[0] <= (ap->xMin_rel + (ix+1.)*(ap->res[*pType])) )  ){

						mXZ[iz][ix] += pt->at(iPart).mass;	
					}
				}   
			}
		}


	}	       		

	hid_t myH5arr = H5Fopen((hp->myHDF5Maps).c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
	hid_t myGroup = H5Gopen(myH5arr,("/PartType"+myString(*pType)).c_str(),H5P_DEFAULT);

	array_dims[0] = nY;
	array_dims[1] = nX;
	arraySpace = H5Screate_simple(RANK_ARR,array_dims,NULL);	
	hid_t myArrXY = H5Dcreate(myGroup,"arrayXY", H5T_NATIVE_FLOAT, arraySpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	stat = H5Dwrite(myArrXY, H5T_NATIVE_FLOAT, arraySpace, arraySpace, H5P_DEFAULT, mXY);	

	array_dims[0] = nY;
	array_dims[1] = nZ;
	arraySpace = H5Screate_simple(RANK_ARR,array_dims,NULL);
	hid_t myArrYZ = H5Dcreate(myGroup,"arrayYZ", H5T_NATIVE_FLOAT, arraySpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	stat = H5Dwrite(myArrYZ, H5T_NATIVE_FLOAT, arraySpace, arraySpace, H5P_DEFAULT, mZY);

	array_dims[0] = nZ;
	array_dims[1] = nX;
	arraySpace = H5Screate_simple(RANK_ARR,array_dims,NULL);
	hid_t myArrZX = H5Dcreate(myGroup,"arrayZX", H5T_NATIVE_FLOAT, arraySpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	stat = H5Dwrite(myArrZX, H5T_NATIVE_FLOAT, arraySpace, arraySpace, H5P_DEFAULT, mXZ);		

	stat = H5Gclose(myGroup);	
	stat = H5Fclose(myH5arr);

	std::cout << "      Array stored in file " << hp->myHDF5Maps+"/PartType"+(myString(*pType)).c_str() << std::endl;

}


