#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include "library_MMazzarini_structs.h"
#include "library_MMazzarini_utilities.h"
#include "library_MMazzarini_density.h"



void read_densities(struct mass_center *dcms,std::vector<struct mass_center> *dcms_Sat, struct header_h5 *hp, int *myTime){

	std::ifstream myDensMWfile;
	std::ifstream myDensSatFile;	
	std::string myDensMWline;
	std::string myDensSatLine;
	int myMWindex = 0;

	

	std::string myMWTitle = "dcms-"+hp->myHDF5Root+"-g01.txt";
	
	myDensMWfile.open(myMWTitle.c_str());
	
	std::cout << "Reading MW file " <<  myMWTitle << std::endl;
	std::cout << "DEBUG-MESSAGE for Milky Way" << std::endl;	

        if (!myDensMWfile){

                std::cout << "ERROR! THE FILE IS NOT FOUND" << std::endl;
                exit(1);
        }

	while(std::getline(myDensMWfile,myDensMWline)){

		float trash;
		if ( myMWindex == (*myTime) ){

			std::istringstream myDensMWstream(myDensMWline.c_str());
			myDensMWstream >> trash;
			myDensMWstream >> dcms->pos[0];			
			myDensMWstream >> dcms->pos[1];			
			myDensMWstream >> dcms->pos[2];			
			myDensMWstream >> dcms->vel[0];			
			myDensMWstream >> dcms->vel[1];			
			myDensMWstream >> dcms->vel[2];			
				
		}

		++myMWindex;

	}

	myDensMWfile.close();

	std::cout << "MW dcms is " << dcms->pos[0] << " " << dcms->pos[1] << " " << dcms->pos[2] << std::endl;

	for (int i=2;i<hp->sat_numb+2;i++){

		int mySatIndex = 0;		
		std::string iStr = myString(i);
		std::cout << "string satellite is " << iStr << std::endl;
		if (i < 10){
			std::string mySatTitle = "dcms-"+hp->myHDF5Root+"-g0"+iStr+".txt";
			myDensSatFile.open(mySatTitle.c_str()); 
			std::cout << "Reading satellite file " << mySatTitle << std::endl;
		}else{

			std::string mySatTitle = "dcms-"+hp->myHDF5Root+"-g"+iStr+".txt";
			myDensSatFile.open(mySatTitle.c_str()); 
			std::cout << "Reading satellite file " <<  mySatTitle << std::endl;

		}
		std::cout << "DEBUG MESSAGE FOR satellite " << i << std::endl;

	        if (!myDensSatFile){

	                std::cout << "ERROR! THE FILE IS NOT FOUND" << std::endl;
        	        exit(1);
	        
		}

		while(std::getline(myDensSatFile,myDensSatLine)){

			float trash;
			if ( mySatIndex == (*myTime) ){

				std::istringstream myDensSatStream(myDensSatLine.c_str());
				myDensSatStream >> trash;
				myDensSatStream >> dcms_Sat->at(i-2).pos[0];			
				myDensSatStream >> dcms_Sat->at(i-2).pos[1];			
				myDensSatStream >> dcms_Sat->at(i-2).pos[2];			
				myDensSatStream >> dcms_Sat->at(i-2).vel[0];			
				myDensSatStream >> dcms_Sat->at(i-2).vel[1];			
				myDensSatStream >> dcms_Sat->at(i-2).vel[2];			

				std::cout << "Reading line " << mySatIndex << std::endl; 
				
			}
			
			++mySatIndex;
			
		}			

		myDensSatFile.close();
	
		std::cout << "satellite dcms is " << dcms_Sat->at(i-2).pos[0] << " " << dcms_Sat->at(i-2).pos[1] << " " << dcms_Sat->at(i-2).pos[2] << std::endl;

	}

}


