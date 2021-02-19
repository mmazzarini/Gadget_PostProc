#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include "library_MMazzarini_structs.h"
#include "library_MMazzarini_readwriteHDF5.h"
#include "library_MMazzarini_utilities.h"


void read_sbox_bin_MW(std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2, std::vector<struct h5_particle> *pp3, struct header_h5 *hp, int *myTime){

	float f;
	int count, ind_disk, ind_bulge, ind_halo, ind_file;	

	std::string myTimeNumb;

	if(*myTime<10){

		myTimeNumb = "0000"+myString(*myTime);
	
	}else if(*myTime >= 10 and *myTime < 100){

		myTimeNumb = "000"+myString(*myTime);

	}else if(*myTime >= 100 and *myTime < 1000){

		myTimeNumb = "00"+myString(*myTime);

	}else if(*myTime >= 1000 and *myTime < 10000){

		myTimeNumb = "0"+myString(*myTime);

	}else if(*myTime >= 10000){

		myTimeNumb = myString(*myTime);

	}	


	system(("rsync -a sfb069@milkyway.zah.uni-heidelberg.de:/home0/sfb042/superbox_tests/REZA/dbhs/"+hp->myHDF5Root+"/"+hp->myHDF5Root+"-g01."+myTimeNumb+" ./").c_str());

        std::ifstream fin((hp->myHDF5Root+"-g01."+myTimeNumb).c_str(), std::ios::binary);

	struct h5_particle pp_provv;

	//ind = 0;
	ind_disk = 0;
	ind_bulge = 0;
	ind_halo = 0;
	ind_file = 0;

	while( !fin.eof() ){

		for (count=0;count<7;count++){

			if(count >= 0 and count <3){

				fin.read(reinterpret_cast<char*>(&f), sizeof(float));
				if (fin.eof()) break;	
					
				pp_provv.pos[count] = f;

			}else if(count >= 3 and count <6){

				fin.read(reinterpret_cast<char*>(&f), sizeof(float));
				if (fin.eof()) break;	
					
				pp_provv.vel[count-3] = f;

			}else if(count == 6){

				fin.read(reinterpret_cast<char*>(&f), sizeof(float));
				if (fin.eof()) break;	
					
				pp_provv.mass = f;
				
				//if(ind_file <11){

					//std::cout << "the provvisory mass is::::::::: " << pp_provv.mass << std::endl;

				//}
			}

		}		

		for (count=0;count<5;count++){

			fin.read(reinterpret_cast<char*>(&f), sizeof(float));
			if (fin.eof()) break;				

		}

		if(ind_file < 1e7){
		
			if( std::abs(pp_provv.pos[0]) < 1e20 and std::abs(pp_provv.pos[1]) < 1e20 and std::abs(pp_provv.pos[2]) < 1e20 ){
			
				pp1->push_back(pp_provv);		
				ind_disk++;
			}	

		}else if(ind_file >= 1e7 and ind_file <10500000){
		
			if( std::abs(pp_provv.pos[0]) < 1e20 and std::abs(pp_provv.pos[1]) < 1e20 and std::abs(pp_provv.pos[2]) < 1e20 ){
			
				pp2->push_back(pp_provv);
				ind_bulge++;		

			}

		}else if(ind_file >= 10500000 and ind_file <14500000){
		
			if( std::abs(pp_provv.pos[0]) < 1e20 and std::abs(pp_provv.pos[1]) < 1e20 and std::abs(pp_provv.pos[2]) < 1e20){
			
				pp3->push_back(pp_provv);	
				ind_halo++;	

			}

		}
	
		ind_file++;	

	}

	std::cout << "trota salmonata" << std::endl;

	fin.close();

	system(("rm "+hp->myHDF5Root+"-g01."+myTimeNumb).c_str());

	hp->nbody1 = pp1->size();
	hp->nbody2 = pp2->size();
	hp->nbody3 = pp3->size();

	std::cout << "Read disk " << hp->nbody1 << std::endl;  
	std::cout << "Read bulge " << hp->nbody2 << std::endl;  
	std::cout << "Read halo " << hp->nbody3 << std::endl;  

}


void read_sbox_bin_sat(std::vector<struct h5_particle> *ppsat, struct header_h5 *hp, int *myTime, int *sat_ind){

	std::cout << "Reading satellite " << *sat_ind << std::endl;

	float f;
	int count, ind, ind_sat, ind_file;

	std::string myTimeNumb, myNumString, myFileName;

	if(*myTime<10){

		myTimeNumb = "0000"+myString(*myTime);
	
	}else if(*myTime >= 10 and *myTime < 100){

		myTimeNumb = "000"+myString(*myTime);

	}else if(*myTime >= 100 and *myTime < 1000){

		myTimeNumb = "00"+myString(*myTime);

	}else if(*myTime >= 1000 and *myTime < 10000){

		myTimeNumb = "0"+myString(*myTime);

	}else if(*myTime >= 10000){

		myTimeNumb = myString(*myTime);

	}

	
	myNumString = myString(*sat_ind);

	if(*sat_ind < 10){

		myFileName = hp->myHDF5Root+"-g0"+myNumString+"."+myTimeNumb; 

	}else{

		myFileName = hp->myHDF5Root+"-g"+myNumString+"."+myTimeNumb; 

	}

	system(("rsync -a sfb069@milkyway.zah.uni-heidelberg.de:/home0/sfb042/superbox_tests/REZA/dbhs/"+hp->myHDF5Root+"/"+myFileName+" ./").c_str());

	std::ifstream fin(myFileName.c_str(), std::ios::binary);

	struct h5_particle pp_provv;

	ind = 0;
	ind_sat = 0;
	ind_file = 0;

	while( !fin.eof() ){

		for (count=0;count<7;count++){

			if(count >= 0 and count <3){

				fin.read(reinterpret_cast<char*>(&f), sizeof(float));
				if (fin.eof()) break;	
					
				pp_provv.pos[count] = f;

			}else if(count >= 3 and count <6){

				fin.read(reinterpret_cast<char*>(&f), sizeof(float));
				if (fin.eof()) break;	
					
				pp_provv.vel[count-3] = f;

			}else if(count == 6){

				fin.read(reinterpret_cast<char*>(&f), sizeof(float));
				if (fin.eof()) break;	
					
				pp_provv.mass = f;
			
				//if(ind_sat <11){

					//std::cout << "the provvisory satellite mass is::::::::: " << pp_provv.mass << std::endl;

				//}

			}

		}		

		for (count=0;count<5;count++){

			fin.read(reinterpret_cast<char*>(&f), sizeof(float));
			if (fin.eof()) break;				

		}

		if(std::abs(pp_provv.pos[0]) < 1e20 and std::abs(pp_provv.pos[1]) < 1e20 and std::abs(pp_provv.pos[2]) < 1e20 and ind_file<50000){

			ppsat->push_back(pp_provv);		
			ind_sat++;

		}

		ind_file++;

	}

	fin.close();

	system(("rm "+myFileName).c_str());
	
	//nbody 

	hp->nbodySat[*sat_ind-2] = ind_sat;	

}


