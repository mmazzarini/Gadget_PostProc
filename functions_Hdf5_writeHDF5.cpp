#include <stdio.h>
#include <stdlib.h>
#include <iostream>
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

/*
*
* Converts tipsy particles to HDF5 arrays
*
*/

void write_Hdf5_postproc(std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2, std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4,std::vector<struct h5_particle> *pp5,struct header_h5 *hp, struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, struct mass_center *dcms_Baryon, struct mass_center *dcms_Dark){
 

	/* Variables and arrays  */

	int dim_short  =  8;
	int dim_long   = 10;
	int dim_large  = 40;
        int Dim2 = 4;   
	int Dim3 = 3;
	herr_t ret;      

	int nPart_high_table[6]; 
	int nPart_table[6];
	int mass_table[6];
	float dcms_MW_arr[6];
	float dcms_Sat_arr[6*hp->sat_numb];
	float dcms_Baryon_arr[6];
	float dcms_Dark_arr[6];	
	float dcms_AngMom_arr[6];
	float InertiaDark_arr[3][3][50];
	float InertiaDark_Unbound_arr[3][3][50];
        float InertiaStar_arr[3][3][50];
        float InertiaStar_Unbound_arr[3][3][50];
        float InertiaHalo_arr[3][3][50];
        float InertiaHalo_Unbound_arr[3][3][50];
	float geomDark_Unbound_arr[3][3][50];
        float geomStar_Unbound_arr[3][3][50];
        float geomHalo_arr[3][3][50];



	float sat_massBound_Arr[hp->sat_numb];
	float sat_massBound_p4_Arr[hp->sat_numb];
	float sat_massBound_p5_Arr[hp->sat_numb];
	float sat_totMass_Arr[hp->sat_numb];
	float sat_totMass_p4_Arr[hp->sat_numb];
	float sat_totMass_p5_Arr[hp->sat_numb];
	float sat_outTidMass_Arr[hp->sat_numb];
	float sat_outTidMass_p4_Arr[hp->sat_numb];
	float sat_outTidMass_p5_Arr[hp->sat_numb];
	float MW_integratedMass_Arr[hp->sat_numb];
	float sat_tidRadius_Arr[hp->sat_numb];
	float sat_tidalMass_Arr[hp->sat_numb];
	float sat_tidalMass_p4_Arr[hp->sat_numb];
	float sat_tidalMass_p5_Arr[hp->sat_numb];
	float sat_ecc_Arr[hp->sat_numb];
	float sat_ax_Arr[hp->sat_numb];

	for (int i=0;i<3;i++){

		for(int j=0;j<3;j++){

			for (int iR=0;iR<50;iR++){		
	
				InertiaDark_arr[i][j][iR] = hp->inertiaTensorDark[i][j][iR];
				InertiaDark_Unbound_arr[i][j][iR] = hp->inertiaTensorDark_Unbound[i][j][iR];
        	                InertiaStar_arr[i][j][iR] = hp->inertiaTensorStar[i][j][iR];
        	                InertiaStar_Unbound_arr[i][j][iR] = hp->inertiaTensorStar_Unbound[i][j][iR];
        	                InertiaHalo_arr[i][j][iR] = hp->inertiaTensorHalo[i][j][iR];
        	                InertiaHalo_Unbound_arr[i][j][iR] = hp->inertiaTensorHalo_Unbound[i][j][iR];

				geomDark_Unbound_arr[i][j][iR] = hp->geomTensorDark_Unbound[i][j][iR];
        	                geomStar_Unbound_arr[i][j][iR] = hp->geomTensorStar_Unbound[i][j][iR];
        	                geomHalo_arr[i][j][iR] = hp->geomTensorHalo[i][j][iR];


			}

		}

	}

	for(int i=0;i<hp->sat_numb;i++){

		sat_massBound_Arr[i] = hp->sat_massBound[i];
		sat_massBound_p4_Arr[i] = hp->sat_massBound_p4[i];
		sat_massBound_p5_Arr[i] = hp->sat_massBound_p5[i];
		sat_totMass_Arr[i] = hp->sat_totMass[i];
		sat_totMass_p4_Arr[i] = hp->sat_totMass_p4[i];
		sat_totMass_p5_Arr[i] = hp->sat_totMass_p5[i];
		sat_outTidMass_Arr[i] = hp->sat_outTidMass[i];
		sat_outTidMass_p4_Arr[i] = hp->sat_outTidMass_p4[i];
		sat_outTidMass_p5_Arr[i] = hp->sat_outTidMass_p5[i];
		MW_integratedMass_Arr[i] = hp->MW_integratedMass[i];
		sat_tidRadius_Arr[i] = hp->sat_tidRadius[i];
		sat_tidalMass_Arr[i] = hp->sat_tidalMass[i];
		sat_tidalMass_p4_Arr[i] = hp->sat_tidalMass_p4[i];
		sat_tidalMass_p5_Arr[i] = hp->sat_tidalMass_p5[i];
		sat_ecc_Arr[i] = hp->sat_ecc[i];
		sat_ax_Arr[i] = hp->sat_ax[i];
	
		for (int j=0;j<3;j++){

			dcms_Sat_arr[6*i+j] = dcms_Sat->at(i).pos[j];
			dcms_Sat_arr[6*i+j+3] = dcms_Sat->at(i).vel[j];

		}

	}	

	int dim_dcms = 6*hp->sat_numb;

        /* Create dataSpaces for Groups */

	hsize_t short_Dims_1D; //8
	hsize_t long_Dims_1D;  //20
	hsize_t large_Dims_1D;  //40
	hsize_t hiRes_Dims_1D; //50
	hsize_t hiRes2_Dims_1D; //250
	hsize_t hiRes0_Dims_1D; //25

	hsize_t multi_dims[2];
	hsize_t dens_dims_1D;
	hsize_t densSat_dims_1D;
	hsize_t sat_dims_1D;
	hsize_t sat_dims_large_1D;
	hsize_t sat_dims_hiRes_1D;
	hsize_t sat_dims_hiRes2_1D;
	hsize_t sat_dims_hiRes0_1D;
	hsize_t inertia_dims_3D[3];

	short_Dims_1D = dim_short;
	long_Dims_1D = dim_long;
	large_Dims_1D = dim_large;
	hiRes_Dims_1D = 50;
	hiRes2_Dims_1D = 250;
	hiRes0_Dims_1D = 25;
	multi_dims[0] = dim_large;
	multi_dims[1] = Dim2;
	dens_dims_1D = 6;
	densSat_dims_1D = dim_dcms;
	sat_dims_1D = hp->sat_numb;
	sat_dims_large_1D = 40*hp->sat_numb;
	sat_dims_hiRes_1D = 50*hp->sat_numb;
	sat_dims_hiRes2_1D = 250*hp->sat_numb;
	sat_dims_hiRes0_1D = 25*hp->sat_numb;
	inertia_dims_3D[0] = Dim3;
	inertia_dims_3D[1] = Dim3;
	inertia_dims_3D[2] = 50;	

       	hid_t short_Dspace_1D = H5Screate_simple(1,&short_Dims_1D,NULL);
	hid_t long_Dspace_1D  = H5Screate_simple(1,&long_Dims_1D,NULL);
       	hid_t large_Dspace_1D = H5Screate_simple(1,&large_Dims_1D,NULL);
	hid_t hRes_Dspace_1D = H5Screate_simple(1,&hiRes_Dims_1D,NULL);
	hid_t hRes2_Dspace_1D = H5Screate_simple(1,&hiRes2_Dims_1D,NULL);
	hid_t hRes0_Dspace_1D = H5Screate_simple(1,&hiRes0_Dims_1D,NULL);
	hid_t density_Dspace_1D = H5Screate_simple(1,&dens_dims_1D,NULL);
	hid_t densitySat_Dspace_1D = H5Screate_simple(1,&densSat_dims_1D,NULL);
	hid_t sat_Dspace_1D = H5Screate_simple(1,&sat_dims_1D,NULL);
	hid_t sat_Dspace_large_1D = H5Screate_simple(1,&sat_dims_large_1D,NULL);
	hid_t sat_Dspace_hRes_1D = H5Screate_simple(1,&sat_dims_hiRes_1D,NULL);
	hid_t sat_Dspace_hRes2_1D = H5Screate_simple(1,&sat_dims_hiRes2_1D,NULL);
	hid_t sat_Dspace_hRes0_1D = H5Screate_simple(1,&sat_dims_hiRes0_1D,NULL);
	hid_t Inertia_Dspace_3D = H5Screate_simple(3,inertia_dims_3D,NULL);
	hid_t multi_Dspace_2D = H5Screate_simple(2,multi_dims,NULL);

       	hid_t short_Buffspace_1D = H5Screate_simple(1,&short_Dims_1D,NULL);
	hid_t long_Buffspace_1D  = H5Screate_simple(1,&long_Dims_1D,NULL);
       	hid_t large_Buffspace_1D = H5Screate_simple(1,&large_Dims_1D,NULL);
	hid_t hRes_Buffspace_1D = H5Screate_simple(1,&hiRes_Dims_1D,NULL);
	hid_t hRes2_Buffspace_1D = H5Screate_simple(1,&hiRes2_Dims_1D,NULL);
	hid_t hRes0_Buffspace_1D = H5Screate_simple(1,&hiRes0_Dims_1D,NULL);
	hid_t density_Buffspace_1D = H5Screate_simple(1,&dens_dims_1D,NULL);
	hid_t densitySat_Buffspace_1D = H5Screate_simple(1,&densSat_dims_1D,NULL);
	hid_t sat_Buffspace_1D = H5Screate_simple(1,&sat_dims_1D,NULL);
	hid_t sat_Buffspace_large_1D = H5Screate_simple(1,&sat_dims_large_1D,NULL);
	hid_t sat_Buffspace_hRes_1D = H5Screate_simple(1,&sat_dims_hiRes_1D,NULL);
	hid_t sat_Buffspace_hRes2_1D = H5Screate_simple(1,&sat_dims_hiRes2_1D,NULL);
	hid_t sat_Buffspace_hRes0_1D = H5Screate_simple(1,&sat_dims_hiRes0_1D,NULL);
	hid_t Inertia_Buffspace_3D = H5Screate_simple(3,inertia_dims_3D,NULL);
	hid_t multi_Buffspace_2D = H5Screate_simple(2,multi_dims,NULL);

	/* Store data          */

	std::cout << "Storing data in " << hp->myHDF5Post << std::endl;
	hid_t my_write_File = H5Fcreate((hp->myHDF5Post).c_str(),H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);			   
	
	/*    Create groups; In this version, only part types from 1 to 5 are stored  */

	hid_t headGroup  = H5Gcreate(my_write_File, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t part0 = H5Gcreate(my_write_File, "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t part1 = H5Gcreate(my_write_File, "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t part2 = H5Gcreate(my_write_File, "/PartType2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t part3 = H5Gcreate(my_write_File, "/PartType3", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t part4 = H5Gcreate(my_write_File, "/PartType4", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t part5 = H5Gcreate(my_write_File, "/PartType5", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t densC = H5Gcreate(my_write_File, "/Density_Centers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

       /* Set arrays for Header Attributes      */

	std::cout << "Storing time: " << hp->time << std::endl;

	for (int i=0;i<6;i++){
	   mass_table[i] = 0.;
	   nPart_high_table[i] = 0;
	}

        nPart_table[0] = hp->nbody0;
        nPart_table[1] = hp->nbody1;
        nPart_table[2] = hp->nbody2;
        nPart_table[3] = hp->nbody3;
        nPart_table[4] = hp->nbody4;
        nPart_table[5] = hp->nbody5;

	/* Set dataSpaces for Header Attributes		*/

	hsize_t attr_dims_vect = 6;
	hsize_t attr_dims_scal = 1;
	
	hid_t vect_Aspace = H5Screate_simple(1,&attr_dims_vect,NULL);
	hid_t scal_Aspace = H5Screate_simple(1,&attr_dims_scal,NULL);

	/* Set arrays for density center */

	for(int i=0;i<3;i++){
		dcms_MW_arr[i] = dcms->pos[i];
		dcms_AngMom_arr[i] = dcms->angM[i];
		/*
		dcms_Sat_arr[i] = dcms_Sat->pos[i];
		dcms_Baryon_arr[i] = dcms_Baryon->pos[i];
		dcms_Dark_arr[i] = dcms_Dark->pos[i];
		 
		*/	
		dcms_MW_arr[i+3] = dcms->vel[i];
		dcms_AngMom_arr[i+3] = 0.;
		/*
		dcms_Sat_arr[i+3] = dcms_Sat->vel[i];
		dcms_Baryon_arr[i+3] = dcms_Baryon->vel[i];
		dcms_Dark_arr[i+3] = dcms_Dark->vel[i];		
		
		*/
	}




	float diskVzed_p1_Arr[40];
	float diskZed_p1_Arr[40];
	float diskVzed_p4_Arr[40];
	float diskZed_p4_Arr[40];
	float diskVzed_p5_Arr[40];
	float diskZed_p5_Arr[40];
	float diskSigVzed_p1_Arr[40];
	float diskSigZed_p1_Arr[40];
	float diskVtan_p1_Arr[40];
	float diskMass_p1_Arr[40];
	float diskMass_p4_Arr[40];
	float diskMass_p5_Arr[40];

	for(int i=0;i<40;i++){

		diskVzed_p1_Arr[i] = hp->disk_Vzed_p1[i];
		diskZed_p1_Arr[i] = hp->disk_Zed_p1[i];
		diskSigZed_p1_Arr[i] = hp->disk_SigZed_p1[i];
		diskSigVzed_p1_Arr[i] = hp->disk_SigVzed_p1[i];
		diskVtan_p1_Arr[i] = hp->disk_Vtan_p1[i];
		diskMass_p1_Arr[i] = hp->disk_Mass_p1[i];		
		diskMass_p4_Arr[i] = hp->disk_Mass_p4[i];		
		diskMass_p5_Arr[i] = hp->disk_Mass_p5[i];		
		diskVzed_p4_Arr[i] = hp->disk_Vzed_p4[i];
		diskZed_p4_Arr[i] = hp->disk_Zed_p4[i];
		diskVzed_p5_Arr[i] = hp->disk_Vzed_p5[i];
		diskZed_p5_Arr[i] = hp->disk_Zed_p5[i];

	}


	float bulgeVrad_p2_Arr[10];
	float bulgeVtheta_p2_Arr[10];
	float bulgeVphi_p2_Arr[10];
	float bulgeSigrad_p2_Arr[10];
	float bulgeSigtheta_p2_Arr[10];
	float bulgeSigphi_p2_Arr[10];
	float bulgeMass_p2_Arr[10];
	float bulgeMass_p4_Arr[10];
	float bulgeMass_p5_Arr[10];
	float bulgeVrad_p4_Arr[10];
	float bulgeVtheta_p4_Arr[10];
	float bulgeVphi_p4_Arr[10];
	float bulgeVrad_p5_Arr[10];
	float bulgeVtheta_p5_Arr[10];
	float bulgeVphi_p5_Arr[10];

	for(int i=0;i<10;i++){

		bulgeVrad_p2_Arr[i] = hp->bulge_Vrad_p2[i];
		bulgeVtheta_p2_Arr[i] = hp->bulge_Vtheta_p2[i];
		bulgeVphi_p2_Arr[i] = hp->bulge_Vphi_p2[i];
		bulgeVrad_p4_Arr[i] = hp->bulge_Vrad_p4[i];
		bulgeVtheta_p4_Arr[i] = hp->bulge_Vtheta_p4[i];
		bulgeVphi_p4_Arr[i] = hp->bulge_Vphi_p4[i];
		bulgeVrad_p5_Arr[i] = hp->bulge_Vrad_p5[i];
		bulgeVtheta_p5_Arr[i] = hp->bulge_Vtheta_p5[i];
		bulgeVphi_p5_Arr[i] = hp->bulge_Vphi_p5[i];
		bulgeSigrad_p2_Arr[i] = hp->bulge_Sigrad_p2[i];
		bulgeSigtheta_p2_Arr[i] = hp->bulge_Sigtheta_p2[i];		
		bulgeSigphi_p2_Arr[i] = hp->bulge_Sigphi_p2[i];		
		bulgeMass_p2_Arr[i] = hp->bulge_Mass_p2[i];
		bulgeMass_p4_Arr[i] = hp->bulge_Mass_p4[i];
		bulgeMass_p5_Arr[i] = hp->bulge_Mass_p5[i];

	}

	float halo_Vrad_p3_Arr[40];
	float halo_Vrad_p4_Arr[40];
	float halo_Vrad_p5_Arr[40];
	float halo_Vtheta_p3_Arr[40];
	float halo_Vtheta_p4_Arr[40];
	float halo_Vtheta_p5_Arr[40];
	float halo_Vphi_p3_Arr[40];
	float halo_Vphi_p4_Arr[40];
	float halo_Vphi_p5_Arr[40];
	float halo_SigVrad_p3_Arr[40];
	float halo_SigVrad_p4_Arr[40];
	float halo_SigVrad_p5_Arr[40];
	float halo_SigVtheta_p3_Arr[40];
	float halo_SigVtheta_p4_Arr[40];
	float halo_SigVtheta_p5_Arr[40];
	float halo_SigVphi_p3_Arr[40];
	float halo_SigVphi_p4_Arr[40];
	float halo_SigVphi_p5_Arr[40];
	float halo_ecc_p3_Arr[40];
	float halo_ecc_p4_Arr[40];
	float halo_ecc_p5_Arr[40];	
	float smj_ax_p3_Arr[40];
	float smj_ax_p4_Arr[40];
	float smj_ax_p5_Arr[40];
	float halo_Mass_p3_Arr[40];
	float halo_Mass_p4_Arr[40*hp->sat_numb];
	float halo_Mass_p5_Arr[40*hp->sat_numb];
	float halo_Mass_hiRes_p3_Arr[50];
	float halo_Mass_hiRes_p4_Arr[50*hp->sat_numb];
	float halo_Mass_hiRes_p5_Arr[50*hp->sat_numb];
	float halo_Mass_hiRes2_p3_Arr[250];
	float halo_Mass_hiRes2_p4_Arr[250*hp->sat_numb];
	float halo_Mass_hiRes2_p5_Arr[250*hp->sat_numb];
	float halo_Mass_hiRes0_p3_Arr[25];
	float halo_Mass_hiRes0_p4_Arr[25*hp->sat_numb];
	float halo_Mass_hiRes0_p5_Arr[25*hp->sat_numb];

	for(int i=0;i<40;i++){

		halo_Vrad_p3_Arr[i] = hp->halo_Vrad_p3[i];
		halo_Vrad_p4_Arr[i] = hp->halo_Vrad_p4[i];
		halo_Vrad_p5_Arr[i] = hp->halo_Vrad_p5[i];
		halo_Vtheta_p3_Arr[i] = hp->halo_Vtheta_p3[i];
		halo_Vtheta_p4_Arr[i] = hp->halo_Vtheta_p4[i];
		halo_Vtheta_p5_Arr[i] = hp->halo_Vtheta_p3[i];
		halo_Vphi_p3_Arr[i] = hp->halo_Vphi_p3[i];
		halo_Vphi_p4_Arr[i] = hp->halo_Vphi_p4[i];
		halo_Vphi_p5_Arr[i] = hp->halo_Vphi_p5[i];
		halo_SigVrad_p3_Arr[i] = hp->halo_SigVrad_p3[i];
		halo_SigVrad_p4_Arr[i] = hp->halo_SigVrad_p4[i];
		halo_SigVrad_p5_Arr[i] = hp->halo_SigVrad_p5[i];
		halo_SigVtheta_p3_Arr[i] = hp->halo_SigVtheta_p3[i];
		halo_SigVtheta_p4_Arr[i] = hp->halo_SigVtheta_p4[i];
		halo_SigVtheta_p5_Arr[i] = hp->halo_SigVtheta_p5[i];
		halo_SigVphi_p3_Arr[i] = hp->halo_SigVphi_p3[i];
		halo_SigVphi_p4_Arr[i] = hp->halo_SigVphi_p4[i];
		halo_SigVphi_p5_Arr[i] = hp->halo_SigVphi_p5[i];
		halo_ecc_p3_Arr[i] = hp->halo_ecc_p3[i];
		halo_ecc_p4_Arr[i] = hp->halo_ecc_p4[i];
		halo_ecc_p5_Arr[i] = hp->halo_ecc_p5[i];
		halo_Mass_p3_Arr[i] = hp->halo_Mass_p3[i];

		for(int j=0;j<hp->sat_numb;j++){	

			halo_Mass_p4_Arr[i+40*j] = hp->halo_Mass_p4[i+40*j];
			halo_Mass_p5_Arr[i+40*j] = hp->halo_Mass_p5[i+40*j];

		}

		smj_ax_p3_Arr[i] = hp->smj_ax_p3[i];
		smj_ax_p4_Arr[i] = hp->smj_ax_p4[i];
		smj_ax_p5_Arr[i] = hp->smj_ax_p5[i];

	}
		
	for(int i=0;i<50;i++){

		halo_Mass_hiRes_p3_Arr[i] = hp->halo_Mass_hiRes_p3[i];

		for(int j=0;j<hp->sat_numb;j++){	

			halo_Mass_hiRes_p4_Arr[i+50*j] = hp->halo_Mass_hiRes_p4[i+50*j];
			halo_Mass_hiRes_p5_Arr[i+50*j] = hp->halo_Mass_hiRes_p5[i+50*j];

		}

	}

	for(int i=0;i<250;i++){

		halo_Mass_hiRes2_p3_Arr[i] = hp->halo_Mass_hiRes2_p3[i];

		for(int j=0;j<hp->sat_numb;j++){	

			halo_Mass_hiRes2_p4_Arr[i+250*j] = hp->halo_Mass_hiRes2_p4[i+250*j];
			halo_Mass_hiRes2_p5_Arr[i+250*j] = hp->halo_Mass_hiRes2_p5[i+250*j];

		}

	}

	for(int i=0;i<25;i++){

		halo_Mass_hiRes0_p3_Arr[i] = hp->halo_Mass_hiRes0_p3[i];

		for(int j=0;j<hp->sat_numb;j++){	

			halo_Mass_hiRes0_p4_Arr[i+25*j] = hp->halo_Mass_hiRes0_p4[i+25*j];
			halo_Mass_hiRes0_p5_Arr[i+25*j] = hp->halo_Mass_hiRes0_p5[i+25*j];

		}

	}


	/* 
	*
	*  Create and write Header attributes	
	*
	*  First, old Gadget attributes
	*
	*/

	hid_t att_time = H5Acreate(headGroup,"Time",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_time,H5T_NATIVE_FLOAT, &hp->time);

        hid_t att_sfr = H5Acreate(headGroup,"Flag_StarFormation",H5T_NATIVE_INT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_sfr,H5T_NATIVE_INT, &hp->sfr);

        hid_t att_cooling = H5Acreate(headGroup,"Flag_Cooling",H5T_NATIVE_INT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_cooling,H5T_NATIVE_INT, &hp->cool);

        hid_t att_feed = H5Acreate(headGroup,"Flag_Feedback",H5T_NATIVE_INT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_feed,H5T_NATIVE_INT, &hp->feed);

        hid_t att_redsh = H5Acreate(headGroup,"Redshift",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_redsh,H5T_NATIVE_FLOAT, &hp->redshift);

        hid_t att_box = H5Acreate(headGroup,"BoxSize",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_box,H5T_NATIVE_FLOAT, &hp->boxSize);

        hid_t att_omega0 = H5Acreate(headGroup,"Omega0",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_omega0,H5T_NATIVE_FLOAT, &hp->omega0);

        hid_t att_omegaL = H5Acreate(headGroup,"OmegaLambda",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_omegaL,H5T_NATIVE_FLOAT, &hp->omegaL);

        hid_t att_Hubble = H5Acreate(headGroup,"HubbleParam",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_Hubble,H5T_NATIVE_FLOAT, &hp->Hubble);

        hid_t att_doublePrec = H5Acreate(headGroup,"Flag_DoublePrecision",H5T_NATIVE_INT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_doublePrec,H5T_NATIVE_INT, &hp->dPrec);

        hid_t att_nFile = H5Acreate(headGroup,"NumFilesPerSnapshot",H5T_NATIVE_INT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_nFile,H5T_NATIVE_INT, &hp->nFile);

	hid_t att_mass = H5Acreate(headGroup,"MassTable",H5T_NATIVE_FLOAT,vect_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_mass,H5T_NATIVE_FLOAT, hp->mass_table);

	hid_t att_nPart_this = H5Acreate(headGroup,"NumPart_ThisFile",H5T_NATIVE_INT,vect_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_nPart_this,H5T_NATIVE_INT, nPart_table);

        hid_t att_nPart_tot = H5Acreate(headGroup,"NumPart_Total",H5T_NATIVE_INT,vect_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_nPart_tot,H5T_NATIVE_INT, nPart_table);
	
        hid_t att_nPart_high = H5Acreate(headGroup,"NumPart_Total_HighWord",H5T_NATIVE_INT,vect_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_nPart_high,H5T_NATIVE_INT, nPart_high_table);

	hid_t att_concentr = H5Acreate(headGroup,"MW_concentration",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(att_concentr,H5T_NATIVE_FLOAT, &hp->MW_concentration);


	/* Now add new header attributes */

	std::cout << "flagging writing 1..." << std::endl;

        hid_t withBaryons = H5Acreate(headGroup,"want_Baryons",H5T_NATIVE_INT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(withBaryons,H5T_NATIVE_INT, &hp->wantBaryons);

        hid_t massBound = H5Acreate(headGroup,"sat_Mass_Bound",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(massBound,H5T_NATIVE_FLOAT, sat_massBound_Arr);

        hid_t massBound4 = H5Acreate(headGroup,"sat_Mass_Bound_p4",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(massBound4,H5T_NATIVE_FLOAT, sat_massBound_p4_Arr);

        hid_t massBound5 = H5Acreate(headGroup,"sat_Mass_Bound_p5",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(massBound5,H5T_NATIVE_FLOAT, sat_massBound_p5_Arr);

        hid_t satTotMass = H5Acreate(headGroup,"sat_Tot_Mass",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satTotMass,H5T_NATIVE_FLOAT, sat_totMass_Arr);

        hid_t satTotMass4 = H5Acreate(headGroup,"sat_Tot_Mass_p4",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satTotMass4,H5T_NATIVE_FLOAT, sat_totMass_p4_Arr);

	//std::cout << "flagging writing 2..." << std::endl;

        hid_t satTotMass5 = H5Acreate(headGroup,"sat_Tot_Mass_p5",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satTotMass5,H5T_NATIVE_FLOAT, sat_totMass_p5_Arr);

	hid_t satOutMass = H5Acreate(headGroup,"sat_Out_Mass",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satOutMass,H5T_NATIVE_FLOAT, sat_outTidMass_Arr);

	hid_t satOutMass4 = H5Acreate(headGroup,"sat_Out_Mass_p4",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satOutMass4,H5T_NATIVE_FLOAT, sat_outTidMass_p4_Arr);

	hid_t satOutMass5 = H5Acreate(headGroup,"sat_Out_Mass_p5",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satOutMass5,H5T_NATIVE_FLOAT, sat_outTidMass_p5_Arr);

	hid_t massOut242_p4 = H5Acreate(headGroup,"mass_Out_of_242_kpc_p4",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(massOut242_p4,H5T_NATIVE_FLOAT,&hp->mass_out_of_242_p4);

	hid_t massOut242_p5 = H5Acreate(headGroup,"mass_Out_of_242_kpc_p5",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(massOut242_p5,H5T_NATIVE_FLOAT,&hp->mass_out_of_242_p5);

	//std::cout << "flagging writing 3..." << std::endl;

        /* Old quantities

	hid_t satDiskMass = H5Acreate(headGroup,"sat_diskMass",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satDiskMass,H5T_NATIVE_FLOAT, &hp->sat_massDisk);

        hid_t satDiskMass4 = H5Acreate(headGroup,"sat_diskMass_p4",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satDiskMass4,H5T_NATIVE_FLOAT, &hp->sat_massDisk_p4);

        hid_t satDiskMass5 = H5Acreate(headGroup,"sat_diskMass_p5",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satDiskMass5,H5T_NATIVE_FLOAT, &hp->sat_massDisk_p5);

        hid_t satBulgeMass = H5Acreate(headGroup,"sat_bulgeMass",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satBulgeMass,H5T_NATIVE_FLOAT, &hp->sat_massBulge);

        hid_t satBulgeMass4 = H5Acreate(headGroup,"sat_bulgeMass_p4",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satBulgeMass4,H5T_NATIVE_FLOAT, &hp->sat_massBulge_p4);

        hid_t satBulgeMass5 = H5Acreate(headGroup,"sat_bulgeMass_p5",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satBulgeMass5,H5T_NATIVE_FLOAT, &hp->sat_massBulge_p5);

        //hid_t satHaloMass = H5Acreate(headGroup,"sat_haloMass",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        //H5Awrite(satHaloMass,H5T_NATIVE_FLOAT, &hp->sat_massHalo);

	*/

        hid_t satHaloMass4 = H5Acreate(headGroup,"sat_haloMass_p4",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satHaloMass4,H5T_NATIVE_FLOAT, &hp->halo_Mass_p4_tot);

        hid_t satHaloMass5 = H5Acreate(headGroup,"sat_haloMass_p5",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satHaloMass5,H5T_NATIVE_FLOAT, &hp->halo_Mass_p5_tot);

        hid_t satTidalRad = H5Acreate(headGroup,"sat_Tidal_Radius",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satTidalRad,H5T_NATIVE_FLOAT, sat_tidRadius_Arr);

        hid_t MWintMass = H5Acreate(headGroup,"MW_Integrated_Mass",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(MWintMass,H5T_NATIVE_FLOAT, MW_integratedMass_Arr);

	hid_t MWtotMass = H5Acreate(headGroup,"MW_Tot_Mass",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
	H5Awrite(MWtotMass,H5T_NATIVE_FLOAT, &hp->tot_halo_mass);

        hid_t satTidalFrac = H5Acreate(headGroup,"sat_Tidal_Mass_Fraction",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satTidalFrac,H5T_NATIVE_FLOAT, sat_tidalMass_Arr);

	std::cout << "flagging writing 4..." << std::endl;

        hid_t satTidalFrac4 = H5Acreate(headGroup,"sat_Tidal_Mass_Fraction_p4",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satTidalFrac4,H5T_NATIVE_FLOAT, sat_tidalMass_p4_Arr);

        hid_t satTidalFrac5 = H5Acreate(headGroup,"sat_Tidal_Mass_Fraction_p5",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satTidalFrac5,H5T_NATIVE_FLOAT, sat_tidalMass_p5_Arr);

        hid_t satEcc = H5Acreate(headGroup,"sat_Eccentricity",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satEcc,H5T_NATIVE_FLOAT, sat_ecc_Arr);

        hid_t satAx = H5Acreate(headGroup,"sat_sMajAxis",H5T_NATIVE_FLOAT,sat_Dspace_1D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(satAx,H5T_NATIVE_FLOAT, sat_ax_Arr);

        hid_t ItensorDark = H5Acreate(headGroup,"Inertia_Tensor_Dark",H5T_NATIVE_FLOAT,Inertia_Dspace_3D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(ItensorDark,H5T_NATIVE_FLOAT,InertiaDark_arr);

        hid_t ItensorDark_Unbound = H5Acreate(headGroup,"Inertia_Tensor_Dark_Unbound",H5T_NATIVE_FLOAT,Inertia_Dspace_3D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(ItensorDark_Unbound,H5T_NATIVE_FLOAT,InertiaDark_Unbound_arr);

        hid_t ItensorStar = H5Acreate(headGroup,"Inertia_Tensor_Star",H5T_NATIVE_FLOAT,Inertia_Dspace_3D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(ItensorStar,H5T_NATIVE_FLOAT,InertiaStar_arr);

        hid_t ItensorStar_Unbound = H5Acreate(headGroup,"Inertia_Tensor_Star_Unbound",H5T_NATIVE_FLOAT,Inertia_Dspace_3D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(ItensorStar_Unbound,H5T_NATIVE_FLOAT,InertiaStar_Unbound_arr);

	hid_t ItensorHalo = H5Acreate(headGroup,"Inertia_Tensor_Halo",H5T_NATIVE_FLOAT,Inertia_Dspace_3D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(ItensorHalo,H5T_NATIVE_FLOAT,InertiaHalo_arr);

        hid_t ItensorHalo_Unbound = H5Acreate(headGroup,"Inertia_Tensor_Halo_Unbound",H5T_NATIVE_FLOAT,Inertia_Dspace_3D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(ItensorHalo_Unbound,H5T_NATIVE_FLOAT,InertiaHalo_Unbound_arr);

        H5Awrite(ItensorDark,H5T_NATIVE_FLOAT,InertiaDark_arr);

        hid_t geomDark_Unbound = H5Acreate(headGroup,"Geom_Tensor_Dark_Unbound",H5T_NATIVE_FLOAT,Inertia_Dspace_3D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(geomDark_Unbound,H5T_NATIVE_FLOAT,geomDark_Unbound_arr);

        hid_t geomStar_Unbound = H5Acreate(headGroup,"Geom_Tensor_Star_Unbound",H5T_NATIVE_FLOAT,Inertia_Dspace_3D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(geomStar_Unbound,H5T_NATIVE_FLOAT,geomStar_Unbound_arr);

	hid_t geomHalo = H5Acreate(headGroup,"Geom_Tensor_Halo",H5T_NATIVE_FLOAT,Inertia_Dspace_3D,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(geomHalo,H5T_NATIVE_FLOAT,geomHalo_arr);

	// Disk 

	hid_t disk_Vz1 = H5Dcreate(part1,"disk_Vzed_p1",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t disk_Z1 = H5Dcreate(part1,"disk_Zed_p1",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t disk_SigVz1 = H5Dcreate(part1,"disk_SigVzed_p1",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t disk_SigZ1 = H5Dcreate(part1,"disk_SigZed_p1",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	ret = H5Dwrite(disk_Vz1, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskVzed_p1_Arr);
	ret = H5Dwrite(disk_Z1, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskZed_p1_Arr);
        ret = H5Dwrite(disk_SigVz1, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskSigVzed_p1_Arr);
        ret = H5Dwrite(disk_SigZ1, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskSigZed_p1_Arr);

	hid_t disk_vTang1 = H5Dcreate(part1,"disk_Vtan_p1",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t disk_Mass1 = H5Dcreate(part1,"disk_Mass_p1",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	//hid_t disk_zMassAll1 = H5Dcreate(part1,"disk_zMass_all_p1",H5T_NATIVE_FLOAT, short_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//hid_t disk_zMassBin1 = H5Dcreate(part1,"disk_zMass_bin_p1",H5T_NATIVE_FLOAT, multi_Dspace_2D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//hid_t disk_epicF1 = H5Dcreate(part1,"disk_Epic_freq_p1",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //hid_t disk_vertF1 = H5Dcreate(part1,"disk_Vert_freq_p1",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ret = H5Dwrite(disk_vTang1, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskVtan_p1_Arr);
        ret = H5Dwrite(disk_Mass1, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskMass_p1_Arr);
        //ret = H5Dwrite(disk_zMassAll1, H5T_NATIVE_FLOAT, short_Buffspace_1D, short_Dspace_1D, H5P_DEFAULT, hp->disk_vertical_total_p1);
        //ret = H5Dwrite(disk_zMassBin1, H5T_NATIVE_FLOAT, multi_Buffspace_2D, multi_Dspace_2D, H5P_DEFAULT, hp->disk_vertical_zone_p1);
        //ret = H5Dwrite(disk_epicF1, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->disk_vertical_zone_p1); 		//epicyc. freq	
        //ret = H5Dwrite(disk_vertF1, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->);					//vert. freq

	// Bulge 
	
	//hid_t bulge_vels2 = H5Dcreate(part2,"bulge_vels_p2",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_vRad2 = H5Dcreate(part2,"bulge_Vrad_p2",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t bulge_vTheta2 = H5Dcreate(part2,"bulge_Vtheta_p2",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_vPhi2 = H5Dcreate(part2,"bulge_Vphi_p2",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_rMass2 = H5Dcreate(part2,"bulge_rMass_p2",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t bulge_sigRad2 = H5Dcreate(part2,"bulge_SigRad_p2",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_sigTheta2 = H5Dcreate(part2,"bulge_SigTheta_p2",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t bulge_sigPhi2 = H5Dcreate(part2,"bulge_SigPhi_p2",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        //ret = H5Dwrite(bulge_vels2, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_vels_p2);

	ret = H5Dwrite(bulge_vRad2, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeVrad_p2_Arr);
        ret = H5Dwrite(bulge_vTheta2, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeVtheta_p2_Arr);
        ret = H5Dwrite(bulge_vPhi2, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeVphi_p2_Arr);
        ret = H5Dwrite(bulge_rMass2, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeMass_p2_Arr);
	ret = H5Dwrite(bulge_sigRad2, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeSigrad_p2_Arr);
        ret = H5Dwrite(bulge_sigTheta2, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeSigtheta_p2_Arr);
        ret = H5Dwrite(bulge_sigPhi2, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeSigphi_p2_Arr);

	// Halo 

	//hid_t halo_vels3 = H5Dcreate(part3,"halo_vels_p3",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vRad3 = H5Dcreate(part3,"halo_Vrad_p3",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t halo_vTheta3 = H5Dcreate(part3,"halo_Vtheta_p3",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vPhi3 = H5Dcreate(part3,"halo_Vphi_p3",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass3 = H5Dcreate(part3,"halo_rMass_p3",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass3_hRes = H5Dcreate(part3,"halo_rMass_hiRes_p3",H5T_NATIVE_FLOAT, hRes_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass3_hRes2 = H5Dcreate(part3,"halo_rMass_hiRes2_p3",H5T_NATIVE_FLOAT, hRes2_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass3_hRes0 = H5Dcreate(part3,"halo_rMass_hiRes0_p3",H5T_NATIVE_FLOAT, hRes0_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_sigRad3 = H5Dcreate(part3,"halo_SigRad_p3",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_sigTheta3 = H5Dcreate(part3,"halo_SigTheta_p3",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_sigPhi3 = H5Dcreate(part3,"halo_SigPhi_p3",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t smj_ax3 = H5Dcreate(part3,"smj_axis3",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_ecc3 = H5Dcreate(part3,"halo_eccentricity3",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        //ret = H5Dwrite(halo_vels3, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_vels_p3);
	ret = H5Dwrite(halo_vRad3, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Vrad_p3_Arr);
        ret = H5Dwrite(halo_vTheta3, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Vtheta_p3_Arr);
        ret = H5Dwrite(halo_vPhi3, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Vphi_p3_Arr);
        ret = H5Dwrite(halo_rMass3, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Mass_p3_Arr);
        ret = H5Dwrite(halo_rMass3_hRes, H5T_NATIVE_FLOAT, hRes_Buffspace_1D, hRes_Dspace_1D, H5P_DEFAULT, halo_Mass_hiRes_p3_Arr);
        ret = H5Dwrite(halo_rMass3_hRes2, H5T_NATIVE_FLOAT, hRes2_Buffspace_1D, hRes2_Dspace_1D, H5P_DEFAULT, halo_Mass_hiRes2_p3_Arr);
        ret = H5Dwrite(halo_rMass3_hRes0, H5T_NATIVE_FLOAT, hRes0_Buffspace_1D, hRes0_Dspace_1D, H5P_DEFAULT, halo_Mass_hiRes0_p3_Arr);
	ret = H5Dwrite(halo_sigRad3, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_SigVrad_p3_Arr);
        ret = H5Dwrite(halo_sigTheta3, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_SigVtheta_p3_Arr);
        ret = H5Dwrite(halo_sigPhi3, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_SigVphi_p3_Arr);
        ret = H5Dwrite(smj_ax3, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, smj_ax_p3_Arr);
        ret = H5Dwrite(halo_ecc3, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_ecc_p3_Arr);

	// Part4 

	hid_t disk_mass_p4 = H5Dcreate(part4,"disk_mass_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ret = H5Dwrite(disk_mass_p4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskMass_p4_Arr);

	hid_t disk_Vzed_p4 = H5Dcreate(part4,"disk_Vzed_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ret = H5Dwrite(disk_Vzed_p4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskVzed_p4_Arr);

	hid_t disk_Zed_p4 = H5Dcreate(part4,"disk_Zed_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ret = H5Dwrite(disk_Zed_p4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskZed_p4_Arr);

	hid_t bulge_mass_p4 = H5Dcreate(part4,"bulge_mass_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ret = H5Dwrite(bulge_mass_p4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeMass_p4_Arr);

	hid_t bulge_Vrad_p4 = H5Dcreate(part4,"bulge_Vrad_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ret = H5Dwrite(bulge_Vrad_p4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeVrad_p4_Arr);

	hid_t bulge_Vtheta_p4 = H5Dcreate(part4,"bulge_Vtheta_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ret = H5Dwrite(bulge_Vtheta_p4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeVtheta_p4_Arr);

	hid_t bulge_Vphi_p4 = H5Dcreate(part4,"bulge_Vphi_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ret = H5Dwrite(bulge_Vphi_p4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeVphi_p4_Arr);


	//hid_t halo_vels4 = H5Dcreate(part3,"halo_vels_p3",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vRad4 = H5Dcreate(part4,"halo_Vrad_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t halo_vTheta4 = H5Dcreate(part4,"halo_Vtheta_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vPhi4 = H5Dcreate(part4,"halo_Vphi_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass4 = H5Dcreate(part4,"halo_rMass_p4",H5T_NATIVE_FLOAT, sat_Dspace_large_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass4_hRes = H5Dcreate(part4,"halo_rMass_hiRes_p4",H5T_NATIVE_FLOAT, sat_Dspace_hRes_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass4_hRes2 = H5Dcreate(part4,"halo_rMass_hiRes2_p4",H5T_NATIVE_FLOAT, sat_Dspace_hRes2_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass4_hRes0 = H5Dcreate(part4,"halo_rMass_hiRes0_p4",H5T_NATIVE_FLOAT, sat_Dspace_hRes0_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_sigRad4 = H5Dcreate(part4,"halo_SigRad_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_sigTheta4 = H5Dcreate(part4,"halo_SigTheta_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_sigPhi4 = H5Dcreate(part4,"halo_SigPhi_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t smj_ax4 = H5Dcreate(part4,"smj_axis4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_ecc4 = H5Dcreate(part4,"halo_eccentricity4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        //ret = H5Dwrite(halo_vels4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_vels_p3);
	ret = H5Dwrite(halo_vRad4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Vrad_p4_Arr);
        ret = H5Dwrite(halo_vTheta4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Vtheta_p4_Arr);
        ret = H5Dwrite(halo_vPhi4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Vphi_p4_Arr);
        ret = H5Dwrite(halo_rMass4, H5T_NATIVE_FLOAT, sat_Buffspace_large_1D, sat_Dspace_large_1D, H5P_DEFAULT, halo_Mass_p4_Arr);
        ret = H5Dwrite(halo_rMass4_hRes, H5T_NATIVE_FLOAT, sat_Buffspace_hRes_1D, sat_Dspace_hRes_1D, H5P_DEFAULT, halo_Mass_hiRes_p4_Arr);
        ret = H5Dwrite(halo_rMass4_hRes2, H5T_NATIVE_FLOAT, sat_Buffspace_hRes2_1D, sat_Dspace_hRes2_1D, H5P_DEFAULT, halo_Mass_hiRes2_p4_Arr);
        ret = H5Dwrite(halo_rMass4_hRes0, H5T_NATIVE_FLOAT, sat_Buffspace_hRes0_1D, sat_Dspace_hRes0_1D, H5P_DEFAULT, halo_Mass_hiRes0_p4_Arr);
	ret = H5Dwrite(halo_sigRad4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_SigVrad_p4_Arr);
        ret = H5Dwrite(halo_sigTheta4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_SigVtheta_p4_Arr);
        ret = H5Dwrite(halo_sigPhi4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_SigVphi_p4_Arr);
        ret = H5Dwrite(smj_ax4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, smj_ax_p4_Arr);
        ret = H5Dwrite(halo_ecc4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_ecc_p4_Arr);

/*
	hid_t disk_vTang4 = H5Dcreate(part4,"disk_Vtan_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t disk_rMass4 = H5Dcreate(part4,"disk_rMass_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t disk_zMassAll4 = H5Dcreate(part4,"disk_zMass_all_p4",H5T_NATIVE_FLOAT, short_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t disk_zMassBin4 = H5Dcreate(part4,"disk_zMass_bin_p4",H5T_NATIVE_FLOAT, multi_Dspace_2D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t disk_epicF4 = H5Dcreate(part4,"disk_Epic_freq_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t disk_vertF4 = H5Dcreate(part4,"disk_Vert_freq_p4",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_vels4 = H5Dcreate(part4,"bulge_vels_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_vRad4 = H5Dcreate(part4,"bulge_Vrad_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t bulge_vTheta4 = H5Dcreate(part4,"bulge_Vtheta_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_vPhi4 = H5Dcreate(part4,"bulge_Vphi_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_rMass4 = H5Dcreate(part4,"bulge_rMass_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t bulge_sigRad4 = H5Dcreate(part4,"bulge_SigRad_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_sigTheta4 = H5Dcreate(part4,"bulge_SigTheta_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t bulge_sigPhi4 = H5Dcreate(part4,"bulge_SigPhi_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vels4 = H5Dcreate(part4,"halo_vels_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vRad4 = H5Dcreate(part4,"halo_Vrad_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t halo_vTheta4 = H5Dcreate(part4,"halo_Vtheta_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vPhi4 = H5Dcreate(part4,"halo_Vphi_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass4 = H5Dcreate(part4,"halo_rMass_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_sigRad4 = H5Dcreate(part4,"halo_SigRad_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_sigTheta4 = H5Dcreate(part4,"halo_SigTheta_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_sigPhi4 = H5Dcreate(part4,"halo_SigPhi_p4",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ret = H5Dwrite(disk_vTang4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->disk_Vtan_p4);
	ret = H5Dwrite(disk_rMass4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->disk_Mass_p4);
        ret = H5Dwrite(disk_zMassAll4, H5T_NATIVE_FLOAT, short_Buffspace_1D, short_Dspace_1D, H5P_DEFAULT, hp->disk_vertical_total_p4);
        ret = H5Dwrite(disk_zMassBin4, H5T_NATIVE_FLOAT, multi_Buffspace_2D, multi_Dspace_2D, H5P_DEFAULT, hp->disk_vertical_zone_p4);
        //ret = H5Dwrite(disk_epicF4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->disk_vertical_zone_p4); 		//epicyc. freq	
        //ret = H5Dwrite(disk_vertF4, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->);				//vert. freq
        ret = H5Dwrite(bulge_vels4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_vels_p4);
	ret = H5Dwrite(bulge_vRad4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Vrad_p4);
        ret = H5Dwrite(bulge_vTheta4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Vtheta_p4);
        ret = H5Dwrite(bulge_vPhi4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Vphi_p4);
        ret = H5Dwrite(bulge_rMass4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Mass_p4);
	ret = H5Dwrite(bulge_sigRad4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Sigrad_p4);
        ret = H5Dwrite(bulge_sigTheta4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Sigtheta_p4);
        ret = H5Dwrite(bulge_sigPhi4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Sigphi_p4);
        ret = H5Dwrite(halo_vels4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_vels_p4);
	ret = H5Dwrite(halo_vRad4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Vrad_p4);
        ret = H5Dwrite(halo_vTheta4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Vtheta_p4);
        ret = H5Dwrite(halo_vPhi4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Vphi_p4);
        ret = H5Dwrite(halo_rMass4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Mass_p4);
	ret = H5Dwrite(halo_sigRad4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Sigrad_p4);
        ret = H5Dwrite(halo_sigTheta4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Sigtheta_p4);
        ret = H5Dwrite(halo_sigPhi4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Sigphi_p4);

*/

	// Part5

	if(hp->wantBaryons == 1){

		hid_t disk_mass_p5 = H5Dcreate(part5,"disk_mass_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		ret = H5Dwrite(disk_mass_p5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskMass_p5_Arr);

		hid_t disk_Vzed_p5 = H5Dcreate(part5,"disk_Vzed_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		ret = H5Dwrite(disk_Vzed_p5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskVzed_p5_Arr);

		hid_t disk_Zed_p5 = H5Dcreate(part5,"disk_Zed_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		ret = H5Dwrite(disk_Zed_p5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, diskZed_p5_Arr);

		hid_t bulge_mass_p5 = H5Dcreate(part5,"bulge_mass_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		ret = H5Dwrite(bulge_mass_p5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeMass_p5_Arr);

		hid_t bulge_Vrad_p5 = H5Dcreate(part5,"bulge_Vrad_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		ret = H5Dwrite(bulge_Vrad_p5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeVrad_p5_Arr);

		hid_t bulge_Vtheta_p5 = H5Dcreate(part5,"bulge_Vtheta_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		ret = H5Dwrite(bulge_Vtheta_p5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeVtheta_p4_Arr);

		hid_t bulge_Vphi_p5 = H5Dcreate(part5,"bulge_Vphi_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		ret = H5Dwrite(bulge_Vphi_p5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, bulgeVphi_p5_Arr);

		//hid_t halo_vels4 = H5Dcreate(part3,"halo_vels_p3",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t halo_vRad5 = H5Dcreate(part5,"halo_Vrad_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
		hid_t halo_vTheta5 = H5Dcreate(part5,"halo_Vtheta_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t halo_vPhi5 = H5Dcreate(part5,"halo_Vphi_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t halo_rMass5 = H5Dcreate(part5,"halo_rMass_p5",H5T_NATIVE_FLOAT, sat_Dspace_large_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t halo_rMass5_hRes = H5Dcreate(part5,"halo_rMass_hiRes_p5",H5T_NATIVE_FLOAT, sat_Dspace_hRes_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t halo_rMass5_hRes2 = H5Dcreate(part5,"halo_rMass_hiRes2_p5",H5T_NATIVE_FLOAT, sat_Dspace_hRes2_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t halo_rMass5_hRes0 = H5Dcreate(part5,"halo_rMass_hiRes0_p5",H5T_NATIVE_FLOAT, sat_Dspace_hRes0_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	        hid_t halo_sigRad5 = H5Dcreate(part5,"halo_SigRad_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t halo_sigTheta5 = H5Dcreate(part5,"halo_SigTheta_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	        hid_t halo_sigPhi5 = H5Dcreate(part5,"halo_SigPhi_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t smj_ax5 = H5Dcreate(part5,"smj_axis5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	        hid_t halo_ecc5 = H5Dcreate(part5,"halo_eccentricity5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	        //ret = H5Dwrite(halo_vels4, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_vels_p3);
		ret = H5Dwrite(halo_vRad5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Vrad_p5_Arr);
	        ret = H5Dwrite(halo_vTheta5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Vtheta_p5_Arr);
	        ret = H5Dwrite(halo_vPhi5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_Vphi_p5_Arr);
   		ret = H5Dwrite(halo_rMass5, H5T_NATIVE_FLOAT, sat_Buffspace_large_1D, sat_Dspace_large_1D, H5P_DEFAULT, halo_Mass_p5_Arr);
        	ret = H5Dwrite(halo_rMass5_hRes, H5T_NATIVE_FLOAT, sat_Buffspace_hRes_1D, sat_Dspace_hRes_1D, H5P_DEFAULT, halo_Mass_hiRes_p5_Arr);
        	ret = H5Dwrite(halo_rMass5_hRes2, H5T_NATIVE_FLOAT, sat_Buffspace_hRes2_1D, sat_Dspace_hRes2_1D, H5P_DEFAULT, halo_Mass_hiRes2_p5_Arr);
	       	ret = H5Dwrite(halo_rMass5_hRes0, H5T_NATIVE_FLOAT, sat_Buffspace_hRes0_1D, sat_Dspace_hRes0_1D, H5P_DEFAULT, halo_Mass_hiRes0_p5_Arr);
		ret = H5Dwrite(halo_sigRad5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_SigVrad_p5_Arr);
	        ret = H5Dwrite(halo_sigTheta5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_SigVtheta_p5_Arr);
	        ret = H5Dwrite(halo_sigPhi5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_SigVphi_p5_Arr);
	        ret = H5Dwrite(smj_ax5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, smj_ax_p5_Arr);
	        ret = H5Dwrite(halo_ecc5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, halo_ecc_p5_Arr);

	}	

/*

	hid_t disk_vTang5 = H5Dcreate(part5,"disk_Vtan_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t disk_rMass5 = H5Dcreate(part5,"disk_rMass_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t disk_zMassAll5 = H5Dcreate(part5,"disk_zMass_all_p5",H5T_NATIVE_FLOAT, short_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t disk_zMassBin5 = H5Dcreate(part5,"disk_zMass_bin_p5",H5T_NATIVE_FLOAT, multi_Dspace_2D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//hid_t disk_epicF5 = H5Dcreate(part5,"disk_Epic_freq_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //hid_t disk_vertF5 = H5Dcreate(part5,"disk_Vert_freq_p5",H5T_NATIVE_FLOAT, large_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_vels5 = H5Dcreate(part5,"bulge_vels_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_vRad5 = H5Dcreate(part5,"bulge_Vrad_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t bulge_vTheta5 = H5Dcreate(part5,"bulge_Vtheta_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_vPhi5 = H5Dcreate(part5,"bulge_Vphi_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_rMass5 = H5Dcreate(part5,"bulge_rMass_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t bulge_sigRad5 = H5Dcreate(part5,"bulge_SigRad_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t bulge_sigTheta5 = H5Dcreate(part5,"bulge_SigTheta_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t bulge_sigPhi5 = H5Dcreate(part5,"bulge_SigPhi_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vels5 = H5Dcreate(part5,"halo_vels_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vRad5 = H5Dcreate(part5,"halo_Vrad_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);		
	hid_t halo_vTheta5 = H5Dcreate(part5,"halo_Vtheta_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_vPhi5 = H5Dcreate(part5,"halo_Vphi_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_rMass5 = H5Dcreate(part5,"halo_rMass_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_sigRad5 = H5Dcreate(part5,"halo_SigRad_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t halo_sigTheta5 = H5Dcreate(part5,"halo_SigTheta_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t halo_sigPhi5 = H5Dcreate(part5,"halo_SigPhi_p5",H5T_NATIVE_FLOAT, long_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ret = H5Dwrite(disk_vTang5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->disk_Vtan_p5);
	ret = H5Dwrite(disk_rMass5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->disk_Mass_p5);
        ret = H5Dwrite(disk_zMassAll5, H5T_NATIVE_FLOAT, short_Buffspace_1D, short_Dspace_1D, H5P_DEFAULT, hp->disk_vertical_total_p5);
        ret = H5Dwrite(disk_zMassBin5, H5T_NATIVE_FLOAT, multi_Buffspace_2D, multi_Dspace_2D, H5P_DEFAULT, hp->disk_vertical_zone_p5);
        //ret = H5Dwrite(disk_epicF5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->disk_vertical_zone_p5); 		//epicyc. freq	
        //ret = H5Dwrite(disk_vertF5, H5T_NATIVE_FLOAT, large_Buffspace_1D, large_Dspace_1D, H5P_DEFAULT, hp->);				//vert. freq
        ret = H5Dwrite(bulge_vels5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_vels_p5);
	ret = H5Dwrite(bulge_vRad5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Vrad_p5);
        ret = H5Dwrite(bulge_vTheta5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Vtheta_p5);
        ret = H5Dwrite(bulge_vPhi5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Vphi_p5);
        ret = H5Dwrite(bulge_rMass5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Mass_p5);
	ret = H5Dwrite(bulge_sigRad5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Sigrad_p5);
        ret = H5Dwrite(bulge_sigTheta5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Sigtheta_p5);
        ret = H5Dwrite(bulge_sigPhi5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->bulge_Sigphi_p5);
        ret = H5Dwrite(halo_vels5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_vels_p5);
	ret = H5Dwrite(halo_vRad5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Vrad_p5);
        ret = H5Dwrite(halo_vTheta5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Vtheta_p5);
        ret = H5Dwrite(halo_vPhi5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Vphi_p5);
        ret = H5Dwrite(halo_rMass5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Mass_p5);
	ret = H5Dwrite(halo_sigRad5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Sigrad_p5);
        ret = H5Dwrite(halo_sigTheta5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Sigtheta_p5);
        ret = H5Dwrite(halo_sigPhi5, H5T_NATIVE_FLOAT, long_Buffspace_1D, long_Dspace_1D, H5P_DEFAULT, hp->halo_Sigphi_p5);

*/

	// Density centers 

	hid_t dens_MW = H5Dcreate(densC,"dcms_MW",H5T_NATIVE_FLOAT, density_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t dens_Sat = H5Dcreate(densC,"dcms_Sat",H5T_NATIVE_FLOAT, densitySat_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//hid_t dens_Baryons = H5Dcreate(densC,"dcms_Baryons",H5T_NATIVE_FLOAT, density_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//hid_t dens_Dark = H5Dcreate(densC,"dcms_Dark",H5T_NATIVE_FLOAT, density_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t dens_MW_AngMom = H5Dcreate(densC,"dcms_AngMom",H5T_NATIVE_FLOAT, density_Dspace_1D,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		
       	ret = H5Dwrite(dens_MW, H5T_NATIVE_FLOAT, density_Buffspace_1D, density_Dspace_1D, H5P_DEFAULT, dcms_MW_arr);
	ret = H5Dwrite(dens_Sat, H5T_NATIVE_FLOAT, densitySat_Buffspace_1D, densitySat_Dspace_1D, H5P_DEFAULT, dcms_Sat_arr);
        //ret = H5Dwrite(dens_Baryons, H5T_NATIVE_FLOAT, density_Buffspace_1D, density_Dspace_1D, H5P_DEFAULT, dcms_Baryon_arr);
        //ret = H5Dwrite(dens_Dark, H5T_NATIVE_FLOAT, density_Buffspace_1D, density_Dspace_1D, H5P_DEFAULT, dcms_Dark_arr);		
        ret = H5Dwrite(dens_MW_AngMom, H5T_NATIVE_FLOAT, density_Buffspace_1D, density_Dspace_1D, H5P_DEFAULT, dcms_AngMom_arr);			

	//Close file

	H5Fclose(my_write_File);

	std::cout << "File " << hp->myHDF5Post << " written!" << std::endl;

}

//Create array file...filling will be done in the functions_make_arrays.cpp module



void create_HDF5_arrays_file(struct header_h5 *hp,struct array_hd *ap){

	std::string myStoreFile = hp->myHDF5Maps; 
	hsize_t attr_dims_vect = 6;
	hsize_t attr_dims_scal = 1;
	
	hid_t vect_Aspace = H5Screate_simple(1,&attr_dims_vect,NULL);
	hid_t scal_Aspace = H5Screate_simple(1,&attr_dims_scal,NULL);

	hid_t arrayFile = H5Fcreate((hp->myHDF5Maps).c_str(),H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t arrayHead  = H5Gcreate(arrayFile,"/Header",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);			
	hid_t arrayPart0 = H5Gcreate(arrayFile,"/PartType0",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t arrayPart1 = H5Gcreate(arrayFile,"/PartType1",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t arrayPart2 = H5Gcreate(arrayFile,"/PartType2",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t arrayPart3 = H5Gcreate(arrayFile,"/PartType3",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t arrayPart4 = H5Gcreate(arrayFile,"/PartType4",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t arrayPart5 = H5Gcreate(arrayFile,"/PartType5",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	//Fill array header
	hid_t headSoft = H5Acreate(arrayHead,"Softening",H5T_NATIVE_FLOAT,vect_Aspace,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(headSoft,H5T_NATIVE_FLOAT, &ap->res);

	hid_t headXgrid = H5Acreate(arrayHead,"x_Grids",H5T_NATIVE_INT,vect_Aspace,H5P_DEFAULT,H5P_DEFAULT);
       	H5Awrite(headXgrid,H5T_NATIVE_INT, ap->nX);

       	hid_t headYgrid = H5Acreate(arrayHead,"y_Grids",H5T_NATIVE_INT,vect_Aspace,H5P_DEFAULT,H5P_DEFAULT);
       	H5Awrite(headYgrid,H5T_NATIVE_INT, ap->nY);

       	hid_t headZgrid = H5Acreate(arrayHead,"z_Grids",H5T_NATIVE_INT,vect_Aspace,H5P_DEFAULT,H5P_DEFAULT);
       	H5Awrite(headZgrid,H5T_NATIVE_INT, ap->nZ);
    
       	hid_t headXmin = H5Acreate(arrayHead,"xMin",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
       	H5Awrite(headXmin,H5T_NATIVE_FLOAT, &ap->xMin_rel);

       	hid_t headXmax = H5Acreate(arrayHead,"xMax",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
       	H5Awrite(headXmax,H5T_NATIVE_FLOAT, &ap->xMax_rel);

       	hid_t headYMin = H5Acreate(arrayHead,"yMin",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
       	H5Awrite(headYMin,H5T_NATIVE_FLOAT, &ap->yMin_rel);

       	hid_t headYmax = H5Acreate(arrayHead,"yMax",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
       	H5Awrite(headYmax,H5T_NATIVE_FLOAT, &ap->yMax_rel);

       	hid_t headZmin = H5Acreate(arrayHead,"zMin",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
       	H5Awrite(headZmin,H5T_NATIVE_FLOAT, &ap->zMin_rel);

       	hid_t headZmax = H5Acreate(arrayHead,"zMax",H5T_NATIVE_FLOAT,scal_Aspace,H5P_DEFAULT,H5P_DEFAULT);
       	H5Awrite(headZmax,H5T_NATIVE_FLOAT, &ap->zMax_rel);

	//close file
	H5Fclose(arrayFile);		

}








