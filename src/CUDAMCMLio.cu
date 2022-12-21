/*	This file is part of CUDAMCML.

    CUDAMCML is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUDAMCML is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUDAMCML.  If not, see <http://www.gnu.org/licenses/>.*/

#define NFLOATS 5
#define NINTS 5

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

int interpret_arg(int argc, char* argv[], unsigned long long* seed)//, int* ignoreAdetection)
{

	int unknown_argument;
	for(int i=2;i<argc;i++)
	{
		unknown_argument=1;
		/*if(!strcmp(argv[i],"-A"))
		{ 
			unknown_argument=0;
			*ignoreAdetection=1;
			printf("Ignoring A-detection!\n");
		}*/
		if(!strncmp(argv[i],"-S",2) && sscanf(argv[i],"%*2c %llu",seed))
		{
		unknown_argument=0;
		printf("Seed=%llu\n",*seed);
		}
		if(unknown_argument)
		{
			printf("Unknown argument %s!\n",argv[i]);
			return 1;
		}
	}
	return 0;
}

 int Open_Output_Stream(SimulationStruct* sim)
 {
   //int k;
   
   float fl_temp[STR_LEN]; //n,g,mus;
   short int in_temp;
   //unsigned int N_phot_thd=5+ceil(sim->number_of_photons/(NUM_THREADS));
	//unsigned long int num_phot=N_phot_thd*NUM_THREADS;
	
  // k=1+ceil(sim->number_of_photons/(NUM_THREADS));
  // num_phot=k*NUM_THREADS;
  //num_phot=sim->number_of_photons+NUM_THREADS;
  
   sim->pout=fopen (sim->outp_filename , "wb+");
   if (sim->pout == NULL){perror ("Error opening output file");return 0;}
   //Write Header
   in_temp=HEADER_DIMENSION;
   fwrite(&in_temp,sizeof(short int),1,sim->pout);
   fputc('3',sim->pout);  						//File Version
   fwrite(&(sim->layers[0].n),sizeof(float),1,sim->pout);		// refr index up
   fwrite(&(sim->n_layers),sizeof(short int),1,sim->pout); 		//num_layers
   
   for (int i=0;i<sim->n_layers;i++){					//opt prop per layers
     fl_temp[4*i]=sim->layers[i+1].n;
     fl_temp[4*i+1]=1.0f/sim->layers[i+1].mutr;
     //printf("pathlength=%.6f\n",sim->layers[i+1].mutr);
     fl_temp[4*i+2]=sim->layers[i+1].g;
     fl_temp[4*i+3]=(sim->layers[i+1].z_max-sim->layers[i+1].z_min);
   }
   fwrite(fl_temp,sizeof(float),sim->n_layers*4,sim->pout);
   fwrite(&(sim->layers[sim->n_layers+1].n),sizeof(float),1,sim->pout);		// refr index down
   fputc(sim->RorT,sim->pout);						//Refl or Trans
   
   fwrite(&(sim->n_detectors),sizeof(short int),1,sim->pout);		//Num_detectors
   
   fwrite(sim->detectors,sizeof(float),2*sim->n_detectors,sim->pout); 	//detectors position
   
   //fwrite(&num_phot,sizeof(unsigned long),1,sim->pout);			//Num_ph_received per detector
   //printf("Number of photons max to receive per detectors=%lu\n",num_phot);
   //fseek(sim->pout,HEADER_DIMENSION,SEEK_SET);
   return 1;
   
   
   
   
   
  
 }
 
 void Close_Output_Strem(SimulationStruct* sim)
 {
  
   fclose(sim->pout);
 }
 
 
 
 
 int Write_Simulation_Results(float** buffer, SimulationStruct* sim)//, clock_t simulation_time)
 {
	//int i;
	unsigned long int N_photons;
	
	N_photons=2*(sim->number_of_photons/NUM_THREADS)+1;
	N_photons=N_photons*NUM_THREADS;
	//printf("k=%d\n",N_photons);
	
	fseek(sim->pout,HEADER_DIMENSION,SEEK_SET);
	  for(int j=0;j<sim->n_detectors*sim->n_layers;j++)
  	   	  fwrite(buffer[j],sizeof(float),N_photons,sim->pout);
   	   
	  
	//fseek(sim->pout,(k-1)*sizeof(float)*N_photons,SEEK_CUR);
	//}
	//fwrite(buffer[i],sizeof(float),N_photons,sim->pout);
	//fseek(sim->pout,-(int)(k*(sim->n_detectors*sim->n_layers-1)*N_photons)*sizeof(float),SEEK_CUR);
	
	
	
	
	
	//fwrite(HostMem->distance,sizeof(float),NUMSTEPS_GPU*NUM_THREADS,sim->pout);
	//fseek(sim->pout,(k-1)*sizeof(float)*NUMSTEPS_GPU*NUM_THREADS,SEEK_CUR);
	
	
 	//fclose(pFile_outp);
	
	//fwrite(HostMem->Tdistance,sizeof(float),NUMSTEPS_GPU*NUM_THREADS,pFile_outpT);
	//fwrite(HostMem->Tpath,sizeof(float),NUMSTEPS_GPU*NUM_THREADS,pFile_outpT);
 	//fclose(pFile_outpT);
 	return 0;
 
 }
int Write_Simulation_Results2(MemStruct* HostMem, SimulationStruct* sim)//, clock_t simulation_time)
 {
	//int i;
	unsigned long int N_photons_per_thd;
	
	N_photons_per_thd=2*ceil(sim->number_of_photons/NUM_THREADS)+1;
	//N_photons=N_photons*NUM_THREADS;
	//printf("k=%d\n",N_photons);
	int zeropath;
	zeropath = 0;
	
	fseek(sim->pout,HEADER_DIMENSION,SEEK_SET);
	  for(int j=0;j<sim->n_detectors*sim->n_layers;j++){
		
	  	for(int i=0;i<NUM_THREADS;i++){
	  		int n_phot=HostMem->received[i];
			
			
	  			//printf("%d\n",n_phot);
	  			fwrite(&HostMem->path[j*NUM_THREADS*N_photons_per_thd+N_photons_per_thd*i],sizeof(float),n_phot,sim->pout);
				//if(i==1){
				for(int ip=0;ip<n_phot;ip++){
				//	printf("%.10f\n",HostMem->path[j*NUM_THREADS*N_photons_per_thd+N_photons_per_thd*i+ip]);
				
				if (HostMem->path[j*NUM_THREADS*N_photons_per_thd+N_photons_per_thd*i+ip]<1e-10) zeropath++;
				}
				//}
	  	}
	  }
	printf("Zero path = %d\n",zeropath);
	  zeropath = 0;
  	  for(int j=0;j<sim->n_detectors*sim->n_layers;j++){
		for(int i=0;i<NUM_THREADS;i++){
			int n_phot=HostMem->received[i];
	  		//printf("%d\n",n_phot);
	  		fwrite(&HostMem->kappa[j*NUM_THREADS*N_photons_per_thd+N_photons_per_thd*i],sizeof(unsigned short int),n_phot,sim->pout);
			//if(i==1){
			for(int ip=0;ip<n_phot;ip++){
				//	printf("%.10f\n",HostMem->path[j*NUM_THREADS*N_photons_per_thd+N_photons_per_thd*i+ip]);
				
				if (HostMem->kappa[j*NUM_THREADS*N_photons_per_thd+N_photons_per_thd*i+ip]<1) zeropath++;
			}
				//}
	  	}
	  }
  	
  	   	  
   	return 0;
 
 }

int isnumeric(char a)
{
	if(a>=(char)48 && a<=(char)57) return 1;
	else return 0;
}

int readfloats(int n_floats, float* temp, FILE* pFile)
{
	int ii=0;
	char mystring [200];

	if(n_floats>NFLOATS) return 0; //cannot read more than NFLOATS floats

	while(ii<=0)
	{
		if(feof(pFile)) return 0; //if we reach EOF here something is wrong with the file!
		fgets(mystring , 200 , pFile);
		memset(temp,0,NFLOATS*sizeof(float));
		ii=sscanf(mystring,"%f %f %f %f %f",&temp[0],&temp[1],&temp[2],&temp[3],&temp[4]);
		if(ii>n_floats) return 0; //if we read more number than defined something is wrong with the file!
	}
	return 1; // Everyting appears to be ok!
}

int readints(int n_ints, int* temp, FILE* pFile) //replace with template?
{
	int ii=0;
	char mystring[STR_LEN];

	if(n_ints>NINTS) return 0; //cannot read more than NFLOATS floats

	while(ii<=0)
	{
		if(feof(pFile)) return 0; //if we reach EOF here something is wrong with the file!
		fgets(mystring , STR_LEN , pFile);
		memset(temp,0,NINTS*sizeof(int));
		ii=sscanf(mystring,"%d %d %d %d %d",&temp[0],&temp[1],&temp[2],&temp[3],&temp[4]);
		if(ii>n_ints) return 0; //if we read more number than defined something is wrong with the file!
		//printf("ii=%d temp=%f %f %f %f %f\n",ii,temp[0],temp[1],temp[2],temp[3],temp[4]);
	}
	return 1; // Everyting appears to be ok!
}

int ischar(char a)
{
	if((a>=(char)65 && a<=(char)90)||(a>=(char)97 && a<=(char)122)) return 1;
	else return 0;
}

int read_simulation_data(char* filename, SimulationStruct** simulations)//, int ignoreAdetection)
{
	int i=0;
	int ii=0;
	unsigned long number_of_photons;
//AF	unsigned int start_weight;
	int n_simulations = 0;
	int n_layers = 0;
	unsigned int n_detectors;
	FILE * pFile;
	char mystring [STR_LEN];
	char str[STR_LEN],buf[STR_LEN];
	float radius,dtot=0;
	
	unsigned int tmax;
	char RorT;


	float ftemp[NFLOATS];//Find a more elegant way to do this...
	int itemp[NINTS];


	pFile = fopen(filename , "r");
	if (pFile == NULL){perror ("Error opening file");return 0;}
	
	// First read the first data line (file version) and ignore
	if(!readfloats(1, ftemp, pFile)){perror ("Error reading file version");return 0;}
	//printf("File version: %f\n",ftemp[0]);

	// Second, read the number of runs
	if(!readints(1, itemp, pFile)){perror ("Error reading number of runs");return 0;}
	n_simulations = itemp[0];
	//printf("Number of runs: %d\n",n_simulations);
	
	// Allocate memory for the SimulationStruct array
	*simulations = (SimulationStruct*) malloc(sizeof(SimulationStruct)*n_simulations);
	if(*simulations == NULL){perror("Failed to malloc simulations.\n");return 0;}//{printf("Failed to malloc simulations.\n");return 0;}

	for(i=0;i<n_simulations;i++)
	{
		// Store the input filename
		strcpy((*simulations)[i].inp_filename,filename);
		// Echo the Filename
		//printf("Input filename: %s\n",filename);

		// Store ignoreAdetection data
		//AF (*simulations)[i].ignoreAdetection=ignoreAdetection;

		// Read the output filename and determine ASCII or Binary output
		ii=0;
		printf("============================ SIMULATION %d ============================\n",i+1);
		while(ii<=0)
		{
			(*simulations)[i].begin=ftell(pFile);
			fgets (mystring , STR_LEN , pFile);
			//printf("mystring=%s\n",mystring);
			ii=sscanf(mystring,"%s",str);
			//printf("str=%s\n",str);
			if(feof(pFile)|| ii>2){perror("Error reading output filename");return 0;}
			if(ii>0)ii=ischar(str[0]);
		}
		// Echo the Filename and AorB
		strcpy(buf,"./");
		strcat(buf,str);
		printf("Output filename: %s\n",buf);
		strcpy((*simulations)[i].outp_filename,buf);
		//(*simulations)[i].AorB=AorB;

		//printf("begin=%d\n",(*simulations)[i].begin);

		// Read the number of photons
		ii=0;
		
		while(ii<=0)
		{
			fgets(mystring , STR_LEN , pFile);
			number_of_photons=0;
			ii=sscanf(mystring,"%lu",&number_of_photons);
			if(feof(pFile) || ii>1){perror("Error reading number of photons");return 0;} //if we reach EOF or read more number than defined something is wrong with the file!
			//printf("ii=%d temp=%f %f %f %f %f\n",ii,temp[0],temp[1],temp[2],temp[3],temp[4]);
		}
		printf("Number of photons: %lu\n",number_of_photons);
		(*simulations)[i].number_of_photons=number_of_photons;

		//AF aggiunta TMax e Rifle or Transmittance
		// Read TMAX
		ii=0;
		 
		while(ii<=0)
		{
			fgets(mystring , STR_LEN , pFile);
			tmax=0;
			ii=sscanf(mystring,"%d",&tmax);
			if(feof(pFile) || ii>1){perror("Error reading time max");return 0;} //if we reach EOF or read more number than defined something is wrong with the file!
			(*simulations)[i].tmax=tmax;
			
			printf("Max time=%d ps\n",tmax);
		}
		// Read if calculate Reflectance or Transmittance
		ii=0;
		while(ii<=0)
		{
			fgets (mystring , STR_LEN , pFile);
			ii=sscanf(mystring,"%c",&RorT);
			if(feof(pFile)|| ii>1){perror("Error reading Reflectance or Transmittance flag");return 0;}
			if(ii>0)ii=ischar(str[0]);
			//printf("RorT=%c\n",RorT);
			(*simulations)[i].RorT = RorT;
			if((*simulations)[i].RorT =='T')
			  printf("Transmittance\n");
			else
			  printf("Reflectance\n");
		}
		
		// Read Cylinder radius
		ii=0;
		 
		while(ii<=0)
		{
			fgets(mystring , STR_LEN , pFile);
			tmax=0;
			ii=sscanf(mystring,"%f",&radius);
			if(feof(pFile) || ii>1){perror("Error reading cylinder radius");return 0;} //if we reach EOF or read more number than defined something is wrong with the file!
			(*simulations)[i].radius=radius;
			
			printf("Radius=%.3f cm\n",radius);
		}
		
		
		//Read No. of detectors
		printf("========================= DETECTORS ========================\n");
		printf("Num\tR_min (cm)\tR_max (cm)\n");
		printf("-----------------------------------\n");
		if(!readints(1, itemp, pFile)){perror ("Error reading No. of detectors");return 0;}
		//printf("No. of detectors=%d\n",itemp[0]);
		n_detectors=itemp[0];
		(*simulations)[i].n_detectors = itemp[0];
		(*simulations)[i].detectors=(float*)calloc(itemp[0]*2,sizeof(float));
		//read detectors boundaries
		
		
		
		for(ii=0;ii<n_detectors;ii++)
		{
			
		  // Read detectors data (2x float)
			if(!readfloats(2, ftemp, pFile)){perror ("Error reading detectors data");return 0;}
			(*simulations)[i].detectors[2*ii]=ftemp[0];
			(*simulations)[i].detectors[2*ii+1]=ftemp[1];
			//printf("r_min=%f, r_max=%f\n",(*simulations)[i].detectors[2*ii],(*simulations)[i].detectors[2*ii+1]);
			printf("%d\t%f\t%f\n",ii+1,(*simulations)[i].detectors[2*ii],(*simulations)[i].detectors[2*ii+1]);
		}
		printf("============================================================\n");
		
		
		
		// Read No. of layers (1xint)
		if(!readints(1, itemp, pFile)){perror ("Error reading No. of layers");return 0;}
		//printf("No. of layers=%d\n",itemp[0]);
		n_layers = itemp[0];
		(*simulations)[i].n_layers = itemp[0];

		// Allocate memory for the layers (including one for the upper and one for the lower)
		(*simulations)[i].layers = (LayerStruct*) malloc(sizeof(LayerStruct)*(n_layers+2));
		if((*simulations)[i].layers == NULL){perror("Failed to malloc layers.\n");return 0;}//{printf("Failed to malloc simulations.\n");return 0;}

		printf("===================== SAMPLE PROPERTIES ====================\n");
		printf("n\tmua (cm-1)\tmus (cm-1)\tg\tthick (cm)\n");
		printf("------------------------------------------------------------\n");

		// Read upper refractive index (1xfloat)
		if(!readfloats(1, ftemp, pFile)){perror ("Error reading upper refractive index");return 0;}
		printf("%.3f\n",ftemp[0]);
		(*simulations)[i].layers[0].n=ftemp[0];

		dtot=0;
		for(ii=1;ii<=n_layers;ii++)
		{
			
		  // Read Layer data (5x float)
			if(!readfloats(5, ftemp, pFile)){perror ("Error reading layer data");return 0;}
			printf("%.3f\t%.3f\t\t%.3f\t\t%.3f\t%.3f\n",ftemp[0],ftemp[1],ftemp[2],ftemp[3],ftemp[4]);
			(*simulations)[i].layers[ii].n=ftemp[0];
//AF			(*simulations)[i].layers[ii].mua=ftemp[1];
			(*simulations)[i].layers[ii].g=ftemp[3];
			(*simulations)[i].layers[ii].z_min=dtot;
			dtot+=ftemp[4];
			(*simulations)[i].layers[ii].z_max=dtot;
			if(ftemp[2]==0.0f) (*simulations)[i].layers[ii].mutr=FLT_MAX; //Glas layer
				
			else(*simulations)[i].layers[ii].mutr=1.0f/(ftemp[1]+ftemp[2]); //in cm
			
			  //printf("mutr=%f\n",(*simulations)[i].layers[ii].mutr);
			//printf("z_min=%f, z_max=%f\n",(*simulations)[i].layers[ii].z_min,(*simulations)[i].layers[ii].z_max);
		}//end ii<n_layers

		// Read lower refractive index (1xfloat)
		if(!readfloats(1, ftemp, pFile)){perror ("Error reading lower refractive index");return 0;}
		printf("%.3f\n",ftemp[0]);
		printf("============================================================\n\n");
		(*simulations)[i].layers[n_layers+1].n=ftemp[0];

		(*simulations)[i].end=ftell(pFile);
		//printf("end=%d\n",(*simulations)[i].end);

		//calculate start_weight
		//double n1=(*simulations)[i].layers[0].n;
		//double n2=(*simulations)[i].layers[1].n;
		//double r = (n1-n2)/(n1+n2);
		//r = r*r;
//AF		start_weight = (unsigned int)((double)0xffffffff*(1-r));
		//printf("Start weight=%u\n",start_weight);
//AF		(*simulations)[i].start_weight=start_weight;

	}//end for i<n_simulations
	return n_simulations;
}

void PrintSimulationData(SimulationStruct sim,int i){
	printf("============================ SIMULATION %d ============================\n",i);
	
	printf("Output filename: %s\n",sim.outp_filename);
	printf("Number of photons: %lu\n",sim.number_of_photons);
	printf("Max time=%d ps\n",sim.tmax);
	if(sim.RorT =='T')
			  printf("Transmittance\n");
			else
			  printf("Reflectance\n");
	
	printf("Radius=%.3f cm\n",sim.radius);
	
	//Read No. of detectors
		printf("========================= DETECTORS ========================\n");
		printf("Num\tR_min (cm)\tR_max (cm)\n");
		printf("-----------------------------------\n");
		for(int ii=0;ii<sim.n_detectors;ii++)
		{
		  printf("%d\t%f\t%f\n",ii+1,sim.detectors[2*ii],sim.detectors[2*ii+1]);
		}
		printf("============================================================\n");
		
			
		printf("===================== SAMPLE PROPERTIES ====================\n");
		printf("n\tmua (cm-1)\tmus (cm-1)\tg\tthick (cm)\n");
		printf("------------------------------------------------------------\n");
		
		printf("%.3f\n",sim.layers[0].n);
		for(int ii=1;ii<=sim.n_layers;ii++)
		{
			float thick;
			thick=sim.layers[ii].z_max-sim.layers[ii].z_min;
		  printf("%.3f\t0.000\t\t%.3f\t\t%.3f\t%.3f\n",sim.layers[ii].n,1/sim.layers[ii].mutr,sim.layers[ii].g,thick);
		}
		printf("%.3f\n",sim.layers[sim.n_layers+1].n);
		printf("============================================================\n\n");
	}

