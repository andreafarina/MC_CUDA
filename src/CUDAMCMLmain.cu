/////////////////////////////////////////////////////////////
//
//		CUDA-based Monte Carlo simulation of photon migration in layered media (CUDAMCML).
//	
//			Some documentation is avialable for CUDAMCML and should have been distrbuted along 
//			with this source code. If that is not the case: Documentation, source code and executables
//			for CUDAMCML are available for download on our webpage:
//			http://www.atomic.physics.lu.se/Biophotonics
//			or, directly
//			http://www.atomic.physics.lu.se/fileadmin/atomfysik/Biophotonics/Software/CUDAMCML.zip
//
//			We encourage the use, and modifcation of this code, and hope it will help 
//			users/programmers to utilize the power of GPGPU for their simulation needs. While we
//			don't have a scientifc publication describing this code, we would very much appreciate
//			if you cite our original GPGPU Monte Carlo letter (on which CUDAMCML is based) if you 
//			use this code or derivations thereof for your own scientifc work:
//			E. Alerstam, T. Svensson and S. Andersson-Engels, "Parallel computing with graphics processing
//			units for high-speed Monte Carlo simulations of photon migration", Journal of Biomedical Optics
//			Letters, 13(6) 060504 (2008).
//
//			To compile and run this code, please visit www.nvidia.com and download the necessary 
//			CUDA Toolkit and SKD. We also highly recommend the Visual Studio wizard 
//			(available at:http://forums.nvidia.com/index.php?showtopic=69183) 
//			if you use Visual Studio 2005 
//			(The express edition is available for free at: http://www.microsoft.com/express/2005/). 
//
//			This code is distributed under the terms of the GNU General Public Licence (see
//			below). 
//
//
///////////////////////////////////////////////////////////////

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

#include <float.h> //for FLT_MAX 
#include <stdio.h>
//#include <cutil.h>
#include "CUDAMCML.h"

//__device__ __constant__ unsigned int num_photons_dc[1];	
__device__ __constant__ unsigned int n_layers_dc[1],n_detectors_dc[1];
__device__ __constant__ unsigned long int N_photons_dc[1];
__device__ __constant__ float detectors_dc[2*MAX_DETECTORS];
//__device__ __constant__ unsigned int start_weight_dc[1];	
__device__ __constant__ LayerStruct layers_dc[MAX_LAYERS];	
//__device__ __constant__ DetStruct det_dc[1];		
__device__ __constant__ unsigned int tmax_dc[1];
__device__ __constant__ char RorT_dc;	
__device__ __constant__ float Radius_dc[1];
//__device__ unsigned long long int pos;
__device__ size_t pitch_dc[1];
//__device__ float * row;

#include "CUDAMCMLmem.cu"
#include "CUDAMCMLio.cu"
#include "CUDAMCMLrng.cu"
#include "CUDAMCMLtransport.cu"




// wrapper for device code
void DoOneSimulation(SimulationStruct* simulation,unsigned long long seed)//AF, unsigned long long* x,unsigned int* a)
{
	MemStruct DeviceMem;
	MemStruct HostMem;
	ThreadStates tstates;	
	//unsigned int threads_active_total=1;
	unsigned int size;
	unsigned int threads_active_total=1;
	
	//unsigned long N_photons;
	
	
	unsigned long long x[NUM_THREADS];//AF
	unsigned int a[NUM_THREADS];
	
	//unsigned long long int n_launched_det[MAX_DETECTORS];
	unsigned long long int num_photons_launched=0;
	unsigned long long int num_photons_received=0;
	//float **buffer;
	//unsigned long int* buffer2;

    cudaError_t cudastat;
    clock_t time1,time2;


	// Start the clock
    
	dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    dim3 dimGrid(NUM_BLOCKS);

	// x and a are already initialised in memory
	//AF HostMem.x=x;
	//HostMem.a=a;
	
	

//	if (!Open_Output_Stream(simulation))
//		return;			

	// seed =1012313ull; 
	
	 init_RNG(x, a, NUM_THREADS, "safeprimes_base32.txt", seed); //AF C'era IF ....return 1, ma la funzione è void
	
//for (i=0;i<n_cycle;i++){
	
	
	HostMem.x=&x[0];
	HostMem.a=&a[0];
	//for(int i=0;i<NUM_THREADS;i++)
	// printf("x=%u \ta=%u\n",HostMem.x[i],HostMem.a[i]);
	 
	InitMemStructs(&HostMem, &DeviceMem, &tstates, simulation);
	InitDCMem(simulation);
	
	printf("=================== RANDOM NUMBER GENERATOR =================\n");
	printf("seed\t\tx[0]\t\t\ta[0]\n");
	printf("-------------------------------------------------------------\n");
	  
	printf("%llu\t%llu\t%u\n",seed,HostMem.x[0],HostMem.a[0]);
	printf("=============================================================\n");
	
	printf("======================= RUN SIMULATIONS =====================\n");
	printf("LaunchPhotonGlobal...");
	LaunchPhoton_Global<<<dimGrid,dimBlock>>>(tstates);
	printf("Ok!\n");
	cudaDeviceSynchronize(); // Wait for all threads to finish
	cudastat=cudaGetLastError(); // Check if there was an error
	if(cudastat)printf("After LaunchPhoton Global Error code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));

	
		printf("\t\t");
		for (int j=0;j<simulation->n_detectors;j++) printf("Detector %d\t\t",j+1);
		printf("\n");
		printf("Active Threads\tReceived Ph.\t");
		for (int j=0;j<simulation->n_detectors;j++) printf("Launched Ph.\t");
		printf("\n");
		printf("-------------------------------------------------------------\n");	
		cudaDeviceSynchronize();	
		time1=clock();
		//while(threads_active_total>0){
		while(num_photons_received<simulation->number_of_photons){  
		//for  (int i=0;i<2;i++){
			//printf("DEBUG: before Mcd \n");
			MCd<<<dimGrid,dimBlock>>>(DeviceMem,tstates);
			cudaDeviceSynchronize(); // Wait for all threads to finish
			//printf("while...");
			threads_active_total=0;
			num_photons_received=0;
			num_photons_launched=0;	
			
			
			
			//copia tstates.thd_active in HostMem.thd_Active
			size=NUM_THREADS*sizeof(unsigned short int);
			cudaMemcpy(HostMem.thd_active,DeviceMem.thd_active,size,cudaMemcpyDeviceToHost);
			//printf("%u\t",HostMem.thd_active[0]);
			//somma sui threads
			for (int i=0;i<NUM_THREADS;i++) threads_active_total+=HostMem.thd_active[i];
			printf("%u\t\t",threads_active_total);
			//copia il conteggio dei fotoni ricevuti x ogni thread
			size=NUM_THREADS*sizeof(unsigned long int);
			cudaMemcpy(HostMem.received,DeviceMem.received,size,cudaMemcpyDeviceToHost);
			
			//somma sui threads
			for (int i=0;i<NUM_THREADS;i++) num_photons_received+=HostMem.received[i];
			printf("%llu\t\t",num_photons_received);
			
			//copia il conteggio dei fotoni lanciati x ogni thread
			size=NUM_THREADS*sizeof(unsigned long int);
			cudaMemcpy(HostMem.counter,DeviceMem.counter,size,cudaMemcpyDeviceToHost);
			//copia ismoving
			size=NUM_THREADS*sizeof(unsigned short int);
			cudaMemcpy(HostMem.ismoving,DeviceMem.ismoving,size,cudaMemcpyDeviceToHost);
					
			
			//somma sui threads
			for (int i=0;i<NUM_THREADS;i++) num_photons_launched+=(HostMem.counter[i]-HostMem.ismoving[i]);
			printf("%llu\n",num_photons_launched);
						
 			//for (int i=0;i<NUM_THREADS;i++){
			//printf("Launched photons in thread %d=%lu\n",i,HostMem.counter[i]);
			

		//}	
		//cudaDeviceSynchronize();
		
		}
		printf("-------------------------------------------------------------\n");
		
		cudaDeviceSynchronize(); // Wait for all threads to finish
		time2=clock();
		cudastat=cudaGetLastError(); // Check if there was an error
		if(cudastat)
			printf("Error code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		else
			printf("ExitFromDevice...Ok!\n");
		CopyDeviceToHostMem(&HostMem, &DeviceMem, simulation);
		
		printf("-------------------------------------------------------------\n");	
			
		// Print terminated photons for each detector
		printf("Det.\tLaunched ph.\tReceived ph.\n");
		printf("------------------------------------------\n");
		for (int i=0;i<simulation->n_detectors;i++){
		  	printf("%i\t%llu\t\t%llu\n",i+1,num_photons_launched,num_photons_received);
		}
		printf("------------------------------------------\n");
		
// DEBUG
// printf("Host Memory\n");
//   	for(int j=0;j<simulation->n_detectors*simulation->n_layers;j++){
//   	   for (int i=0;i<N_photons;i++) 
//   	     //{if(HostMem.path[j][i]<0.6f)
//   	     printf("%f  \t",HostMem.path[j][i]);
//   	     printf("\n");
//    	   }
// 	

	
	
	//free(buffer);
	
	
	if (!Open_Output_Stream(simulation))
		return;
	//write detected photons for each detector
	//for (int i=0;i<simulation->n_detectors;i++)
	fwrite(&num_photons_received,sizeof(unsigned long long int),simulation->n_detectors,simulation->pout);
	//for (int i=0;i<simulation->n_detectors;i++)
	fwrite(&num_photons_launched,sizeof(unsigned long long int),simulation->n_detectors,simulation->pout);
	
	
	Write_Simulation_Results2(&HostMem, simulation);
	FreeDeviceMem(&DeviceMem,&tstates);
	FreeHostMem(&HostMem);
//fseek(simulation->pout,0,SEEK_END);
//fwrite(&HostMem.ph_launched,sizeof(unsigned long long int),1,simulation->pout);
Close_Output_Strem(simulation);
//time2=clock();
printf("Simulation time: %.2f sec\n",(double)(time2-time1)/CLOCKS_PER_SEC);
printf("Simulation done!\n");

}



int main(int argc,char* argv[])
{
	int i;
	SimulationStruct* simulations;
	int n_simulations;
	unsigned long long seed = (unsigned long long) time(NULL);// Default, use time(NULL) as seed
	//int ignoreAdetection = 0;
	char* filename;
	//show cuda
	int nDevices;

  cudaGetDeviceCount(&nDevices);
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
  }
	
	if(argc<2){printf("Not enough input arguments!\n");return 1;}
	else{filename=argv[1];}

	if(interpret_arg(argc,argv,&seed)) return 1;
	printf("============================================================\n");
	printf("+                     MONTE CARLO CUDA                     +\n");
	printf("============================================================\n\n");
	n_simulations = read_simulation_data(filename, &simulations);
	
	

	if(n_simulations == 0)
	{
		printf("Something wrong with read_simulation_data!\n");
		return 1;
	}
	else
	{
		
		printf("Read %d simulations\n",n_simulations);
	}

	// Allocate memory for RNG's
	//AFunsigned long long x[NUM_THREADS];
	//unsigned int a[NUM_THREADS];

	//file = fopen("outp2.txt", "w");
	//fclose(file);
	
	//Init RNG's
	//AF if(init_RNG(x, a, NUM_THREADS, "safeprimes_base32.txt", seed)) return 1;

	
	//perform all the simulations
	for(i=0;i<n_simulations;i++)
	{
		//seed = (unsigned long long) time(NULL);// Default, use time(NULL) as seed
	 // if(init_RNG(x, a, NUM_THREADS, "safeprimes_base32.txt", seed)) return 1;
	  // Run a simulation
	 
	  PrintSimulationData(simulations[i],i+1); //Write Simulation data to display
		DoOneSimulation(&simulations[i],seed);//AF,x,a);
	printf("========================== END SIMULATION %d ==========================\n\n",i+1);
	}

	FreeSimulationStruct(simulations, n_simulations);

	return 0; 
}
