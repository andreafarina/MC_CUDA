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



float** falloc2D(int dim2, int dim1) { 

	float** m;
	float* tmp;
	int i;
	
	m = (float**) calloc(dim2,sizeof(float*));
	tmp = (float*) calloc(dim2*dim1,sizeof(float));

	for (i=0;i<dim2;i++) m[i] = tmp + i*dim1;
	//for (i=0;i<dim2;i++) m[i] = (float*)calloc(dim1,sizeof(float));
	return m;

	}

int CopyDeviceToHostMem(MemStruct* HostMem, MemStruct* DeviceMem, SimulationStruct* sim)
{ //Copy data from Device to Host memory

	cudaError_t cudastat;
	
	unsigned int size,height;
	unsigned int N_phot_thd=2*ceil(sim->number_of_photons/(NUM_THREADS))+1;
	unsigned long int N_photons=N_phot_thd*NUM_THREADS;
	
	
	//Also copy the state of the RNG's
	size=NUM_THREADS*sizeof(unsigned long long);
	printf("CopyDeviceToHostMem...\n");
	cudaMemcpy(HostMem->x,DeviceMem->x,size,cudaMemcpyDeviceToHost);
	cudastat=cudaGetLastError();
	if(cudastat){
		printf("DeviceMem->x to HostMem->x copy error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		exit(1);
	}
	
	
	size=N_photons*sizeof(float);
	height=sim->n_detectors*sim->n_layers;
	
	cudaMemcpy(HostMem->path,DeviceMem->path,size*height,cudaMemcpyDeviceToHost);
	size=N_photons*sizeof(unsigned short int);
	cudaMemcpy(HostMem->kappa,DeviceMem->kappa,size*height,cudaMemcpyDeviceToHost);
	size=N_photons*sizeof(float);
	cudaMemcpy(HostMem->zmax,DeviceMem->zmax,size*sim->n_detectors,cudaMemcpyDeviceToHost);
	cudaMemcpy(HostMem->sumz,DeviceMem->sumz,size*sim->n_detectors,cudaMemcpyDeviceToHost);
		
	//cudaDeviceSynchronize();	
	cudastat=cudaGetLastError();
	if(cudastat){
		printf("DeviceMem->path to HostMem->path or kappa copy error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		exit(1);
	}
	printf("Ok!\n");
	return 0;
}


// Initialize the Device Constant Memory: 0->OK, 1->Error
int InitDCMem(SimulationStruct* sim)
{
	printf("InitConstMem...");
	cudaError_t cudastat;
	unsigned int N_phot_thd=2*ceil(sim->number_of_photons/(NUM_THREADS))+1;
	
	// Copy num_photons_dc to constant device memory
	 cudaMemcpyToSymbol(n_layers_dc,&(sim->n_layers),sizeof(short int));

	// Copy layer data to constant device memory
	cudaMemcpyToSymbol(layers_dc,sim->layers,(sim->n_layers+2)*sizeof(LayerStruct));
	
	// Copy num detectors to constant memory
	cudaMemcpyToSymbol(n_detectors_dc,&(sim->n_detectors),sizeof(short int));
	
	// Copy detector coordinates to constant memory
	cudaMemcpyToSymbol(detectors_dc,sim->detectors,(sim->n_detectors)*2*sizeof(float));
	
	// Copy max time to constant memory
	cudaMemcpyToSymbol(tmax_dc,&(sim->tmax),sizeof(unsigned int));
	
	// Copy RorT (reflectance or transmittance) flag to constant memory
	cudaMemcpyToSymbol(RorT_dc,&(sim->RorT),sizeof(char));
	
	// Copy Cylinder radius to constant memory
	cudaMemcpyToSymbol(Radius_dc,&(sim->radius),sizeof(float));
	
	// Copy num photons per thread to receive in constant memory
	cudaMemcpyToSymbol(N_photons_dc,&(N_phot_thd),sizeof(unsigned int));

	cudastat=cudaGetLastError();
	if(cudastat){
		printf("Device Constant Memory initialization error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		exit(1);
	}
	printf("Ok!\n");
	return 0;
}

// Initialize Memory: 0->OK, 1->Error
int InitMemStructs(MemStruct* HostMem, MemStruct* DeviceMem, ThreadStates* tstates, SimulationStruct* sim)
{
	
	cudaError_t cudastat;
	unsigned int size,height;
	unsigned int N_phot_thd=2*ceil(sim->number_of_photons/(NUM_THREADS))+1;
	unsigned long int N_photons=N_phot_thd*NUM_THREADS;//sim->number_of_photons+NUM_THREADS;
	
	printf("=========================== MEMORY =========================\n");
	printf("Num Blocks\tNum Threads per block\tNum threads\n");	
	printf("------------------------------------------------------------\n");
	printf("%d\t\t%d\t\t\t%d\n",NUM_BLOCKS,NUM_THREADS_PER_BLOCK,NUM_THREADS);
	printf("------------------------------------------------------------\n");
	
	printf("Photons per thread=%i\n",N_phot_thd);
	printf("Effective number of photons=%lu\n",N_photons);
	printf("------------------------------------------------------------\n");
	
	printf("InitMem...");
	// Allocate x and a on the device (For MWC RNG)
    	size=NUM_THREADS*sizeof(unsigned long long);
    	cudaMalloc((void**)&DeviceMem->x,size);
    	cudaMemcpy(DeviceMem->x,HostMem->x,size,cudaMemcpyHostToDevice);
	
	size=NUM_THREADS*sizeof(unsigned int);
    	cudaMalloc((void**)&DeviceMem->a,size);
    	cudaMemcpy(DeviceMem->a,HostMem->a,size,cudaMemcpyHostToDevice);
	cudastat=cudaGetLastError();
	if(cudastat){
		printf("DeviceMem->x or DeviceMem->a allocation error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		exit(1);
	}
		
	// Allocate 2D matrix in Device & Host Memory for pathlength: row->photons; height: rec1*n_layers + rec2*n_layers
	HostMem->path=(float*) calloc(sim->n_detectors*sim->n_layers*N_photons,sizeof(float));
	HostMem->kappa=(unsigned short int*) calloc(sim->n_detectors*sim->n_layers*N_photons,sizeof(unsigned short int));
	HostMem->zmax=(float*) calloc(sim->n_detectors*N_photons,sizeof(float));
	HostMem->sumz=(float*) calloc(sim->n_detectors*N_photons,sizeof(float));
	if(HostMem->path==NULL){printf("Error allocating HostMem->path"); exit (1);}
	if(HostMem->kappa==NULL){printf("Error allocating HostMem->kappa"); exit (1);}
	if(HostMem->zmax==NULL){printf("Error allocating HostMem->zmax"); exit (1);}
	if(HostMem->sumz==NULL){printf("Error allocating HostMem->sumz"); exit (1);}
	
	size=N_photons*sizeof(float);	//width
	height=sim->n_detectors*sim->n_layers;
	
	// 2D linearized vector
	cudaMalloc((void**)&DeviceMem->path,size*height);
	cudaMemset(DeviceMem->path,0.0,size*height);
	
	size=N_photons*sizeof(unsigned short int);	//width
	cudaMalloc((void**)&DeviceMem->kappa,size*height);
	cudaMemset(DeviceMem->kappa,0.0,size*height);
	
	size=N_photons*sizeof(float);	//width
	cudaMalloc((void**)&DeviceMem->zmax,size*sim->n_detectors);
	cudaMemset(DeviceMem->zmax,0.0,size*sim->n_detectors);
	
	cudaMalloc((void**)&DeviceMem->sumz,size*sim->n_detectors);
	cudaMemset(DeviceMem->sumz,0.0,size*sim->n_detectors);
	
	cudastat=cudaGetLastError();
	if(cudastat){
		printf("DeviceMem->path or kappa or zmax or sumz allocation error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		exit(1);
	}
	
	
	
	// Allocate counter of photon launched on Device Memory
	size=sizeof(unsigned long int);
	HostMem->counter=(unsigned long int *)calloc(NUM_THREADS,size);
	if(HostMem->counter==NULL){printf("Error allocating HostMem->counter"); exit (1);}
	
	cudaMalloc((void**)&DeviceMem->counter,NUM_THREADS*size); //conta i fotoni lanciati per ogni thread per evitare di usare una AtomicAdd
	cudaMemset(DeviceMem->counter,0,NUM_THREADS*size);
	cudastat=cudaGetLastError();
	if(cudastat){
		printf("DeviceMem->counter allocation error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		exit(1);
	}
	
	size=sizeof(unsigned long int);
	HostMem->received=( int *)calloc(NUM_THREADS,size);
	if(HostMem->received==NULL){printf("Error allocating HostMem->received"); exit (1);}
	
	cudaMalloc((void**)&DeviceMem->received,NUM_THREADS*size); //conta i fotoni lanciati per ogni thread per evitare di usare una AtomicAdd
	cudaMemset(DeviceMem->received,0,NUM_THREADS*size);
	cudastat=cudaGetLastError();
	if(cudastat){
		printf("DeviceMem->received allocation error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		exit(1);
	}
	
	// NEW Allocate flag for thd_active  on Device Memory
	size=sizeof(unsigned short int);
	HostMem->thd_active=(unsigned short int *)calloc(NUM_THREADS,size);
	if(HostMem->thd_active==NULL){printf("Error allocating HostMem->thd_active"); exit (1);}
	
	for (int i=0;i<NUM_THREADS;i++) HostMem->thd_active[i]=1;

	// inizialize al valore 1
	cudaMalloc((void**)&DeviceMem->thd_active,NUM_THREADS*size); 
	cudaMemcpy(DeviceMem->thd_active,HostMem->thd_active,NUM_THREADS*size,cudaMemcpyHostToDevice);
	cudastat=cudaGetLastError();
	if(cudastat){
		printf("DeviceMem->thd_active allocation error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		exit(1);
	}
	
	// Allocate Host memory for ismoving
	size=sizeof(unsigned short int);
 	HostMem->ismoving=(unsigned short int *)calloc(NUM_THREADS,size);
 	if(HostMem->ismoving==NULL){printf("Error allocating HostMem->ismoving"); exit (1);}
 	
	 cudaMalloc((void**)&DeviceMem->ismoving,NUM_THREADS*size); 
 	cudaMemset(DeviceMem->ismoving,0,NUM_THREADS*size);
 	cudastat=cudaGetLastError();
 	if(cudastat){
 		printf("DeviceMem->ismoving allocation error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
 		exit(1);
 	}
	
	//Allocate tstates photon 
	size=NUM_THREADS*sizeof(float);
	cudaMalloc((void**)&tstates->x,size);
	cudaMalloc((void**)&tstates->y,size); 
	cudaMalloc((void**)&tstates->z,size);
	cudaMalloc((void**)&tstates->dx,size);
	cudaMalloc((void**)&tstates->dy,size);
	cudaMalloc((void**)&tstates->dz,size);
	cudaMalloc((void**)&tstates->t,sim->n_layers*size);
	cudaMalloc((void**)&tstates->time_tot,size);
	cudaMalloc((void**)&tstates->t_left,size);
	cudaMalloc((void**)&tstates->zmax,size);
	cudaMalloc((void**)&tstates->sumz,size);
	size=NUM_THREADS*sizeof(unsigned short int);
	cudaMalloc((void**)&tstates->k,sim->n_layers*size);
	cudaMalloc((void**)&tstates->layer,size);
	cudaMalloc((void**)&tstates->weight,size);
	
	cudastat=cudaGetLastError();
	if(cudastat){
		printf("ThreadStates->p allocation error, code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));
		exit(1);
	}	
	
	printf("Ok!\n");
	    
	return 0;
}

void FreeDeviceMem(MemStruct* DeviceMem, ThreadStates* tstates)
{
	cudaFree(DeviceMem->x); DeviceMem->x = NULL;
  	cudaFree(DeviceMem->a); DeviceMem->a = NULL;
  	cudaFree(DeviceMem->path);
  	cudaFree(DeviceMem->kappa);
  	cudaFree(tstates->x); tstates->x = NULL;
  	cudaFree(tstates->y); tstates->y = NULL;
  	cudaFree(tstates->z); tstates->z = NULL;
  	cudaFree(tstates->dx); tstates->dx = NULL;
  	cudaFree(tstates->dy); tstates->dy = NULL;
  	cudaFree(tstates->dz); tstates->dz = NULL;
  	cudaFree(tstates->weight); tstates->weight = NULL;
  	cudaFree(tstates->layer); tstates->layer = NULL;
  	cudaFree(tstates->time_tot); tstates->time_tot = NULL;
  	cudaFree(tstates->t); tstates->t = NULL;
  	cudaFree(tstates->k); tstates->k = NULL;
  	cudaFree(tstates->t_left); tstates->t_left = NULL;
  	cudaFree(tstates->zmax); tstates->zmax = NULL;
  	cudaFree(tstates->sumz); tstates->sumz = NULL;
}

void FreeHostMem(MemStruct* HostMem){  
    	free(HostMem->path);
    	free(HostMem->kappa);
    	free(HostMem->zmax);
    	free(HostMem->sumz);
   	free(HostMem->counter);
   	free(HostMem->received);
    	free(HostMem->thd_active);
    	free(HostMem->ismoving);
}   

void FreeSimulationStruct(SimulationStruct* sim, int n_simulations)
{
	for(int i=0;i<n_simulations;i++)free(sim[i].layers);
	free(sim);
}

