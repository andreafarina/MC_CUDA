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

// DEFINES

//The register usage varies with platform. 64-bit Linux and 32.bit Windows XP have been tested.

//#ifdef __linux__ //uses 25 registers per thread (64-bit)
//	#define NUM_THREADS_PER_BLOCK 16//320 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
//	#define NUM_THREADS 256//20480//MAX 61440
//#endif

#ifdef __linux__ //uses 25 registers per thread (64-bit)
	#define NUM_BLOCKS 2048         //256//256// 330 //64//MAX 192//64//256  //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
	#define NUM_THREADS_PER_BLOCK 256   //256 //224//384//176//384//320 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
	#define NUM_THREADS  (NUM_BLOCKS*NUM_THREADS_PER_BLOCK) //81920//98304//49152//16384//5632//13440//24576//11264//22528//11264//24576//MAX 61440
#endif
//DEBUG
// #ifdef __linux__ //uses 25 registers per thread (64-bit)
// 	#define NUM_BLOCKS  1//64//MAX 192//64//256  //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
// 	#define NUM_THREADS_PER_BLOCK 1//384//176//384//320 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
// 	#define NUM_THREADS 1//98304//49152//16384//5632//13440//24576//11264//22528//11264//24576//MAX 61440
// #endif
//FINE DEBUG

#ifdef _WIN32 //uses 26 registers per thread
	#define NUM_BLOCKS  256//32//64//MAX 192//64//256  //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
	#define NUM_THREADS_PER_BLOCK 48 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
	#define NUM_THREADS 12288
#endif
#ifdef __APPLE__ //uses 25 registers per thread (64-bit)
#define NUM_BLOCKS  128//32//64//MAX 192//64//256  //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
#define NUM_THREADS_PER_BLOCK 128//384//176//384//320 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
#define NUM_THREADS 16384//5632//13440//24576//11264//22528//11264//24576//MAX 61440
#endif



#define NUMSTEPS_GPU 10000 //MAX 1000
#define PI 3.141592654f
//#define RPI 0.318309886f
#define MAX_LAYERS 20
#define MAX_DETECTORS 1
#define STR_LEN 200
#define HEADER_DIMENSION 512 //Byte
#define LIGHT_SPEED 0.0299792458f //0.03f //cm/ps
// Cilindro e fibre
#define NA 0.38f		//fiber_NA in/out
#define RAD_FIBER 0.05f		//input fiber_NA
//#define N_LATERAL 1.0f	//refr index laterale
//#define RADIUS 0.595f		//cylinder radius

//#define WEIGHT 0.0001f
//#define WEIGHTI 429497u //0xFFFFFFFFu*WEIGHT
//#define CHANCE 0.1f


// TYPEDEFS
typedef struct //__align__(16)
{
	float z_min;		// Layer z_min [cm]
	float z_max;		// Layer z_max [cm]
	float mutr;		// Reciprocal mu_total [cm]
	float g;		// Anisotropy factor [-]
	float n;		// Refractive index [-]
}LayerStruct;

typedef struct //__align__(16) 
{
	float x;		// Global x coordinate [cm]
	float y;		// Global y coordinate [cm]
	float z;		// Global z coordinate [cm]
	float dx;		// (Global, normalized) x-direction
	float dy;		// (Global, normalized) y-direction
	float dz;		// (Global, normalized) z-direction
	float t[MAX_LAYERS];	// time of flight for each layer
	float time_tot;		// total time of flight
	float t_left;	//residual step
	unsigned short int weight;	// Photon weight
	unsigned short int layer;		// Current layer
}PhotonStruct;

typedef struct //__align__(16) 
{
	float* x;		// Global x coordinate [cm]
	float* y;		// Global y coordinate [cm]
	float* z;		// Global z coordinate [cm]
	float* dx;		// (Global, normalized) x-direction
	float* dy;		// (Global, normalized) y-direction
	float* dz;		// (Global, normalized) z-direction
	float* t;//t[MAX_LAYERS];	// time of flight for each layer coalesced t1:NUM_THD:t2:NUM_THD:t3...
	float* time_tot;		// total time of flight
	float* t_left;
	unsigned short int* weight;	// Photon weight
	unsigned short int* layer;		// Current layer
	//unsigned long int* counter; 	// Pointer all'array contenente i fotoni lanciati in ogni thread
	//unsigned long int* received;    //d1:NUM_THD:t2:NUM_THD.....
	//unsigned short int* thd_active;  // Thd is active or all detector at saturated 
}ThreadStates;
typedef struct 
{
	unsigned long int number_of_photons;
	//int ignoreAdetection;
	short int n_layers;
	short int n_detectors;
	float* detectors; 	//rho_first and rho_last for each detector
	//unsigned int time_max;
	//unsigned int start_weight;
	char outp_filename[STR_LEN];
	FILE* pout; //pointer to the output file
	char inp_filename[STR_LEN];
	long begin,end;
	//char AorB;
	float radius;
	unsigned int tmax; //TMAX ps
	char RorT; //1=Reflectance or 0=transmittance 
	//DetStruct det;
	LayerStruct* layers;
	size_t pitch;
}SimulationStruct;


typedef struct //__align__(8)
{
	//PhotonStruct* p;			// Pointer to structure array containing all the photon data
	unsigned long long *x;			// Pointer to the array containing all the WMC x's
	unsigned int *a;			// Pointer to the array containing all the WMC a's
	
	unsigned short int* ismoving;	// Controlla quanti fotoni sono in movimento x poi sottrarli, alla fine della sim, a quelli lanciati
	int* received; 	//num_detectors array containing  detected photons for each thread
	unsigned long int* counter;		// Pointer all'array contenente i fotoni lanciati in ogni thread
	//unsigned short int* sat_det;
	unsigned short int* thd_active;
	
	float* path;				// path for each layers and detecors
	//size_t pitch;				// pitch for 2D matrix
	//unsigned long long int *ph_launched;	// launched photons
}MemStruct;

