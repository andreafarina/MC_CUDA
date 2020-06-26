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

// forward declaration of the device code
//AF template <int ignoreAdetection> __global__ void MCd(MemStruct);
__global__ void MCd(MemStruct DeviceMem,ThreadStates tstates);
__device__ float rand_MWC_oc(unsigned long long*,unsigned int*);
__device__ float rand_MWC_co(unsigned long long*,unsigned int*);

__device__ void LaunchPhoton(PhotonStruct*);//,unsigned long long*, unsigned int*);
__global__ void LaunchPhoton_Global(ThreadStates);
__device__ void RestoreStates(MemStruct* DeviceMem,ThreadStates* tstates,PhotonStruct* p,unsigned long long* x,unsigned int* a);
__device__ void SaveStates(MemStruct* DeviceMem,ThreadStates* tstates,PhotonStruct* p,unsigned long long x);
__device__ void ComputeStepSize(PhotonStruct* p, float* s,unsigned long long* x, unsigned int* a);
__device__ void Hop(PhotonStruct* p,float s);
__device__ void Spin(PhotonStruct*, float,unsigned long long*,unsigned int*);
__device__ int Intersection(PhotonStruct* p, float* s);
__device__ unsigned int Reflect(PhotonStruct*, int, unsigned long long*, unsigned int*);
__device__ void DetectAndSavePath(MemStruct* DeviceMem,PhotonStruct* p, int new_layer,int thd_id);

__device__ void AtomicAddULL(unsigned long long* address, unsigned int add);
__device__ void AtomicAddUL(unsigned long* address, unsigned int add);

//********************************************************************************

//*********************************** MCd ****************************************
__global__ void MCd(MemStruct DeviceMem, ThreadStates tstates)
 {
	int thd_id=NUM_THREADS_PER_BLOCK*blockIdx.x+threadIdx.x;
    	unsigned int ii=0,iexit=NUMSTEPS_GPU; 
 	unsigned long long int x;	//coherent
 	unsigned int a;			//coherent
 	PhotonStruct p;
 	
 	float s;			//step length
 	short int new_layer;
	
 	RestoreStates(&DeviceMem,&tstates,&p,&x,&a);
 	
	if (!DeviceMem.thd_active[thd_id]) iexit=0;
	for(ii = 0;ii < iexit;ii++) //this is the main while loop
 	{
 		ComputeStepSize(&p, &s, &x, &a);
 		new_layer = Intersection(&p, &s);
 		Hop(&p, s);
 		if (p.time_tot >= tmax_dc[0]) p.weight=0u;
 			
 		if (p.weight == 1u){
			if(new_layer != p.layer){
				if(Reflect(&p, new_layer, &x, &a) == 0u){	//Check for reflection
					DetectAndSavePath(&DeviceMem, &p, new_layer, thd_id);
				} 
			} else Spin(&p, layers_dc[p.layer].g, &x, &a);
		} 
		
 		if ((p.weight == 0u)){ //NEW PHOTON
 			if (DeviceMem.received[thd_id] < N_photons_dc[0]){  
 		    		LaunchPhoton(&p);
 		  		DeviceMem.counter[thd_id]++;
 		  	}else{
				DeviceMem.thd_active[thd_id]=0;
				ii=NUMSTEPS_GPU;
			}
		 }
 		
 	}//end for   
 	SaveStates(&DeviceMem,&tstates,&p,x);
}//end MCd

//*********************************** LaunchPhoton ****************************************
__device__ void LaunchPhoton(PhotonStruct* p)//per fibra con NA, unsigned long long* x, unsigned int* a)
{
	// We are currently not using the RNG but might do later
	p->x  = 0.0f;
	p->y  = 0.0f;
	p->z  = 0.0f;
	p->dx = 0.0f;
	p->dy = 0.0f;
	p->dz = 1.0f;
		
	for (short int i = 0;i < n_layers_dc[0];i++) p->t[i] = 0.0f;
	p->time_tot=0.0f;
	p->layer = 1;
	p->weight=1u;
	p->t_left=0.0f;
	//AF	p->weight = *start_weight_dc; //specular reflection!

}

//*********************************** LaunchPhoton_Global ****************************************
__global__ void LaunchPhoton_Global(ThreadStates tstates)//PhotonStruct* pd, unsigned long long* x, unsigned int* a)
{
	int thd_id=NUM_THREADS_PER_BLOCK*blockIdx.x+threadIdx.x;
    	PhotonStruct p;

	LaunchPhoton(&p);
	
	tstates.x[thd_id]=p.x;
	tstates.y[thd_id]=p.y;
	tstates.z[thd_id]=p.z;
	tstates.dx[thd_id]=p.dx;
	tstates.dy[thd_id]=p.dy;
	tstates.dz[thd_id]=p.dz;
	for (short int i=0;i<n_layers_dc[0];i++) tstates.t[i*NUM_THREADS+thd_id]=p.t[i];
	tstates.time_tot[thd_id]=p.time_tot;
	tstates.layer[thd_id]=p.layer;
	tstates.weight[thd_id]=p.weight;
	tstates.t_left[thd_id]=p.t_left;	
}

//************************************* RestoreStates ******************************************
__device__ void RestoreStates(MemStruct* DeviceMem,ThreadStates* tstates,PhotonStruct* p,unsigned long long* x, unsigned int* a)
{
	int thd_id=NUM_THREADS_PER_BLOCK*blockIdx.x+threadIdx.x;
	*x = DeviceMem->x[thd_id];//coherent
 	*a = DeviceMem->a[thd_id];//coherent
 	p->x = tstates->x[thd_id];
 	p->y = tstates->y[thd_id];
 	p->z = tstates->z[thd_id];
 	p->dx = tstates->dx[thd_id];
 	p->dy = tstates->dy[thd_id];
 	p->dz = tstates->dz[thd_id];
 	for (short int i = 0;i < n_layers_dc[0];i++) p->t[i] = tstates->t[i*NUM_THREADS+thd_id];
 	p->time_tot = tstates->time_tot[thd_id];
 	p->weight = tstates->weight[thd_id];
 	p->layer = tstates->layer[thd_id];
 	p->t_left = tstates->t_left[thd_id];
}	

//************************************* SaveStates ******************************************
__device__ void SaveStates(MemStruct* DeviceMem,ThreadStates* tstates,PhotonStruct* p,unsigned long long x)
{
	int thd_id=NUM_THREADS_PER_BLOCK*blockIdx.x+threadIdx.x;
	//save the state of the MC simulation in global memory before exiting
 	//DeviceMem->p[thd_id] = p;	//This one is incoherent!!!
 	DeviceMem->x[thd_id] = x;  //this one also seems to be coherent
 	DeviceMem->ismoving[thd_id] = p->weight;
 	tstates->x[thd_id] = p->x;
	tstates->y[thd_id] = p->y;
	tstates->z[thd_id] = p->z;
	tstates->dx[thd_id] = p->dx;
	tstates->dy[thd_id] = p->dy;
	tstates->dz[thd_id] = p->dz;
	for (short int i = 0;i < n_layers_dc[0];i++) tstates->t[i*NUM_THREADS+thd_id] = p->t[i];
	tstates->time_tot[thd_id] = p->time_tot;
	tstates->layer[thd_id] = p->layer;
	tstates->weight[thd_id] = p->weight;
	tstates->t_left[thd_id] = p->t_left;
}	

//************************************* ComputeStepSize ******************************************
__device__ void ComputeStepSize(PhotonStruct* p, float* s,unsigned long long* x, unsigned int* a)
{
	if(layers_dc[p->layer].mutr != FLT_MAX){
 			if (p->t_left == 0.0f){
 				float rand = rand_MWC_oc(x,a);
 				*s = -__logf(rand)*layers_dc[p->layer].mutr;//sample step length [cm] //HERE AN OPEN_OPEN FUNCTION WOULD BE APPRECIATED
 			}else{
 				*s = p->t_left*layers_dc[p->layer].mutr;
 				p->t_left = 0.0f;
 			}
 	}else{
 		*s=100.0f;
 	}
 }

//************************************* Hop ******************************************
__device__ void Hop(PhotonStruct *p,float s)
{
  	p->x += s * p->dx;
  	p->y += s * p->dy;
  	p->z += s * p->dz;

  	p->t[p->layer-1] += s;
  	p->time_tot += s/LIGHT_SPEED*layers_dc[p->layer].n;
}

//************************************* Spin ******************************************
__device__ void Spin(PhotonStruct* p, float g, unsigned long long* x, unsigned int* a)
{
	float cost, sint;	// cosine and sine of the polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the azimuthal angle psi. 
	float temp;
	float tempdir=p->dx;
		
	//GLASS LAYER: se si muove parallelo ai layer non scatterare
	
	//if (layers_dc[p->layer].mutr==FLT_MAX){
	//	printf("glass\n");
	//	 return; //we are in glass...so no spin
	//}
	
	//This is more efficient for g!=0 but of course less efficient for g==0
	temp = __fdividef((1.0f-(g)*(g)),(1.0f-(g)+2.0f*(g)*rand_MWC_co(x,a)));//Should be close close????!!!!!
	cost = __fdividef((1.0f+(g)*(g) - temp*temp),(2.0f*(g)));
	if(g==0.0f)
		cost = 2.0f*rand_MWC_co(x,a) -1.0f;//Should be close close??!!!!!

	sint = sqrtf(1.0f - cost*cost);

	__sincosf(2.0f*PI*rand_MWC_co(x,a),&sinp,&cosp);// spin psi [0-2*PI)
	
	temp = sqrtf(1.0f - p->dz*p->dz);

	if(temp==0.0f) //normal incident.
	{
		p->dx = sint*cosp;
		p->dy = sint*sinp;
		p->dz = copysignf(cost,p->dz*cost);
	}
	else // regular incident.
	{
		p->dx = __fdividef(sint*(p->dx*p->dz*cosp - p->dy*sinp),temp) + p->dx*cost;
		p->dy = __fdividef(sint*(p->dy*p->dz*cosp + tempdir*sinp),temp) + p->dy*cost;
		p->dz = -sint*cosp*temp + p->dz*cost;
	}

	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp=rsqrtf(p->dx*p->dx+p->dy*p->dy+p->dz*p->dz);
	p->dx = p->dx*temp;
	p->dy = p->dy*temp;
	p->dz = p->dz*temp;
}// end Spin

//************************************* Intersection ******************************************
// 29-11-2011 MODIFICHE CYLINDER
__device__ int Intersection(PhotonStruct* p, float* s)
{
	float A,B,C,ss=*s;
	short int new_layer=p->layer;
		 
 	if((p->x+(*s)*p->dx)*(p->x+(*s)*p->dx)+(p->y+(*s)*p->dy)*(p->y+(*s)*p->dy)>Radius_dc[0]*Radius_dc[0]) //Check for lateral boundary hit
 	{
		A=1.0f-(p->dz)*(p->dz);
	  	B=(p->dx)*(p->x)+(p->dy)*(p->y);
	  	C=(p->x)*(p->x)+(p->y)*(p->y)-Radius_dc[0]*Radius_dc[0];
	  	(*s)=fdividef(-B+sqrtf(B*B-A*C),A);
	  	//if (s_cyl<0) s_cyl=(-B-sqrtf(B*B-A*C))/A;
	  	p->t_left=(ss-(*s))/layers_dc[p->layer].mutr;
	  	new_layer=-10;
	}
	if(p->z+(*s)*p->dz < layers_dc[p->layer].z_min){
		*s = __fdividef(layers_dc[p->layer].z_min-p->z,p->dz);
		//printf("down\n");
		p->t_left=(ss-(*s))/layers_dc[p->layer].mutr;
		return p->layer-1;
	} //Check for upwards reflection/transmission & calculate new s
 	if(p->z+(*s)*p->dz > layers_dc[p->layer].z_max){
 		*s = __fdividef(layers_dc[p->layer].z_max-p->z,p->dz);
 		//printf("up\n");
 		p->t_left=(ss-(*s))/layers_dc[p->layer].mutr;
 		return p->layer+1;
 	} //Check for downward reflection/transmission
	 
	  //QUALE soluz prendere??
// 	  if (s_cyl<(*s))
// 	  {
// 	    *s=s_cyl;//fdividef(-B+sqrtf(B*B-A*C),A);
// 	    new_layer=-10;
// 	  }
	return new_layer;
}

//************************************* Reflect ******************************************
__device__ unsigned int Reflect(PhotonStruct* p, int new_layer, unsigned long long* x, unsigned int* a)
{
	//Calculates whether the photon is reflected (returns 1) or not (returns 0)
	// Reflect() will also update the current photon layer (after transmission) and photon direction (both transmission and reflection)


	float n1 = layers_dc[p->layer].n;
	float rad=sqrtf(p->x*p->x+p->y*p->y);
	float normx=0.0f,normy=0.0f,normz=0.0f;
	
	
	float n2;//n2 = layers_dc[new_layer].n;
	float r;
	float cos_angle_i;//cos_angle_i = fabsf(p->dz);

	if (new_layer==-10){
		n2=layers_dc[0].n; //N_LATERAL; //external refr_index
		normx=__fdividef(-(p->x),rad);
		normy=__fdividef(-(p->y),rad);
		cos_angle_i=fabsf(p->dx*normx+p->dy*normy);
	}else{
		n2=layers_dc[new_layer].n;
		normz=-copysignf(1.0f,p->dz);
		cos_angle_i = fabsf(p->dz);
	}
	
	if(n1==n2){//refraction index matching automatic transmission and no direction change
		p->layer = new_layer;
		return 0u;
	}

	if(n1>n2 && n2*n2<n1*n1*(1-cos_angle_i*cos_angle_i)){//total internal reflection, no layer change but z-direction mirroring
		p->dx+=2.0f*cos_angle_i*normx;
		p->dy+=2.0f*cos_angle_i*normy;
		p->dz+=2.0f*cos_angle_i*normz;
	  	p->t_left=0.0f;  
		//p->dz *= -1.0f;
		return 1u; 
	}

	if(cos_angle_i==1.0f){//normal incident
		r = __fdividef((n1-n2),(n1+n2));
		if(rand_MWC_co(x,a)<=r*r)
		{
			//reflection, no layer change but z-direction mirroring
			//p->dx+=2*cos_angle_i*normx;
			//p->dy+=2*cos_angle_i*normy;
			//p->dz+=2*cos_angle_i*normz;
			p->dx *= -1.0f;
			p->dy *= -1.0f;
			p->dz *= -1.0f;
			p->t_left=0.0f;
			return 1u;
		}
		else
		{	//transmission, no direction change but layer change
			//if(new_layer!=-10)
			//{
			//	p->layer = new_layer;
			//	return 0u;
			//}else 
			//{
			p->layer=new_layer;
			p->t_left=0.0f;
			return 0u;
			//}
		}
	}
	
	//gives almost exactly the same results as the old MCML way of doing the calculation but does it slightly faster
	// save a few multiplications, calculate cos_angle_i^2;
	float e = __fdividef(n1*n1,n2*n2)*(1.0f-cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
	r=2*sqrtf((1.0f-cos_angle_i*cos_angle_i)*(1.0f-e)*e*cos_angle_i*cos_angle_i);//use r as a temporary variable
	e=e+(cos_angle_i*cos_angle_i)*(1.0f-2.0f*e);//Update the value of e
	r = e*__fdividef((1.0f-e-r),((1.0f-e+r)*(e+r)));//Calculate r	

	if(rand_MWC_co(x,a)<=r)
	{ 
		// Reflection, mirror z-direction!
		//p->dz *= -1.0f;
		p->dx+=2.0f*cos_angle_i*normx;
		p->dy+=2.0f*cos_angle_i*normy;
		p->dz+=2.0f*cos_angle_i*normz;
		p->t_left=0.0f;
		return 1u;
	}
	else
	{
	  if(new_layer!=-10)
	  {
		// Transmission, update layer and direction
		r = __fdividef(n1,n2);
		e = r*r*(1.0f-cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
		p->dx *= r;
		p->dy *= r;
		p->dz = copysignf(sqrtf(1-e) ,p->dz);
		p->layer = new_layer;
		p->t_left=0.0f;
		return 0u;
	  }
	  else 
	  	{
	  	p->layer = new_layer;
	  	p->t_left=0.0f;
	  	return 0u;
	  	}
	}

}

//************************************* DetectAndSavePath ******************************************
__device__ void DetectAndSavePath(MemStruct* DeviceMem,PhotonStruct* p, int new_layer,int thd_id )
{
	float distance=sqrtf(p->x*p->x+p->y*p->y);
	unsigned int pos;
	
	// ******** DIFFUSE REFLECTANCE ***********
	if(new_layer == 0){ 		
		p->weight = 0u;// Set the remaining weight to 0, effectively killing the photon
		
		// Do I want to save Reflectance ?
		if (RorT_dc=='R'){ 
			// Does the photon hit a detector?
			if  ((distance>=detectors_dc[0]) && (distance<=detectors_dc[1])){
				pos=DeviceMem->received[thd_id];
				DeviceMem->received[thd_id]++; //incrementa n phot ricevuti nel thd
				for(short i_lay=0;i_lay<n_layers_dc[0];i_lay++)
				{ //salva cammini
					*(DeviceMem->path+(i_lay)*NUM_THREADS*N_photons_dc[0]+thd_id*N_photons_dc[0]+pos) = p->t[i_lay];
				}
			}
		} 
	}
	
	// ********* DIFFUSE TRANSMITTANCE *********					
	if(new_layer > *n_layers_dc){
		p->weight = 0u; // Set the remaining weight to 0, effectively killing the photon 
		// Do I want to save transmittance ?
		if (RorT_dc=='T'){
			if  ((distance>=detectors_dc[0]) && (distance<=detectors_dc[1])){
				pos=DeviceMem->received[thd_id];
				DeviceMem->received[thd_id]++; 
				for(short i_lay=0;i_lay<n_layers_dc[0];i_lay++)
				{ //salva cammini
					*(DeviceMem->path+(i_lay)*NUM_THREADS*N_photons_dc[0]+thd_id*N_photons_dc[0]+pos) = p->t[i_lay];
				}
			}							
 		} 
	}
	
	// PHOTON OUT LATERAL
	if(new_layer==-10)
	{
		p->weight=0u;  //uccide il fotone che esce lateralmente
	}
}

//***************************************************************************************************************
//Device function to add an unsigned integer to an unsigned long long using CUDA Compute Capability 1.1
__device__ void AtomicAddULL(unsigned long long* address, unsigned int add)
{
	if(atomicAdd((unsigned int*)address,add)+add<add)
		atomicAdd(((unsigned int*)address)+1,1u);
}

__device__ void AtomicAddUL(unsigned long* address, unsigned int add)
{
	if(atomicAdd((unsigned int*)address,add)+add<add)
		atomicAdd(((unsigned int*)address)+1,1u);
}
