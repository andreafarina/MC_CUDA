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
__global__ void MCd(MemStruct);
__device__ float rand_MWC_oc(unsigned long long*,unsigned int*);
__device__ float rand_MWC_co(unsigned long long*,unsigned int*);
__device__ void LaunchPhoton(PhotonStruct*,unsigned long long*, unsigned int*);
__global__ void LaunchPhoton_Global(MemStruct);
__device__ void Spin(PhotonStruct*, float,unsigned long long*,unsigned int*);
__device__ int Intersection(PhotonStruct* p, float* s);
__device__ unsigned int Reflect(PhotonStruct*, int, unsigned long long*, unsigned int*);
__device__ short int FindDetector(float distance);


//__device__ unsigned int PhotonSurvive(PhotonStruct*, unsigned long long*, unsigned int*);
__device__ void AtomicAddULL(unsigned long long* address, unsigned int add);
__device__ void AtomicAddUL(unsigned long* address, unsigned int add);


// __global__ void MCd(MemStruct DeviceMem)
//{
//	//Block index
//	int bx=blockIdx.x;
//
//	//Thread index
//	int tx=threadIdx.x;
//  //float* row;
//	
//
//	//First element processed by the block
//	int begin=NUM_THREADS_PER_BLOCK*bx;
//	
//
//	
//	//unsigned long long int x=DeviceMem.x[begin+tx];//coherent
//	//unsigned int a=DeviceMem.a[begin+tx];//coherent
//
//	//float s,distance;	//step length, distance
//	
//	//unsigned int i_det,i_lay;
//	
//	unsigned int index;//AF, w, index_old;
//	//AF index_old = 0;
//	//AF w = 0;
//	//AF unsigned int w_temp;
//	
//	//PhotonStruct p = DeviceMem.p[begin+tx];
//
//
//	//int new_layer;
//	
//	//First, make sure the thread (photon) is active
//	//AF unsigned int ii = 0;
//	//int term_photons_in_thread=0;
//	index=(begin+tx)*NUMSTEPS_GPU;
//	
//	
//	float* row=(float*)((char*)DeviceMem.path);
//	row[index]=100.;
//	float* row2=(float*)((char*)DeviceMem.path+DeviceMem.pitch);
//	row2[index]=200.;       
//	  
//}

	




__global__ void MCd(MemStruct DeviceMem)
 {
     //Block index
     int bx=blockIdx.x;
 
     //Thread index
     int tx=threadIdx.x;
	 
//unsigned long long int pos; //posizione salvataggio fotone
   
    
 
     //First element processed by the block
     int begin=NUM_THREADS_PER_BLOCK*bx;
 	
 
     
 	unsigned long long int x=DeviceMem.x[begin+tx];//coherent
 	unsigned int a=DeviceMem.a[begin+tx];//coherent
 
 	float s;	//step length
 	float *row;
	float cosNA=sqrtf(1.0f-NA*NA);
	
 	//unsigned long int PhotonThread=0;
	//short int i_lay;
 	//short int i_det;
 	//unsigned int index;//AF, w, index_old;
 	//AF index_old = 0;
 	//AF w = 0;
 	//AF unsigned int w_temp;
 	//unsigned int count_det[MAX_DETECTORS];
 	
	PhotonStruct p = DeviceMem.p[begin+tx];
 	
 
 	short int new_layer;
 	
 	
	
	for(int ii=0;ii<NUMSTEPS_GPU;ii++) //this is the main while loop
 	//do
 	{
 		if(layers_dc[p.layer].mutr!=FLT_MAX)
 			s = -__logf(rand_MWC_oc(&x,&a))*layers_dc[p.layer].mutr;//sample step length [cm] //HERE AN OPEN_OPEN FUNCTION WOULD BE APPRECIATED
 		else
 			s = 100.0f;//temporary, say the step in glass is 100 cm.
 		//printf("step_RND=%.5f \t",s);
 		//Check for layer transitions and in case, calculate s
 		//new_layer = p.layer;
 		//if(p.z+s*p.dz<layers_dc[p.layer].z_min){new_layer--; s = __fdividef(layers_dc[p.layer].z_min-p.z,p.dz);} //Check for upwards reflection/transmission & calculate new s
 		//if(p.z+s*p.dz>layers_dc[p.layer].z_max){new_layer++; s = __fdividef(layers_dc[p.layer].z_max-p.z,p.dz);} //Check for downward reflection/transmission
		new_layer=Intersection(&p,&s);
 		p.x += p.dx*s;
 		p.y += p.dy*s;
 		p.z += p.dz*s;
 		if ((bx==0)&&(tx==0)) printf("%.3f \t %.3f \t %.3f\n",p.x,p.y,p.z);
		//if (new_layer==-10) printf("cylinder\n");
 		//printf("step=%.5f \n",s);
		//p.t[p.layer-1] += s/LIGHT_SPEED*layers_dc[p.layer].n; //AF
 		p.t[p.layer-1] += s; //AF
 		
 		p.time_tot += s/LIGHT_SPEED*layers_dc[p.layer].n; //AF
 		if (p.time_tot>tmax_dc[0]) {
 		printf("Out of time\n");
 		p.weight=0u;
 		}
 		//if(p.z>layers_dc[p.layer].z_max)p.z=layers_dc[p.layer].z_max;//needed? mah...
 		//if(p.z<layers_dc[p.layer].z_min)p.z=layers_dc[p.layer].z_min;//needed?
 		if (p.weight==1u){
			if(new_layer!=p.layer){
				// set the remaining step length to 0
				s = 0.0f;  
				
				//Check for reflection
				if(Reflect(&p,new_layer,&x,&a)==0u){
					 			 
					// Photon is transmitted
					if(new_layer == 0){
						//Diffuse reflectance
								
						p.weight = 0u;// Set the remaining weight to 0, effectively killing the photon
						printf("%.3f \t %.3f \t %.3f\n",p.x,p.y,p.z);
						printf("Reflected\n");
						// Do I want to save Reflectance ?
						if (RorT_dc=='R'){ 
							short i_det;
							// Does the photon hit a detector?
							
							i_det=FindDetector(sqrtf(p.x*p.x+p.y*p.y));
												
							if (i_det>-1){
								//__syncthreads();
								printf("Detected\n");
								if (DeviceMem.photons[i_det]<N_photons_dc[0]){//N_photons_dc[0]){//occhio!rischio di trascurare gia fotoni utili
									//AtomicAddUL(&DeviceMem.photons[i_det], 1);
									//atomicExch(&pos,DeviceMem.photons[i_det]);
						 
									pos=DeviceMem.photons[i_det];
						 
									DeviceMem.photons[i_det]++;
						 
									for(short i_lay=0;i_lay<n_layers_dc[0];i_lay++){
					
										//float* row=(float*)((char*)DeviceMem.path+(n_layers_dc[0]*i_det+i_lay)*DeviceMem.pitch);
										row=(float*)((char*)DeviceMem.path+(n_layers_dc[0]*i_det+i_lay)*DeviceMem.pitch);
										//atomicExch(&row[pos],p.t[i_lay]);
										row[pos]=p.t[i_lay];
										//row[index+count_det[i_det]]=p.t[i_lay];//*LIGHT_SPEED/layers_dc[i_lay+1].n;
									}
								}
							}
								// set the remaining step length to 0
								//s = 0.0f;
						} // end if (RorT=='R')
					} //end if (new_layer == 0)
						
					
					//Photon is transmitted
					if(new_layer > *n_layers_dc){
						// Diffuse Transmittance
 						printf("%.3f \t %.3f \t %.3f\n",p.x,p.y,p.z);
 						printf("Transmitted\n");
						p.weight = 0u; // Set the remaining weight to 0, effectively killing the photon 
						// Do I want to save transmittance ?
						if (RorT_dc=='T'){
							short i_det;
							// Does the photon hit a detector?
							//if(fabsf(p.dz)>cosNA)
							//{
							i_det=FindDetector(sqrtf(p.x*p.x+p.y*p.y));
				   
							if (i_det>-1){//(i_det>=0)){
								printf("Detected\n");
								if (DeviceMem.photons[i_det]<N_photons_dc[0]){//occhio!rischio di trascurare gia fotoni utili
																	
									 //AtomicAddUL(&DeviceMem.photons[i_det], 1);
									
									pos=DeviceMem.photons[i_det];
									DeviceMem.photons[i_det]++;
								
									for(short i_lay=0;i_lay<n_layers_dc[0];i_lay++){
										float* row=(float*)((char*)DeviceMem.path+(n_layers_dc[0]*i_det+i_lay)*DeviceMem.pitch);
										row[pos]=p.t[i_lay];
										//row[index+count_det[i_det]]=p.t[i_lay];//*LIGHT_SPEED/layers_dc[i_lay+1].n;
									}
								}
							}
							//}
						} //end if (RorT_dc=='T')
					} //end if(new_layer > *n_layers_dc)
					if(new_layer==-10){
					printf("Out lateral\n");
					p.weight=0u;  //uccide il fotone che esce lateralmente
					}
				} // end if(Reflect(&p,new_layer,&x,&a)==0u)
			} //end if(new_layer!=p.layer)
			else Spin(&p,layers_dc[p.layer].g,&x,&a);
		} // end if (p.weight==1)
 		//w=0;
 		//if(s > 0.0f) 
 		//  Spin(&p,layers_dc[p.layer].g,&x,&a);
		  
 		if ((p.weight==0u)){//&&(DeviceMem.photons[0]<N_photons_dc[0])){ 
 		  printf("New Photon \n");
 		  LaunchPhoton(&p,&x,&a);
 		  //AtomicAddULL(DeviceMem.ph_launched, 1);
		  //DeviceMem.ph_launched=DeviceMem.ph_launched+1;
		  //(*DeviceMem.ph_launched)=(*DeviceMem.ph_launched)+1;
		  
		  //DeviceMem.counter[begin+tx]=DeviceMem.counter[begin+tx]+1;
		  //AtomicAddUL(&DeviceMem.counter[begin+tx], 1);
 		}
 		
 	}//end for....NUMSTEP_GPU   //while(term_photons_in_thread<n_detectors_dc[0]*NUMSTEPS_GPU);
 	
 	//save the state of the MC simulation in global memory before exiting
 	DeviceMem.p[begin+tx] = p;	//This one is incoherent!!!
 	DeviceMem.x[begin+tx] = x;  //this one also seems to be coherent
 	
 
 }//end MCd




__device__ void LaunchPhoton(PhotonStruct* p, unsigned long long* x, unsigned int* a)
{
	// We are currently not using the RNG but might do later
	//float input_fibre_radius = 0.03;//[cm]
	//p->x=input_fibre_radius*sqrtf(rand_MWC_co(x,a));

	p->x  = 0.0f;
	p->y  = 0.0f;
	p->z  = 0.0f;
	p->dx = 0.0f;
	p->dy = 0.0f;
	p->dz = 1.0f;
	//float RAD_FIBER=0.05f; //cm
	//float NA=0.38f; //num_aperture
// 	float rad;
// 	float cosa,sina,sint,cost,sinp,cosp;
// 	
// 	rad=RAD_FIBER*rand_MWC_co(x,a);
// 	__sincosf(2.0f*PI*rand_MWC_co(x,a),&sina,&cosa); //coordinate polari nella testa della fibra
// 
// 	//Photon starting coordinates inside the fiber
// 	p->x  = rad*cosa;
// 	p->y  = rad*sina;
//  	p->z  = 0.0f;
// 	__sincosf(asinf(NA)*rand_MWC_co(x,a),&sint,&cost);// spin psi [0-2*PI)
// 	__sincosf(2.0f*PI*rand_MWC_co(x,a),&sinp,&cosp);
// 	
// 	p->dx = sint*cosp;
// 	p->dy = sint*sinp;
// 	p->dz = cost;
	
	for (short int i=0;i<n_layers_dc[0];i++)
	p->t[i] = 0.0f;
	p->time_tot=0.0f;
	p->layer = 1;
	p->weight=1u;
	
	
//AF	p->weight = *start_weight_dc; //specular reflection!

}

__global__ void LaunchPhoton_Global(MemStruct DeviceMem)//PhotonStruct* pd, unsigned long long* x, unsigned int* a)
{
	int bx=blockIdx.x;
    int tx=threadIdx.x;	

    //First element processed by the block
    int begin=NUM_THREADS_PER_BLOCK*bx;

	PhotonStruct p;
	unsigned long long int x=DeviceMem.x[begin+tx];//coherent
	unsigned int a=DeviceMem.a[begin+tx];//coherent

	LaunchPhoton(&p,&x,&a);
	//AtomicAddULL(DeviceMem.ph_launched, 1);

	//__syncthreads();//necessary?
	DeviceMem.p[begin+tx]=p;//incoherent!?
	//DeviceMem.counter[begin+tx]=DeviceMem.counter[begin+tx]+1;
	//DeviceMem.ph_launched++;
}


__device__ void Spin(PhotonStruct* p, float g, unsigned long long* x, unsigned int* a)
{
	float cost, sint;	// cosine and sine of the 
						// polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the 
						// azimuthal angle psi. 
	float temp;

	float tempdir=p->dx;

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

			

// 29-11-2011 MODIFICHE CYLINDER
__device__ int Intersection(PhotonStruct* p, float* s)
{
	//float RADIUS=2000.54f/2.0f; //cylinder radius (cm)
	//float s_cyl;
	float A,B,C;
	short int new_layer=p->layer;
	//if(p->z+(*s)*p->dz<layers_dc[p->layer].z_min){new_layer--; *s = __fdividef(layers_dc[p->layer].z_min-p->z,p->dz);} //Check for upwards reflection/transmission & calculate new s
 	//if(p->z+(*s)*p->dz>layers_dc[p->layer].z_max){new_layer++; *s = __fdividef(layers_dc[p->layer].z_max-p->z,p->dz);} //Check for downward reflection/transmission
	
	 
 	if((p->x+(*s)*p->dx)*(p->x+(*s)*p->dx)+(p->y+(*s)*p->dy)*(p->y+(*s)*p->dy)>RADIUS*RADIUS) //Check for lateral boundary hit
 	{
	  A=1.0f-(p->dz)*(p->dz);
	  B=(p->dx)*(p->x)+(p->dy)*(p->y);
	  C=(p->x)*(p->x)+(p->y)*(p->y)-RADIUS*RADIUS;
	  //printf("A=%.3f,B=%.3f,C=%.3f\n",A,B,C);
	  (*s)=fdividef(-B+sqrtf(B*B-A*C),A);
	  //if (s_cyl<0) s_cyl=(-B-sqrtf(B*B-A*C))/A;
	  new_layer=-10;
	}
	  
	// if(p->z+(*s)*p->dz<layers_dc[p->layer].z_min){new_layer--; *s = __fdividef(layers_dc[p->layer].z_min-p->z,p->dz);} //Check for upwards reflection/transmission & calculate new s
//  	if(p->z+(*s)*p->dz>layers_dc[p->layer].z_max){new_layer++; *s = __fdividef(layers_dc[p->layer].z_max-p->z,p->dz);} //Check for downward reflection/transmission
	if(p->z+(*s)*p->dz<layers_dc[p->layer].z_min){
		*s = __fdividef(layers_dc[p->layer].z_min-p->z,p->dz);
		return p->layer-1;
	} //Check for upwards reflection/transmission & calculate new s
 	if(p->z+(*s)*p->dz>layers_dc[p->layer].z_max){
 		*s = __fdividef(layers_dc[p->layer].z_max-p->z,p->dz);
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

	if (new_layer==-10)
	{
		n2=N_LATERAL; //external refr_index
		normx=__fdividef(-(p->x),rad);
		normy=__fdividef(-(p->y),rad);
		cos_angle_i=fabsf(p->dx*normx+p->dy*normy);
	}else
	{
		n2=layers_dc[new_layer].n;
		normz=-copysignf(1.0f,p->dz);
		cos_angle_i = fabsf(p->dz);
	}
	
	if(n1==n2)//refraction index matching automatic transmission and no direction change
	{	
		p->layer = new_layer;
		return 0u;
	}

	if(n1>n2 && n2*n2<n1*n1*(1-cos_angle_i*cos_angle_i))//total internal reflection, no layer change but z-direction mirroring
	{
		p->dx+=2.0f*cos_angle_i*normx;
		p->dy+=2.0f*cos_angle_i*normy;
		p->dz+=2.0f*cos_angle_i*normz;
	  	  
		//p->dz *= -1.0f;
		return 1u; 
	}

	if(cos_angle_i==1.0f)//normal incident
	{		
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
		return 0u;
	  }
	  else 
	  	{
	  	p->layer = new_layer;
	  	return 0u;
	  	}
	}

}

__device__ short int FindDetector(float distance) //Return -1 if the photon doesn't hit the detector, return detector number 0,1,2,...
{
	char detected='0';
	for(short int i=0;(i<n_detectors_dc[0])&&(detected=='0');i++){
		if  ((distance>=detectors_dc[i+i]) && (distance<=detectors_dc[i+i+1])){
			detected='1';
			return i; //Hit detector number
		}
	}
	return -1;

}


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
