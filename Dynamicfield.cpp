///*Dynamic SCF of polymer adsorption*/ 
#include <math.h>
#include <stdio.h>
//COMPILE: g++ PolAds.cpp -lm -o PolAds     RUN: ./PolAds

double tolerance = 1e-5, eta=0.2, error = 1;
int it = 0, itmax=1000;
int z, n_layers = 50, jz = n_layers+2;
int s, N = 100;
double phib_solvent, phib = 0.01, chi_s = -1;

int main() {
  double phi[jz], phi_solvent[jz];
  double rho[jz];
  double u[jz],us[jz];
  double G1[jz], alpha[jz], alpha_s[jz], G[N][jz];
  double L[jz],bU1[jz],bU2[jz];

  phib_solvent = 1.0-phib;
  for (z = 1; z <= n_layers; z++) {u[z] = 0; phi[z]=phib; phi_solvent[z]=phib_solvent;}
  G[0][0]=0; 
  double time=0.0, endtime=0.002;
  while(time < endtime){
	time=time+0.001;
	for (z = 1; z <= n_layers; z++) us[z]=-log(phi_solvent[z]/phib_solvent);
	error=1; it =0;

  	while (error>tolerance && it < itmax) {
		error=0; it++;
	
		for (z = 1; z <= n_layers; z++) {G1[z]=exp(-u[z]); G[0][z]=G1[z];}
		G[0][n_layers+1]=G[0][n_layers];

   		for (s = 1; s < N; s++) {
   	   		for (z=1; z<=n_layers; z++) G[s][z]=G1[z]*(G[s-1][z-1]+4*G[s-1][z]+G[s-1][z+1])/6;
    	  		G[s][0]=0; G[s][n_layers+1]=G[s][n_layers];
   		 }
   		 for (z = 1; z <= n_layers; z++){
			 rho[z]=0; for (s = 0; s < N; s++) rho[z] += phib/N * G[s][z]*G[N-s-1][z]/G1[z];
  			 u[z] += eta*(rho[z]-phi[z]);
 			 error +=  pow((rho[z]-phi[z]),2);    
  			 error = sqrt(error); 
			// printf("z = %i error = %1e \n",z,error);
   		 }
		 printf("it = %i error = %1e \n",it,error);	
  	}  
			
	for (z=0;z<=n_layers+1; z++){
		L[z]=phi[z]*phi_solvent[z]; 
		alpha_s[z]=us[z];
		if(z == 1) alpha[z]=u[z]+chi_s; else alpha[z]=u[z]; 
		bU1[z]=alpha[z]-alpha_s[z];
		bU2[z]= -1.0*bU1[z];
	}
	double J1[jz],J2[jz];
	for(z=1;z<=n_layers; z++){
		J1[z]=0.50*(((L[z]+L[z+1])*(bU1[z+1]-bU1[z]))-((L[z-1]+L[z])*(bU1[z]-bU1[z-1])));
		J2[z]=0.50*(((L[z]+L[z+1])*(bU2[z+1]-bU2[z]))-((L[z-1]+L[z])*(bU2[z]-bU2[z-1])));
	} 
	for (z=1; z<=n_layers;  z++){
		phi[z]+=0.001*J1[z];
		phi_solvent[z]+=0.001*J2[z];
	}

  }
  for (z = 1; z <= n_layers; z++)
  printf("z = %i phi = %1f phi_sol = %1f \n",z, phi[z], phi_solvent[z]);

/*  
//*******************************************************************************************
	//begin time loop
		

	double time=0.0, endtime=1;
	while(time < endtime){
		time=time+0.001;
		//get density
		//Density should be obtained as initial guess, satisfying compressibility. 
		//This density should be used to obtain the segment potentials
		//compute the potential
		for(z=1; z<=n_layers; z++) {
			u[z]=chi_s*((phi_solvent[z-1]+phi_solvent[z]+phi_solvent[z+1])/3-phib_solvent);
			us[z]=chi_s*((phi[z-1]+phi[z]+phi[z+1])/3-phib);
		}
		//compute the flux (Use Onsager kinetic coeffecients)
		double L[jz],bU1[jz],bU2[jz];
		for (z=0;z<=n_layers+1; z++){
			L[z]=phi[z]*phi_solvent[z]; //L=D*phia*phib (What is the Diffusion coeffecient?)
			bU1[z]=u[z]-us[z];
			bU2[z]= -1.0*bU1[z];
			}
		double J1[jz],J2[jz];
		for(z=1;z<=n_layers; z++){
			J1[z]=-0.50*(((L[z]+L[z+1])*(bU1[z+1]-bU1[z]))-((L[z-1]+L[z])*(bU1[z]-bU1[z-1])));
			J2[z]=-0.50*(((L[z]+L[z+1])*(bU2[z+1]-bU2[z]))-((L[z-1]+L[z])*(bU2[z]-bU2[z-1])));
		}
		//update density (using Cahn-Hilliard langevin dynamcis)
		//There are two density fields for polymer and solvent respectively. 
		for (z=1; z<=n_layers;  z++){
			phi[z]+=0.001*J1[z];
			phi_solvent[z]+=0.001*J2[z];
		}
	printf("time: %1f \t phi[1]: %1f \t phis[1]: %1f \n", time, phi[1], phi_solvent[1]);
		//end time loop
	}
//********************************************************************************************

*/

  return 0;
};
