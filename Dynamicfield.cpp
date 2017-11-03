#include <math.h>
#include <stdio.h>

double tolerance = 1e-5, eta=0.1, error = 1;
int it = 0, itmax=10000;
int z, n_layers = 50, jz = n_layers+2;
int s, N = 10;
double phib_solvent, phib = 0.01, chi_s = -1;

int main() {
  double phi[jz], phi_solvent[jz];
  double rho[jz];
  double u[jz],us[jz];
  double G1[jz], alpha[jz], alpha_s[jz], G[N][jz];
  double L[jz],bU1[jz],bU2[jz];
  double J1[jz],J2[jz];

  phib_solvent = 1.0-phib;
  for (z = 1; z <= n_layers; z++) {u[z] = 0; phi[z]=phib; phi_solvent[z]=phib_solvent;}
  G[0][0]=0; 

  // boundary conditions for densities. Adsorption problem has zero density at wall.
  phi[0]=0; phi_solvent[0]=0;
  phi[n_layers+1]=phi[n_layers]; phi_solvent[n_layers+1]=phi_solvent[n_layers];	

  int count=0;
  double time=0.0, endtime=1000000, theta;

  /* Time evolution*/
  while(time < endtime){
	time=time+0.1;
        count++;
	for (z = 1; z <= n_layers+1; z++) us[z]=-log(phi_solvent[z]/phib_solvent);
	error=1; it =0;
        /* Inverse problem of u to match guessed phi */
  	while (error>tolerance && it<itmax) {
		error=0; it++;
	
		for (z = 1; z <= n_layers; z++) {G1[z]=exp(-u[z]); G[0][z]=G1[z];}
		G[0][n_layers+1]=G[0][n_layers];

   		for (s = 1; s < N; s++) {
   	   		for (z=1; z<=n_layers; z++) G[s][z]=G1[z]*(G[s-1][z-1]+4*G[s-1][z]+G[s-1][z+1])/6;
    	  		G[s][0]=0; G[s][n_layers+1]=G[s][n_layers];
   		 }
   		 for (z = 1; z <= n_layers; z++){
			 rho[z]=0; for (s = 0; s < N; s++) rho[z] += phib/N * G[s][z]*G[N-s-1][z]/G1[z];
		 }
		 rho[0]=0; rho[n_layers+1]=rho[n_layers];
		 for (z=1; z<=n_layers; z++){
  			 u[z] += eta*(rho[z]-phi[z]);
 			 error +=  pow((rho[z]-phi[z]),2);
		//	 printf("z = %i error = %1e phi: %1e rho: %1e \n",z,error, phi[z], rho[z]);
		 }    
  		 error = sqrt(error); 
//		 if(it%1000 ==0) printf("it = %i errr = %1e \n",it,error);	
  	}  
	theta=0;
	for (z=1; z<=n_layers; z++){
        theta += phi[z] ; 
  	}
   	if(count%1000 ==0) printf("time: %1e error: %1e theta: %1e\n", time, error, theta);
	for (z=1;z<=n_layers; z++){
		L[z]=phi[z]*phi_solvent[z]; 
		alpha_s[z]=us[z];
		if(z == 1) alpha[z]=u[z]-chi_s; else alpha[z]=u[z]; 
		bU1[z]=alpha[z]-alpha_s[z];
		bU2[z]= -1.0*bU1[z];
	}

	L[n_layers+1]=phib*phib_solvent; bU1[n_layers+1]=0; bU2[n_layers+1]=0;
        L[0]=L[1]; bU1[0]=bU1[1]; bU2[0]=bU2[1];
	for(z=1;z<=n_layers; z++){
		J1[z]=-0.50*(((L[z]+L[z+1])*(bU1[z+1]-bU1[z]))-((L[z-1]+L[z])*(bU1[z]-bU1[z-1])));
		J2[z]=-0.50*(((L[z]+L[z+1])*(bU2[z+1]-bU2[z]))-((L[z-1]+L[z])*(bU2[z]-bU2[z-1])));
	} 
	//flux at boundaries
//	J1[1]=0.5*(L[2]-L[1])*(bU1[2]-bU1[1]); J2[1]=-1*J1[1];

	for (z=1; z<=n_layers;  z++){
		phi[z]+=0.1*J1[z];
		phi_solvent[z]+=0.1*J2[z];
	}
  }
  for (z = 1; z <= n_layers; z++)
  printf("z = %i phi = %1f phi_sol = %1f \n",z, phi[z], phi_solvent[z]);
  printf("alpha: %1e uzero: %1e Lzero: %1e \n", alpha[0], u[0], L[0]);
  return 0;
};
