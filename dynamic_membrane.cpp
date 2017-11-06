//Standard header files
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;

//variable declaration
int M = 50; //system size in lattice layers, equivalent to n\_layers =50
int N = 37; //2x16 + 5 //two tails + 1 head: total length of the lipid.
double theta = 4.15165700; //amount of lipid in the system for tensionless bilayer
double chi=1.6; //interaction parameter for polar-apolar interactions
double tolerance = 1e-5;
double q, phib, phib_S;

//Main program starts here
main() {
	int z,s,jz=M+2;
	double phi_t[jz], phi_h[jz], phi_S[jz];
	double rho_t[jz], rho_h[jz], rho_S[jz];
	double et[jz],eh[jz],eS[jz], Lht[jz], Lhs[jz], Lst[jz], bU1[jz],bU2[jz], bU3[jz], bU4[jz], J1[jz], J2[jz], J3[jz];
	double u_t[jz], u[jz], u_h[jz], w_t[jz], w[jz], G1_t[jz], G1[jz], alpha[jz];
	int delta[38]; //Kronicker delta’s for specifying chain architecture.310

	for (s=1; s<=37; s++) delta[s]=0; for (s=17; s<=21; s++) delta[s]=1; 

	ifstream fp, pp;
	double data[150];
	fp.open("memguess.dat"); pp.open("potentialguess.dat");
	for (z = 1; z<=M; z++){
		for(s=1; s<=3; s++){ 
			if(s==1) fp >> phi_t[z]; 
			else if(s==2) fp >> phi_h[z]; 
			else  fp >> phi_S[z]; 
		}
	}
	for (z=1; z<=M; z++){
		for(s=1; s<=3; s++){
			if(s==1) pp >> u_t[z];
			else if(s==2) pp >> u[z];
			else pp >> alpha[z];
		}
	}
	phi_t[0]=phi_t[1]; phi_t[M+1]=phi_t[M]; //bc
	phi_h[0]=phi_h[1]; phi_h[M+1]=phi_h[M]; //bc
	phi_S[0]=phi_S[1]; phi_S[M+1]=phi_S[M]; //bc
	fp.close(), pp.close();

	double G_f[jz*N], G_b[jz*N];
	double a,b,c,gz;
	for (z=0; z<=M; z++){
	alpha[z] = 0; u_t[z]=0; u[z]=0;}

	// rough guess the density
	// begin time loop
	int time=0, endtime=100;
	while (time < endtime){
	time++;
		//begin loop to find the chemical potential corresponding to density k
		double error = 1, eta=0.082;
		int it=0;
		while (error > tolerance && it < 10000) {
			it++;
			s=1;
			for (z=1; z<=M; z++) {
				w_t[z] = u_t[z]; w[z] = u[z];
				G1_t[z]=exp(-alpha[z] - u_t[z]); G_f[z] = G1_t[z];
				G1[z] = exp(-alpha[z] - u[z]);
				}
				G_f[0]=G_f[1]; G_f[M+1]=G_f[M]; //bc
			for (s=2; s<=N; s++) {
				b=G_f[jz*(s-2)]; c=G_f[jz*(s-2)+1];
				for (z=1; z<=M; z++) {
				if (delta[s]==0) gz=G1_t[z]; else gz=G1[z]; //select segment weight
				a=b; b=c; c=G_f[jz*(s-2)+z+1]; 
				G_f[jz*(s-1)+z] = gz*(a+b+c)/3; //propagator with lambda = 1/3.
				}
				G_f[jz*(s-1)]=G_f[jz*(s-1)+1];//bc
				G_f[jz*(s-1)+M+1]=G_f[jz*(s-1)+M]; //bc
			}
			s=N;
			for (z=1; z<=M; z++) G_b[jz*(s-1)+z]=G1_t[z];
				G_b[jz*(s-1)+0]=G_b[jz*(s-1)+1]; G_b[jz*(s-1)+M+1]=G_b[jz*(s-1)+M]; //bc
			for (s=N-1; s>=1; s--) {
				b=G_b[jz*s]; c=G_b[jz*s+1];
				for (z=1; z<=M; z++) {
				if (delta[s]==0) gz=G1_t[z]; else gz=G1[z];
				a=b; b=c; c=G_b[jz*s+z+1];
				G_b[jz*(s-1)+z] = gz * (a + b + c)/3;
				}
				G_b[jz*(s-1)]=G_b[jz*(s-1)+1]; G_b[jz*(s-1)+M+1]=G_b[jz*(s-1)+M]; //bc
			}
			q = 0; for (z=1; z<=M; z++) q+=G_b[z];
			phib = theta/q; phib_S = 1-phib;
			for (z=1; z<=M; z++){
				rho_S[z] = phib_S*G1[z]; rho_t[z]=0; rho_h[z]=0;
				for (int s=1; s<=N; s++) { //composition law: collect phi per segment type
				if (delta[s]==0) rho_t[z]+=theta/(q*N)*G_f[jz*(s-1)+z]*G_b[jz*(s-1)+z]/G1_t[z];
				else rho_h[z]+=theta/(q*N)*G_f[jz*(s-1)+z]*G_b[jz*(s-1)+z]/G1[z];
				}
				}
				rho_t[0]=rho_t[1]; rho_t[M+1]=rho_t[M]; //bc
				rho_h[0]=rho_h[1]; rho_h[M+1]=rho_h[M]; //bc
				rho_S[0]=rho_S[1]; rho_S[M+1]=rho_S[M]; //bc
				error=0;
			for (z=1; z<=M; z++) {
				u_t[z]+=eta*(rho_t[z]-phi_t[z]);
				error += pow((rho_t[z]-phi_t[z]),2);
				u[z]+=eta*(rho_S[z]-phi_S[z]);
				error += pow((rho_S[z]-phi_S[z]),2);
				u_h[z]+=eta*(rho_h[z]-phi_h[z]);
				error += pow((rho_h[z]-phi_h[z]),2);
				//alpha[z] += eta*(rho_t[z]+rho_h[z]+rho_S[z]-1);
				//error += pow((rho_t[z]+rho_h[z]+rho_S[z]-1),2);
			}
			error=sqrt(error);
			//printf("it = %i error = %1e \n",it,error);
		}
	// compute alpha for all segment types
	// we do not update the density of head group as it is solved based on compressibility relation. Initial guess also satisfies compressibilty
	for (z=1; z<=M; z++){
		et[z]=u_t[z]-chi*((phi_h[z-1]+phi_h[z]+phi_h[z+1])/3-phib*5/37+(phi_S[z-1]+phi_S[z]+phi_S[z+1])/3-phib_S);
		eS[z]=u[z]-chi*((phi_t[z-1]+phi_t[z]+phi_t[z+1])/3-phib*32/37);
		eh[z]=u_h[z]-chi*((phi_t[z-1]+phi_t[z]+phi_t[z+1])/3-phib*32/37);
	}

	//onsagers coefficients
	for (z=1; z<=M; z++){
		//Diffusion coefficients are probably wrong
		Lhs[z]=phi_h[z]*phi_S[z];// + phi_h[z]*phi_t[z];
		Lht[z]=phi_h[z]*phi_t[z];// + phi_t[z]*phi_h[z];
		Lst[z]=phi_S[z]*phi_t[z];// + phi_S[z]*phi_t[z];
		//this seems ok as u3 = -(u2+u1). but should be checked too. 
		bU1[z]=et[z]-eS[z];// + et[z]-eh[z];
		bU2[z]=et[z]-eh[z];// + eS[z]-eh[z];
		bU3[z]=eS[z]-et[z];// + eh[z]-eS[z];
		bU4[z]=eS[z]-eh[z];//
		}
	// boundary conditions to onsager coeff. required for flux calculation
	Lhs[M+1]=Lhs[M]; Lht[M+1]=Lht[M]; Lst[M+1]=Lst[M]; bU1[M+1]=0; bU2[M+1]=0; bU3[M+1]=0; bU4[M+1]=0;
	Lhs[0]=Lhs[1]; Lht[0] = Lht[1]; Lst[0]=Lst[1]; bU1[0]=bU1[1]; bU2[0]=bU2[1]; bU3[0]=bU3[1]; bU4[0]=bU4[1];
	// flux calculation and density update.
	for (z=1; z<=M; z++) {
		J1[z]=-0.50*(((Lst[z]+Lst[z+1])*(bU1[z+1]-bU1[z]))-((Lst[z-1]+Lst[z])*(bU1[z]-bU1[z-1])));
		J1[z] += -0.50*(((Lht[z]+Lht[z+1])*(bU2[z+1]-bU2[z]))-((Lht[z-1]+Lht[z])*(bU2[z]-bU2[z-1]))); 
		J2[z]=-0.50*(((Lst[z]+Lst[z+1])*(bU3[z+1]-bU3[z]))-((Lst[z-1]+Lst[z])*(bU3[z]-bU3[z-1])));
		J2[z] +=-0.50*(((Lhs[z]+Lhs[z+1])*(bU4[z+1]-bU4[z]))-((Lhs[z-1]+Lhs[z])*(bU4[z]-bU4[z-1])));//J3 violates compressibility, probably because of wrong onsager coeff. 
		phi_t[z]+=0.001*J1[z];
		phi_S[z]+=0.001*J2[z];
		phi_h[z]-=0.001*(J1[z]+J2[z]);
		}
	if(time == 3 )for (z=1; z<=M; z++) printf("z = %i phi_t = %1f phi_h = %1f phi_S = %1f\n",z, phi_t[z], phi_h[z], phi_S[z]);
	}
	for (z=1; z<=M; z++) printf("z = %i phi_t = %1f phi_h = %1f phi_S = %1f\n",z, phi_t[z], phi_h[z], phi_S[z]);
	return(0);
};
