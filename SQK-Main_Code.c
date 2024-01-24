/*
In practice it is often advantageous to work with a lattice where the total number of sites is power of 2 (Page - 334, Newman);
*/


#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

float Sx(float, float);
float Sy(float, float);
float Sz(float, float);
float Hxxz(float, float, float, float, float);
float Hdm(float, float, float, float);


#define PI 3.14159265


int main () 
{	
	int N, L, MCS, Niter, i, j, k,  mm, nn, a1, b1, c1, d1, a2, b2, c2, d2, pp, qq;
	float sign_a, sign_b, sign_c, sign_d;
	float  J_ising = 1.0, DD, delta, T, E, Ei, Ef, Mx, My, Mz, M, dE, u;
	float r, dummy_theta, dummy_phi, dtheta, dphi;
	float DD_list[11], DELTA_list[11];
	int seed = time(NULL);                 // This time function returns the number of seconds ellapsed since Jan 1, 1970
										     //  So this will provide different seed value at each run
	int count;
	
	float Temp[8] = {2.0, 1.5, 1.0, 0.5, 0.2, 0.1, 0.05, 0.01};
	
	for(i = 0; i<11; i++)
	{
		DD_list[i] = (float)(i-5);
		DELTA_list[i] = (float)(i-5); 
	}
										    
	FILE *file01;
                
	file01 = fopen("/home/sumank/project codes/square_kagome_angles.txt", "w+");
	
	srand(seed);
	
	
   	L = 10;
   	Niter = 70000;

   	N = 6*L*L ;                   // Number of sites in a unit cell is 6 in square-kagome lattice
   	dtheta = 0.1*PI;
   	dphi = 0.1*PI;
   	
   	float theta[L][ 6*L], phi[L][6*L] ;
   	float av_theta[L][6*L], av_phi[L][6*L] ;
   	
   	
   	
   	for(pp = 0; pp<11; pp++)
   	{
   		for(qq = 0; qq<11; qq++)
   		{
   			
   			delta = DELTA_list[pp]; DD = DD_list[qq];
   			
   			count = 0;

			for(i = 0; i < L; i++)
   			{	
   				for(j=0; j<6*L; j++)
   				{
   					r = (float) rand() / RAND_MAX;
   					theta[i][j] = r*PI;
					r = (float) rand() / RAND_MAX;
   					phi[i][j] = 2.0*r*PI;	
   				}
   				
   			}
   			//phi[0][0] = 0.0;
   	
   	
   			for(i=0; i<L; i++)
			{	
				for(j=0; j<6*L; j++)
				{
					av_theta[i][j] = 0.0; av_phi[i][j] = 0.0;
				}
			}
 	
 	
			E = 0.0 ; Mx = 0.0 ; My = 0.0; Mz = 0.0;

			for(i=0; i<L; i++)
			{	
				for(j=0; j<6*L; j++)
				{
					if(j%6 ==0)											 // Defining the Nearest Neighbours with PBC
					{
						a1 = i; a2 = (j+1);                                               
						b1 = i; b2 = (j+3); 
						c1 = i; c2= (j+4); 
						d1 = i; d2 = (j + 6*L -1)%(6*L);   
						sign_a = -1.0;
						sign_b = 1.0;
						sign_c = 1.0;
						sign_d = -1.0; 
					} 
					else if(j%6 == 1)
					{
						a1 = i; a2 = (j+1);                                                
						b1 = i; b2 = (j+3); 
						c1 = i; c2 = (j+4); 
						d1 = i; d2 = (j-1);    
						sign_a = -1.0;
						sign_b = -1.0;
						sign_c = 1.0;
						sign_d = 1.0; 
					}
					else if(j%6 == 2)
					{
						a1 =i; a2 = (j+1);                                                
						b1 = i; b2 = (j + 3); 
						c1 = (i+1)%L; c2 = (j+2);    
						d1 = i; d2 = (j-1); 
						sign_a = -1.0;
						sign_b = -1.0;
						sign_c = 1.0;
						sign_d = 1.0;    
		   			}
		   			else if(j%6 == 3)
					{
						a1 = i; a2 = (j-1);                                                
						b1 = i; b2 = (j-3); 
						c1 = (i+1)%L; c2 = (j+1);
						d1 = i; d2 = (j+6*L -4)%(6*L);
						sign_a = 1.0;
						sign_b = -1.0;
						sign_c = -1.0;
						sign_d = 1.0;    
		   			}
		   			else if(j%6 == 4)
					{
						a1 = i; a2 = (j-3);                                                
						b1 = i; b2 = (j-4);
						c1 = (i-1+L)%L; c2 = (j-1);
						d1 = (i-1+L)%L; d2 = (j-2);
						sign_a = 1.0;
						sign_b = -1.0;
						sign_c = 1.0;
						sign_d = -1.0;    
		   			}
		   			else if(j%6 == 5)
					{
						a1 = i; a2 = (j-3) ;                                          
						b1 = i; b2 = (j-4);
						c1 = i; c2 = (j+1)%(6*L);
						d1 = i; d2 = (j+4)%(6*L);
						sign_a = 1.0;
						sign_b = -1.0;
						sign_c = 1.0;
						sign_d = -1.0;    
		   			}
			
				}
		
				E = E - J_ising*(    Hxxz(theta[i][j], phi[i][j], theta[a1][a2], phi[a1][a2], delta) 
							+ Hxxz(theta[i][j], phi[i][j], theta[b1][b2], phi[b1][b2], delta) 
							+ Hxxz(theta[i][j], phi[i][j], theta[c1][c2], phi[c1][c2], delta) 
							+ Hxxz(theta[i][j], phi[i][j], theta[d1][d2], phi[d1][d2], delta) )
					  		+DD*( sign_a*Hdm(theta[i][j], phi[i][j], theta[a1][a2], phi[a1][a2]) 
							+ sign_b*Hdm(theta[i][j], phi[i][j], theta[b1][b2], phi[b1][b2]) 
							+ sign_c*Hdm(theta[i][j], phi[i][j], theta[c1][c2], phi[c1][c2]) 
							+ sign_d*Hdm(theta[i][j], phi[i][j], theta[d1][d2], phi[d1][d2]) )  ; 
							
				Mx = Mx + Sx(theta[i][j], phi[i][j]);   
				My = My + Sy(theta[i][j], phi[i][j]);
				Mz = Mz + Sz( theta[i][j], phi[i][j]);
		               
			}
	
			E = E*0.5 ;
			M = sqrt( Mx*Mx + My*My + Mz*Mz );

			// INITIALIZATION COMPLETED 
	
	
			//EVOLVING THE SYSTEM TO LOWEST TEMPERATURE
	
			for(k =0; k<8; k++)
			{
				T = Temp[k];


				for(MCS = 0; MCS < Niter; MCS++)
				{	
					dummy_theta = 0.0; dummy_phi = 0.0;
		
					for(mm=0; mm < L; mm++)
					{
						for(nn=0; nn<6*L; nn++)
						{
							i = rand() % L ;  j = rand() % (6*L);
				
							if(j%6 ==0)											 // Defining the Nearest Neighbours with PBC
							{
								a1 = i; a2 = (j+1);                                               
								b1 = i; b2 = (j+3); 
								c1 = i; c2= (j+4); 
								d1 = i; d2 = (j + 6*L -1)%(6*L);   
								sign_a = -1.0;
								sign_b = 1.0;
								sign_c = 1.0;
								sign_d = -1.0; 
							} 
							else if(j%6 == 1)
							{
								a1 = i; a2 = (j+1);                                                
								b1 = i; b2 = (j+3); 
								c1 = i; c2 = (j+4); 
								d1 = i; d2 = (j-1);    
								sign_a = -1.0;
								sign_b = -1.0;
								sign_c = 1.0;
								sign_d = 1.0; 
							}
							else if(j%6 == 2)
							{
								a1 =i; a2 = (j+1);                                                
								b1 = i; b2 = (j + 3); 
								c1 = (i+1)%L; c2 = (j+2);    
								d1 = i; d2 = (j-1); 
								sign_a = -1.0;
								sign_b = -1.0;
								sign_c = 1.0;
								sign_d = 1.0;    
		   					}
		   					else if(j%6 == 3)
							{
								a1 = i; a2 = (j-1);                                                
								b1 = i; b2 = (j-3); 
								c1 = (i+1)%L; c2 = (j+1);
								d1 = i; d2 = (j+6*L -4)%(6*L);
								sign_a = 1.0;
								sign_b = -1.0;
								sign_c = -1.0;
								sign_d = 1.0;    
		   					}
		   					else if(j%6 == 4)
							{
								a1 = i; a2 = (j-3);                                                
								b1 = i; b2 = (j-4);
								c1 = (i-1+L)%L; c2 = (j-1);
								d1 = (i-1+L)%L; d2 = (j-2);
								sign_a = 1.0;
								sign_b = -1.0;
								sign_c = 1.0;
								sign_d = -1.0;    
		   					}
		   					else if(j%6 == 5)
							{
								a1 = i; a2 = (j-3) ;                                          
								b1 = i; b2 = (j-4);
								c1 = i; c2 = (j+1)%(6*L);
								d1 = i; d2 = (j+4)%(6*L);
								sign_a = 1.0;
								sign_b = -1.0;
								sign_c = 1.0;
								sign_d = -1.0;    
		   					}
				
							Ei =    - J_ising*(     Hxxz(theta[i][j], phi[i][j], theta[a1][a2], phi[a1][a2], delta) 
												  + Hxxz(theta[i][j], phi[i][j], theta[b1][b2], phi[b1][b2], delta) 
												  + Hxxz(theta[i][j], phi[i][j], theta[c1][c2], phi[c1][c2], delta) 
												  + Hxxz(theta[i][j], phi[i][j], theta[d1][d2], phi[d1][d2], delta) )
					  							  + DD*( sign_a*Hdm(theta[i][j], phi[i][j], theta[a1][a2], phi[a1][a2]) 
						 				           	  + sign_b*Hdm(theta[i][j], phi[i][j], theta[b1][b2], phi[b1][b2]) 
										              + sign_c*Hdm(theta[i][j], phi[i][j], theta[c1][c2], phi[c1][c2]) 
										              + sign_d*Hdm(theta[i][j], phi[i][j], theta[d1][d2], phi[d1][d2]) )  ; 
		
							dummy_theta = theta[i][j]; dummy_phi = phi[i][j];
		
			
							r = (float) rand() / RAND_MAX;
							theta[i][j] = theta[i][j] + (r - 0.5)*dtheta;
							if(theta[i][j] < 0.0) theta[i][j] = (float) abs(theta[i][j]) ;
							theta[i][j] = fmod(theta[i][j] + 2.0*PI, PI);
			
							r = (float) rand() / RAND_MAX;
							phi[i][j] = phi[i][j] + (r - 0.5)*dphi;
							if(phi[i][j] < 0.0) phi[i][j] = (float) abs(phi[i][j]) ;
							phi[i][j] = fmod(phi[i][j] + 2.0*PI, 2.0*PI);
							//phi[0][0] = 0.0;
			
			
							Ef =    - J_ising*(     Hxxz(theta[i][j], phi[i][j], theta[a1][a2], phi[a1][a2], delta) 
												  + Hxxz(theta[i][j], phi[i][j], theta[b1][b2], phi[b1][b2], delta) 
												  + Hxxz(theta[i][j], phi[i][j], theta[c1][c2], phi[c1][c2], delta) 
												  + Hxxz(theta[i][j], phi[i][j], theta[d1][d2], phi[d1][d2], delta) )
								  				  + DD*( sign_a*Hdm(theta[i][j], phi[i][j], theta[a1][a2], phi[a1][a2]) 
												  + sign_b*Hdm(theta[i][j], phi[i][j], theta[b1][b2], phi[b1][b2]) 
												  + sign_c*Hdm(theta[i][j], phi[i][j], theta[c1][c2], phi[c1][c2]) 
												  + sign_d*Hdm(theta[i][j], phi[i][j], theta[d1][d2], phi[d1][d2]) )  ; 
						
							dE = Ef - Ei;
							
							if(dE < 0.0)
							{
								E = E + dE;
								Mx = Mx + Sx(theta[i][j], phi[i][j]) - Sx(dummy_theta, dummy_phi);
								My = My + Sy(theta[i][j], phi[i][j]) - Sy(dummy_theta, dummy_phi);
								Mz = Mz + Sz(theta[i][j], phi[i][j]) - Sz(dummy_theta, dummy_phi);
				
							}
							else
							{
								r = (float) rand() / RAND_MAX;
								u = exp(-dE/T);
					
								if(r < u)
								{
									E = E + dE;
									Mx = Mx + Sx(theta[i][j], phi[i][j]) - Sx(dummy_theta, dummy_phi);
									My = My + Sy(theta[i][j], phi[i][j]) - Sy(dummy_theta, dummy_phi);
									Mz = Mz + Sz(theta[i][j], phi[i][j]) - Sz(dummy_theta, dummy_phi);
				
								}
								else
								{
									theta[i][j] = dummy_theta; phi[i][j] = dummy_phi;
								}       
					
							}    
				
						} 
		
					}
					
		
					if(k==7 && MCS >= 65000)
					{
			
						for(i=0; i<L; i++)
						{
							for(j=0; j<6*L; j++)
							{
								av_theta[i][j] =  av_theta[i][j] + theta[i][j];
								av_phi[i][j] = av_phi[i][j] +  phi[i][j];
							}
						
						}
					
						count = count +1;
					}
		
				}

			} // temp bracket
	
	
			for(i=0; i<L; i++)
			{
				for(j=0; j<6*L; j++)
				{	
					fprintf(file01, "%f \t %f \t%f \t %f \n", delta, DD,  av_theta[i][j]/(float)count,  av_phi[i][j]/(float)count);
				}	
			}   
							
   		}   //parameter grid brackect
   	} // parameter grid bracket 
    	
	
	fclose(file01);
	
	
   	return(0);
}

float Sx( float theta, float phi )
{
	return sin(theta)*cos(phi);
}

float Sy(float theta, float phi)
{
	return sin(theta)*sin(phi);
}

float Sz(float theta, float phi)
{
	return cos(theta);
}

float Hxxz( float theta1, float phi1, float theta2, float phi2, float Delta )
{
	return Sx(theta1, phi1)*Sx(theta2, phi2) + Sy(theta1, phi1)*Sy(theta2, phi2) + (Delta)*Sz(theta1, phi1)*Sz(theta2, phi2);
}

float Hdm(float theta1, float phi1, float theta2, float phi2)
{
	return ( Sx(theta1, phi1)*Sy(theta2, phi2) - Sy(theta1, phi1)*Sx(theta2, phi2) ) ;
}





