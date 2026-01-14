/**
   @Program : D2min.c
   @Author: Himangsu Bhaumik, 2018 
   @Note    : For details see PRE.57.7192(1998)
   @Compilation : gcc D2min.c -L/Data/.gsl-2.2.1/lib/ -lgsl -lgslcblas -lm

  

   execution: ./a.out N initframe finalframe

 
**/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_coupling.h>
#include<gsl/gsl_sf_log.h>



#define maxNN 64

#define e0str "10"
#define str "e0-10_S0-1_B-0.01"




void load_frame(int N,FILE **fp, int iframe, double *L,double *tild, double ***xyz, int **type){

  int i,id, itype,tstep; double t,Li,Lf,x_cor,y_cor,z_cor,vx,vy,vz;

  //FILE **fp=fopen(FILEname,"r");;
  
  rewind(*fp);
  int frame_read=0;
  //Skip frist (iframe-1) frame data
  for(i=0;i<iframe*(N+9);i++)fscanf(*fp,"%*[^\n]\n");
  while (frame_read<1){
    fscanf(*fp,"%*[^\n]\n");
    fscanf(*fp,"%d\n",&tstep);
    printf("time step=%d\n",tstep);
    for(i=0;i<3;i++)fscanf(*fp,"%*[^\n]\n");
    fscanf(*fp,"%*f %*f %lf\n",&t);*tild=t;
    fscanf(*fp,"%lf %lf %*f\n",&Li,&Lf);(*L)=Lf-Li;

    for(i=0;i<2;i++)fscanf(*fp,"%*[^\n]\n");
    for(i=0;i<N;i++){
      fscanf(*fp,"%d %d %lf %lf %lf %lf %lf %lf\n",&id,&itype,&x_cor,&y_cor,&z_cor,&vx,&vy,&vz);
      id-=1;
      (*xyz)[id][0]=x_cor;
      (*xyz)[id][1]=y_cor;
      (*xyz)[id][2]=z_cor;
      (*type)[id]=itype;
    }
    frame_read++;
  }
  
  printf("loaded\n");


}

int delta(int i, int j){
  int value;
  if(i==j) value=1;
  else value=0;
  return value;
}

#include <gsl/gsl_linalg.h>
void inverse (double *a_data, double *b_data){
     
     double inva[9];
     int s, i, j;
     
     gsl_matrix_view m
          = gsl_matrix_view_array(a_data, 3, 3);
     gsl_matrix_view inv
          = gsl_matrix_view_array(inva,3,3);
     gsl_permutation * p = gsl_permutation_alloc (3);
 
     gsl_linalg_LU_decomp (&m.matrix, p, &s);    
     gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);
     for (i = 0; i < 3; ++i)
       for (j = 0; j < 3; ++j)
	 b_data[i*3+j]=gsl_matrix_get(&inv.matrix,i,j);
     gsl_permutation_free (p);
     //return 0;
     }



int main(int argc, char *argv[])     
{
  FILE *fp;
  static char FILE[1024];
  static double L,iL;
  static int N;  
  static long i,j,k;
  
  
  static double tild0,tild1;
  static double xr,yr,zr;
  static int initframe,finalframe,del_frame;
    
  static double dr;

      
 
  static double X[3][3],Y[3][3],Yinv[3][3],epsilon[3][3];
  static double a_data[9],b_data[9];
  static int ipart,ineighbour,ilist;
  static double dx1,dx2;
  static double D2,sumj,tot_neighbour,totD2;
  
    
  static float Rc=1.568;
  
 


  //--------------------------------------------------------------------
  
  //====================================================================
    
  
  if(argc>1){
    N=atof(argv[1]);
    initframe=atof(argv[2]);
    finalframe=atof(argv[3]);
    
    
  }
  else {
    printf("Give N and iframe as arguments\n");
    return 0;
  }
  
  del_frame=finalframe-initframe;
  if(del_frame<0) {printf("Should be: finalframe>initframe\n");return 0;}

  char command[128];
  sprintf(command, "mkdir -p %s",str);
  system(command);
  


  
    //--------------------------------------------------------------------
    //=================defining various array ============================
  double **xyz_0,**xyz_1;
      
    xyz_0=malloc(N*sizeof(double *));
    xyz_1=malloc(N*sizeof(double *));

    for(j=0;j<N;j++){
      xyz_0[j]=malloc(3*sizeof(double));
      xyz_1[j]=malloc(3*sizeof(double));

    }
    
    for(j=0;j<N;j++)
      for(k=0;k<3;k++){
	xyz_0[j][k]=0.0;
	xyz_1[j][k]=0.0;

      }
    
    
    int *type;
    type=malloc(N*sizeof(int));

    int **neighbourLIST;
    neighbourLIST=malloc(N*sizeof(int *));
    for(j=0;j<N;j++)neighbourLIST[j]=malloc(maxNN*sizeof(int));
    int *nlist;
    nlist=malloc(N*sizeof(int));

    int **neighbourLIST1;
    neighbourLIST1=malloc(N*sizeof(int *));
    for(j=0;j<N;j++)neighbourLIST1[j]=malloc(maxNN*sizeof(int));
    int *nlist1;
    nlist1=malloc(N*sizeof(int));



    
    double *D2Narray;
    D2Narray=malloc(N*sizeof(double));
    totD2=0;
    for(i=0;i<N;i++)D2Narray[i]=0;
    

    //-------------------------------------------------------------------
    


 
 
    
    //===================================================================

    
     sprintf(FILE,"../../../LAMMPS_script/data/set1/Data_full_e0-10_S0-1_B-0.01_steps-100000.dump");
    fp=fopen(FILE,"r");
    printf("file open %s\n",FILE);
    

    load_frame(N,&fp,initframe,&L,&tild0,&xyz_0,&type);
    
    load_frame(N,&fp,finalframe,&L,&tild1,&xyz_1,&type);

  
    
    fclose(fp);
    
    iL=1.0/L;
    
    //==========Transform co-ordinate to wrapped one========== 
      for(i=0;i<N;i++){
	for(j=0;j<3;j++){
	 while(xyz_0[i][j]<0)xyz_0[i][j]+=L;
	 while(xyz_0[i][j]>L)xyz_0[i][j]-=L;
	 while(xyz_1[i][j]<0)xyz_1[i][j]+=L;
	 while(xyz_1[i][j]>L)xyz_1[i][j]-=L;

	 
	 
	}
      }
      //========================================================= 


      //printf("N=%d L=%e tild0=%e tild1=%e\n",N,L,tild0,tild1);
     
      //========== Get the neighbour list of frame 0 ======================/
      
      for(i=0;i<N;i++)nlist[i]=0;
      
      for(ipart=0;ipart<N-1;ipart++){
	
	for(ineighbour=ipart+1;ineighbour<N;ineighbour++){

  	  xr=xyz_0[ipart][0]-xyz_0[ineighbour][0];
  	  yr=xyz_0[ipart][1]-xyz_0[ineighbour][1];
  	  zr=xyz_0[ipart][2]-xyz_0[ineighbour][2];

	  xr=xr-tild0*round(yr*iL);
  	  xr=xr-L*round(xr*iL);
  	  yr=yr-L*round(yr*iL);
  	  zr=zr-L*round(zr*iL);
  	  dr=sqrt(xr*xr+yr*yr+zr*zr);
	  
	  Rc=(0.88+(type[ipart]-1)*0.02+(type[ineighbour]-1)*0.02)*1.4;
	  if(dr<Rc){
	    neighbourLIST[ipart][nlist[ipart]++]=ineighbour;
	    neighbourLIST[ineighbour][nlist[ineighbour]++]=ipart;
	    if(nlist[ineighbour]>maxNN || nlist[ipart]> maxNN){
		printf("nlist exeeds maxNN, increase maxNN\n");
		return(0);
	      }
	  }
	}
      }
      //=========================================================       
      



        //========== Get the neighbour list of frame 1 ======================/
      // We getting this only for shortest path caluclation, no dr required.
      for(i=0;i<N;i++)nlist1[i]=0;
      
      for(ipart=0;ipart<N-1;ipart++){
	
	for(ineighbour=ipart+1;ineighbour<N;ineighbour++){

  	  xr=xyz_1[ipart][0]-xyz_1[ineighbour][0];
  	  yr=xyz_1[ipart][1]-xyz_1[ineighbour][1];
  	  zr=xyz_1[ipart][2]-xyz_1[ineighbour][2];

	  xr=xr-tild1*round(yr*iL);
  	  xr=xr-L*round(xr*iL);
  	  yr=yr-L*round(yr*iL);
  	  zr=zr-L*round(zr*iL);
  	  dr=sqrt(xr*xr+yr*yr+zr*zr);
	  
	  //if(dr<Rc[type[ipart]][type[ineighbour]]){
	  Rc=(0.88+(type[ipart]-1)*0.02+(type[ineighbour]-1)*0.02)*1.4;
	  if(dr<Rc){
	    
	    
	    
	    neighbourLIST1[ipart][nlist1[ipart]++]=ineighbour;
	    neighbourLIST1[ineighbour][nlist1[ineighbour]++]=ipart;
	    if(nlist1[ineighbour]>maxNN || nlist1[ipart]> maxNN){
		printf("nlist exeeds maxNN, increase maxNN\n");
		return(0);
	      }
	  }
	}
      }
      //========================================================= 






      //============= D2 calculation begin==============
      
      for(ipart=0;ipart<N;ipart++){
	if(nlist[ipart]>2){
	  
  	for(i=0;i<3;i++)
  	  for(j=0;j<3;j++){
  	    X[i][j]=0.0;Y[i][j]=0.0;Yinv[i][j]=0.0;
  	    epsilon[i][j]=0.0;
  	  }
	
	for(ilist=0;ilist<nlist[ipart];ilist++){
	  ineighbour=neighbourLIST[ipart][ilist];
	  
	  
	  for(i=0;i<3;i++)
	    for(j=0;j<3;j++){
	      
	      dx1=(xyz_1[ineighbour][i]-xyz_1[ipart][i]);
	      dx1=dx1-L*round(dx1*iL);
	      
	      if(i==0){ // considering tild factor
		xr=xyz_1[ineighbour][0]-xyz_1[ipart][0];
		//zr=xyz_1[ineighbour][2]-xyz_1[ipart][2];
		yr=xyz_1[ineighbour][1]-xyz_1[ipart][1];
		xr=xr-tild1*round(yr*iL);
		dx1=xr-L*round(xr*iL);
	      }
	      
	      
	      dx2=(xyz_0[ineighbour][j]-xyz_0[ipart][j]);
	      dx2=dx2-L*round(dx2*iL);
	      
	      if(j==0){ // considering tild factor
		xr=xyz_0[ineighbour][0]-xyz_0[ipart][0];
		//zr=xyz_0[ineighbour][2]-xyz_0[ipart][2];
		yr=xyz_0[ineighbour][1]-xyz_0[ipart][1];
		xr=xr-tild0*round(yr*iL);
		dx2=xr-L*round(xr*iL);
	      }

	      
	      
	      
	      X[i][j]+=dx1*dx2;
		
	      dx1=(xyz_0[ineighbour][i]-xyz_0[ipart][i]);
	      dx1=dx1-L*round(dx1*iL);
	      
	      if(i==0){ // considering tild factor
		xr=xyz_0[ineighbour][0]-xyz_0[ipart][0];
		//zr=xyz_0[ineighbour][2]-xyz_0[ipart][2];
		yr=xyz_0[ineighbour][1]-xyz_0[ipart][1];
		xr=xr-tild0*round(yr*iL);
		dx1=xr-L*round(xr*iL);
	      }
	      
	      
	      
	      Y[i][j]+=dx1*dx2;
	      
	    }
	  
	  
	}//ineighbour
	
	
	for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
	    a_data[i*3+j]=Y[i][j];
	
	inverse(a_data,b_data);
	
	for (i=0;i<9;++i)Yinv[i/3][i%3]=b_data[i];
	
	for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
	    for(k=0;k<3;k++)
	      epsilon[i][j]+=X[i][k]*Yinv[k][j];
	
	
	//Epsilon has been calculated till now for ipart particle
	
	
	D2=0.0;	  tot_neighbour=0.0;
	// For calculation of D2 again loop over neighbour
	
	//for(ineighbour=0;ineighbour<N;ineighbour++){
	for(ilist=0;ilist<nlist[ipart];ilist++){
	  
	  ineighbour=neighbourLIST[ipart][ilist];
	  
	  
	  
	  tot_neighbour++;
	  for(i=0;i<3;i++){
	    
	    sumj=0.0;
	    for(j=0;j<3;j++){
	      dx1=(xyz_0[ineighbour][j]-xyz_0[ipart][j]);
	      dx1=dx1-L*round(dx1*iL);
	      
	      if(j==0){ // considering tild factor
		xr=xyz_0[ineighbour][0]-xyz_0[ipart][0];
		//zr=xyz_0[ineighbour][2]-xyz_0[ipart][2];
		yr=xyz_0[ineighbour][1]-xyz_0[ipart][1];
		xr=xr-tild0*round(yr*iL);
		dx1=xr-L*round(xr*iL);
	      }
	      
	      sumj+=epsilon[i][j]*dx1;
	    }//j
	    
	    dx1=xyz_1[ineighbour][i]-xyz_1[ipart][i];
	    dx1=dx1-L*round(dx1*iL);
	    
	    if(i==0){ // considering tild factor
	      xr=xyz_1[ineighbour][0]-xyz_1[ipart][0];
	      //zr=xyz_1[ineighbour][2]-xyz_1[ipart][2];
	      yr=xyz_1[ineighbour][1]-xyz_1[ipart][1];
	      xr=xr-tild1*round(yr*iL);
	      dx1=xr-L*round(xr*iL);
	    }
	    
	    dx1-=sumj;
	    D2+=pow(dx1,2.0);
	    
	  }//i
	  
	  
	  
	}//ineighbour
	
	D2=D2/tot_neighbour;
	totD2+=D2;
	D2Narray[ipart]=D2;


	}//if neighour>1
      }//ipart
    

      

      //==========================================================
      
      
    
      
	sprintf(FILE,"%s/D2min_frame%d.dat",str,finalframe);
	fp=fopen(FILE,"w");
	for(i=0;i<N;i++){
	  fprintf(fp,"%f\n",D2Narray[i]);
	}
	fclose(fp);



     

    
    
    

 free(xyz_0);
 free(xyz_1);

    
    return 0;
} //end main



