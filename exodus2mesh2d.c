/*
This program converts the ASCII exodus file exported from the CUBIT to 
several input files required by the SPECFEM2D program. Currently, this 
program only handles the 2D quadrilateral elements with four nodes.
The binary exodus file (e.g., .e file) needs to be converted into ASCII file,
generally using a free console application "ncdump" which is a part of the
netCDF library, and can be downloaded from 
http://www.unidata.ucar.edu/downloads/netcdf/index.jsp. This is already
installed on clover. The basic steps starting from the CUBIT:
step1: export mesh file as exodus file say "mesh.e". 
  If the element types are SHELL or SHELL4, "Default" or 3D option should
  be selected during export. If the element type is QUAD or QUAD4, 3D 
  option should be selected. With default or 2D data, it saves only
  X and Y coordinates which is not always correct.
  Make sure that the node ordering is strictly anticlockwise
  (no longer necessary!) for all the elements in CUBIT.

step2: convert "mesh.e" to ASCII file using
  ncdump mesh.e > mesh.txt

step3: compile exodus2mesh2d.c if not already compiled (see below).

step4: produce SPECFEM2D files using
  exodus2mesh2d mesh.txt
  OR
  exodus2mesh2d mesh.txt 1000.0 
There will be several files, file names starting with mesh_. 
------------------------------------------------------------------------------
DEPENDENCY:
  stringmanip.c
COMPILE:
  gcc exodus2mesh2d.c -o exodus2mesh2d -lm
USAGE: 
  exodus2mesh2d <inputfile> OPTIONS
  Example: exodus2mesh2d mest.txt
  Example: exodus2mesh2d mest.txt -check_order=1
OPTIONS:
  -fac: use this option to multiply coordinates. This is important for unit 
    conversion, e.g., to convert m to km use -fac=0.001.
  -check_order: use this option to check node ordering, and convert clockwise
    ordering to counterclockwise ordering as required by specfem2d.
    0: NO, 1: YES (Default). 
TODO:
  - make efficient similar to exodus2mesh.c or better to integrate both 2D
    and 3D. For that, it would be very easy if we had similar formats for
    SPECFEM2D and SPECFEM3D
  - update to be applicable for any elements in 2D at least all TRI, QUAD
    elements!
NOTES:
  - CUBIT by default saves the 2D elements as SHELL elements which are 3D, in
    principle. Therefore, we need to specify QUAD4 in block element type during
    export to EXODUS file.		
REVISION:
  - HNG, May 05,2014
  - HNG, Sep 29,2010 
  - HNG, Feb 08,2009
----------------------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "stringmanip.c"

#define maxnod 8
#define ON 1
#define OFF 0
#define TEST_JACOBIAN 0 /* 0 and 1; developer version */
/* void removeExtension(char *, char *); */
/* auxiliary routines */
void  removeExtension(char *, char *);
int   look_float(float *, char *, char *);
int   shape(float,float,float**);
int   main(int argc,char **argv)
{
int   e,i,j,k,ndim;
int   i1,i2,nod1,nod2,n1,n2,n3,n4,enode[4];
int   ielmt,nnode,nelmt;
int   nblk,nss;
/* int node_count; *//* element, node count */
int   blk_count,ns_count; /* block, node set count */
int   nnode_id,nelmt_id,nblk_id,mat_stat,con_stat,nss_id,ss_id;
int   fmat_id,fcon_id,fcoor_id;
int   blk_att[10],**blk_enode;
int   ibelmt,ibside,nbelmt,nbside;
int   *ss_elmt,*ss_side;
int   *blk_nelmt,*blk_nenod,blk_natt[100],ss_nelmt[100],ss_ndf[100];
/* ss_nelmt=ss_nside */
int   side1,side2,side3,side4;
int   switch_coord[3];
float **coord,*xp,*zp;
float x[4],z[4],s[4],t[4],**lag4;
float dx_ds,dx_dt,dz_ds,dz_dt,detJ;
int   itmp,check_order;
float fac,ftmp; /* multiplication factor for coordinates, temporary float */
char  token[62],dumc[62],etype[12],stag[62];
char  fonly[62],outfname[62];
char  *blk_etype[100];
char  **ss_name; /* side set names */
FILE  *inf,*outf_bnd,*outf_mat,*outf_con,*outf_coor;
int   norder[4]={0,3,2,1};
int   inode,isdone,isflag,ncount,ncw,nccw;

printf("--------------------------------\n");

check_order=1; fac=1.0;
ielmt=ibelmt=ibside=0;
/* Initialize IDs */
nnode_id=0; nelmt_id=0, nblk_id=0, mat_stat=0, con_stat=0,nss_id=0,ss_id=0; /* OFF */
fmat_id=fcon_id=fcoor_id=0;
nbelmt=nblk=nbside=0;
blk_count=0;
for(i=0;i<3;i++)switch_coord[i]=ON;

/* Open input file */
if(argc<2){
  fprintf(stderr,"ERROR: input file not entered!\n");
  exit(-1);
}

/* scan command */
if(argc>2){
  for(i=2;i<argc;i++){
    if(look_float(&ftmp,"-fac=",argv[i])==0){
      fac=ftmp;		  
      continue;
    }else if(look_int(&itmp,"-check_order=",argv[i])==0){
      check_order=itmp;
      if(check_order!=0 && check_order!=1){
        fprintf(stderr,"ERROR: invalid -check_order value! Valid values are 0: NO and 1: YES!\n");
        exit(-1);
      }
      continue;
    }else{
      fprintf(stderr,"ERROR: unrecognized option \"%s\"\n",argv[i]);
      exit(-1);
    }
  }
}

inf=fopen(argv[1],"r");
if(inf==NULL){
  fprintf(stderr,"ERROR: file \"%s\" not found!\n",argv[1]);
  exit(-1);
}

removeExtension(argv[1],fonly);

while(!feof(inf)){ /* Scan token */
  fscanf(inf,"%s",token);

  /* Number of dimensions */
  if(!strcmp(token,"num_dim")){			
    fscanf(inf,"%s",dumc);
    fscanf(inf,"%d",&ndim);			
  }
	
  /* Number of nodes */
  if(!strcmp(token,"num_nodes")){			
    fscanf(inf,"%s",dumc);
    fscanf(inf,"%d",&nnode);
    nnode_id=1;
  }
	
  /* Number of elements */
  if(!strcmp(token,"num_elem")){			
    fscanf(inf,"%s",dumc);
    fscanf(inf,"%d",&nelmt);
    nelmt_id=1;
		
    /* Allocate array for connectivity */
    blk_enode = (int **)malloc(nelmt * sizeof(int *));
    for(i = 0; i < nelmt; i++)
      blk_enode[i] = (int *)malloc(maxnod * sizeof(int));
    }
	
    /* Number of blocks */		
    if(!strcmp(token,"num_el_blk")){
      fscanf(inf,"%s",dumc);
      fscanf(inf,"%d",&nblk);
      nblk_id=1;

      /* allocate memory */
      blk_nelmt=malloc(nblk*sizeof(int));
      blk_nenod=malloc(nblk*sizeof(int));
    }
	
    /* Number of side sets */		
    if(!strcmp(token,"num_side_sets")){
      fscanf(inf,"%s",dumc);
      fscanf(inf,"%d",&nss);
      nss_id=1;
      if(nss>0){
        /* allocate memory */
        ss_name=malloc(nss*sizeof(char *));
        for(i=0;i<nss;i++){
          ss_name[i]=malloc(62*sizeof(char));
          /* each name has maximum of 62 characters */
        }
      }
    }
    if(nss_id==1){
      /* This segment has a significance only if nss has legitimate value */
      for(i=0;i<nss;i++){
        sprintf(stag,"num_side_ss%d",i+1);
        if(!strcmp(token,stag)){
          fscanf(inf,"%s",dumc);
          fscanf(inf,"%d",&ss_nelmt[i]);
          nbelmt+=ss_nelmt[i];
          nbside+=ss_nelmt[i];
        }				
		
        sprintf(stag,"num_df_ss%d",i+1);
        if(!strcmp(token,stag)){
          fscanf(inf,"%s",dumc);
          fscanf(inf,"%d",&ss_ndf[i]);
          if(i==nss-1){
            ss_id=1; /* All side sets pre-information available */
            ss_elmt = (int *)malloc(nbelmt * sizeof(int));
            ss_side = (int *)malloc(nbside * sizeof(int));
          }
        }
      }
    }
	
    /* Elements type */		
    if(!strcmp(token,"connect1:elem_type")){
      fscanf(inf,"%s",dumc);
      fscanf(inf,"%s",etype);
      printf("element type: %s\n",etype);
    }
    /* Read boundary conditions */
    if(ss_id==1){
      /* read coordinate names */
      if(!strcmp(token,"ss_names")){    
        fscanf(inf,"%s",dumc); /* = */      
        for (i=0; i<nss; i++){
          fscanf(inf,"%s",dumc);    
          getfirstquote(dumc,ss_name[i]);
          /*printf("ss name: %s\n",ss_name[i]);*/
        } 
        /*continue;*/
      }/* read coordinate names */			  	
      for(i=0;i<nss;i++){
        sprintf(stag,"elem_ss%d",i+1);
        if(!strcmp(token,stag)){
          fscanf(inf,"%s",dumc);
          for(j=0;j<ss_nelmt[i];j++){
            fscanf(inf,"%d,",&ss_elmt[ibelmt]);
            ibelmt++;
          }
        }				
			
        sprintf(stag,"side_ss%d",i+1);
        if(!strcmp(token,stag)){
          fscanf(inf,"%s",dumc);
          for(j=0;j<ss_nelmt[i];j++){
            fscanf(inf,"%d,",&ss_side[ibside]);
            ibside++;
          }
        }
      }
      if(!strcmp(token,"elem_map")){
        fscanf(inf,"%s",dumc);
        for(i=0;i<nelmt;i++){
          fscanf(inf,"%d,",&itmp);
          /* printf("%d %d\n",i,itmp); */
          if(i!=itmp-1){
            printf("%d %d\n",i,itmp);
            exit(-1);
          }
        }
      }
    }
	
    if(nblk_id==1){
      /* This segment has a significance only if nblk has legitimate value */
      for(i=0;i<nblk;i++){
        sprintf(stag,"num_el_in_blk%d",i+1);
        if(!strcmp(token,stag)){
          fscanf(inf,"%s",dumc);
          fscanf(inf,"%d",&blk_nelmt[i]);
        }				
			
        sprintf(stag,"num_nod_per_el%d",i+1);
        if(!strcmp(token,stag)){
          fscanf(inf,"%s",dumc);
          fscanf(inf,"%d",&blk_nenod[i]);
        }
			
        sprintf(stag,"num_att_in_blk%d",i+1);
        if(!strcmp(token,stag)){
          fscanf(inf,"%s",dumc);
          fscanf(inf,"%d",&blk_natt[i]);
        }
      }				
    }
	
    /* Attributes */	
    /*
    if(nblk_id==1 && !mat_stat){
      if(!fmat_id){
        sprintf(outfname,"%s_materials",fonly);
        outf_mat=fopen(outfname,"w");
        fmat_id=1;
      }
		
      for(i=0;i<nblk;i++){				
        sprintf(stag,"attrib%d",i+1);
        if(!strcmp(token,stag)){
          if(i==nblk-1)mat_stat=1;
          fscanf(inf,"%s",dumc);
          for(j=0;j<blk_nelmt[i];j++){
            for(k=0;k<blk_natt[i];k++){
            fscanf(inf,"%d,",&blk_att[k]);					
          }
          fprintf(outf_mat,"%d\n",blk_att[blk_natt[i]-1]); // Print only the last attribute
        }
      }
    }
  }
  */
	
  /* Connectivity */
  /*
  if(nblk_id==1 && !con_stat){
    if(!fcon_id){
      sprintf(outfname,"%s_connectivity",fonly);
      outf_con=fopen(outfname,"w");
      fcon_id=1;
      fprintf(outf_con,"%d\n",nelmt);
    }
			
    for(i=0;i<nblk;i++){
      sprintf(stag,"connect%d",i+1);
      if(!strcmp(token,stag)){
        if(i==nblk-1)con_stat=1;
        fscanf(inf,"%s",dumc);
        for(j=0;j<blk_nelmt[i];j++){
          for(k=0;k<blk_nenod[i];k++){
            fscanf(inf,"%d,",&blk_enode[ielmt][k]);
            //fprintf(outf_con,"%d ",blk_enode[ielmt][k]);
          }
          fprintf(outf_con,"%d %d %d %d\n",blk_enode[ielmt][0],blk_enode[ielmt][3],blk_enode[ielmt][2],blk_enode[ielmt][1]);
          //fprintf(outf_con,"%d %d %d %d\n",blk_enode[ielmt][3],blk_enode[ielmt][0],blk_enode[ielmt][1],blk_enode[ielmt][2]);
          ielmt++;
          //fprintf(outf_con,"\n"); // New token 
        }
      }
    }
  }
  */

  /* Connectivity */
  if(nblk>0 && con_stat!=ON){   
    /* write connectivity and material id */ 		
    for(i=0;i<nblk;i++){
      sprintf(stag,"connect%d",i+1);
      if(strcmp(token,stag)==0){
        blk_count++;
	  
        /* open connectivity and material files */
        if(blk_count==1){          
          printf("storing connectivity and saving materials..."); 
          /* sprintf(outfname,"%s_connectivity_noncorrected",fonly); 
          outf_con=fopen(outfname,"w");          
          fprintf(outf_con,"%d\n",nelmt); */
	      
          sprintf(outfname,"%s_materials",fonly);
          outf_mat=fopen(outfname,"w");				
          /*fprintf(outf_mat,"%d\n",nelmt);*/
        }
	    
        fscanf(inf,"%s",dumc); /* = */
        for(j=0;j<blk_nelmt[i];j++){
          for(k=0;k<blk_nenod[i];k++){
            fscanf(inf,"%d,",&itmp);
            blk_enode[ielmt][k]=itmp;
            /* fprintf(outf_con,"%d ",itmp); */		
          }
          /* fprintf(outf_con,"%d %d %d %d\n",blk_enode[ielmt][0],
          blk_enode[ielmt][3],blk_enode[ielmt][2],blk_enode[ielmt][1]);*/
          ielmt++;
          /* fprintf(outf_con,"\n"); *//* new token */
          fprintf(outf_mat,"%d\n",i+1);
        } 
        /* printf("i: %d nblk: %d\n",i,blk_count);
        exit(-1); */
        if(blk_count==nblk){			  
          con_stat=ON;
          mat_stat=ON;
          /*free(blk_nelmt);
          free(blk_nenod);*/
          printf("complete!\n");
          /* fclose(outf_con); */
          fclose(outf_mat);			  
        }        
        continue;
      }        
    }
  }

  /* Coordinates */
  if(!strcmp(token,"coord")){
    if(ndim==3)switch_coord[1]=OFF;
    /* for SHELL element make Y coordinates OFF. TODO: add option for this. */
    fscanf(inf,"%s",dumc);
    coord = (float **)malloc(ndim * sizeof(float *));
    for(i = 0; i < ndim; i++)coord[i] = (float *)malloc(nnode * sizeof(float));
		
    if(!fcoor_id){
      sprintf(outfname,"%s_coordinates",fonly);
      outf_coor=fopen(outfname,"w");
      fcoor_id=1;
    }
		
    fprintf(outf_coor,"%d\n",nnode);
		
    for(i=0;i<ndim;i++){
      for(j=0;j<nnode;j++){
        fscanf(inf,"%f,",&coord[i][j]); /* Comma separated data */
        /* printf("%f\n",coord[i][j]); */
      }
    }
		
    for(j=0;j<nnode;j++){
      for(i=0;i<ndim;i++){
        /* printf("%d\n",switch_coord[1]);exit(-1); */              
        if(switch_coord[i]!=OFF)fprintf(outf_coor,"%.6f ",fac*coord[i][j]); /* Do not write Y coordinates */
          /* printf("%f %f %f\n",coord[i][j],fac*coord[i][j],fac); */
      }
      fprintf(outf_coor,"\n");
    }
    fclose(outf_coor);			
  } /* (!strcmp(token,"coord")) */
}/* while(!feof(inf)) */
/*fclose(outf_mat);
fclose(outf_con);
free(coord); */

/* Connectivity */
if(nblk>0 && con_stat==ON){   
  xp=(float *)malloc(4 * sizeof(float));
  zp=(float *)malloc(4 * sizeof(float));

  if(check_order==0){
    printf("connectivity order: preserve!\n");
  }else if(check_order==1){
    printf("connectivity order: check for clockwise ordering!\n");
  }
  printf("saving connectivity..."); 
  sprintf(outfname,"%s_connectivity",fonly);
  outf_con=fopen(outfname,"w");          
  fprintf(outf_con,"%d\n",nelmt);		
			  
  /* WARNING: only for 4 noded elements */
  ncount=0;  
  ncw=0; nccw=0;                   
  for(i=0; i<nelmt; i++){
    /*printf("Before coordX: %6.2f %6.2f %6.2f %6.2f\n",xp[0],xp[1],xp[2],xp[3]);
    printf("Before coordZ: %6.2f %6.2f %6.2f %6.2f\n",zp[0],zp[1],zp[2],zp[3]);*/
    n1=blk_enode[i][0]; n2=blk_enode[i][1];
    n3=blk_enode[i][2]; n4=blk_enode[i][3];
    if(check_order==1){
      /* check ordering, and if found clockwise convert them to counterclockwise */
      for(j=0; j<4; j++){
        inode=blk_enode[i][j]-1;
        xp[j]=coord[0][inode];
        zp[j]=coord[2][inode];
      }
      isflag=isclockwise(4,xp,zp);
      if(isflag>0){
        /* found clockwise */
        ncw+=1;
        /* convert to counterclockwise */
        n1=blk_enode[i][0]; n2=blk_enode[i][3]; n3=blk_enode[i][2]; n4=blk_enode[i][1];
      }
      if(isflag<0)nccw+=1;         
    }
    /* preserve the order */
    fprintf(outf_con,"%d %d %d %d\n",n1,n2,n3,n4);  
  }		
			  
  con_stat=ON;
  printf("complete!\n");
  fclose(outf_con);  
  free(xp);
  free(zp);
  free(blk_nelmt);
  free(blk_nenod);
}
printf("clockwise-ordered elements: %d\n",ncw);
printf("anticlockwise-ordered elements: %d\n",nccw);

/* Write boundary files */
/* Side numbering: 
for QUAD elements ndim=2 or 3, and side numbering - 1,2,3,4 but 
for SHELL element ndim=3, and side numbering - 3,4,5,6 */
side1=1; side2=2; side3=3; side4=4;
/* side1=1; side2=4; side3=3; side4=2; */
if(ndim==3 & strstr(etype,"SHELL")!=NULL){
  /* SHELL elements. Need to check carefully whether true for all SHELL
     elements, and how they number. */
  side1+=2; side2+=2; side3+=2; side4+=2;
}
if(ss_id==1){
  i1=0; i2=0;
  for(i=0;i<nss;i++){                        
    sprintf(outfname,"%s_%s",fonly,ss_name[i]); /*,i+1);*/
    outf_bnd=fopen(outfname,"w");
    fprintf(outf_bnd,"%d\n",ss_nelmt[i]);
    i2+=ss_nelmt[i];
    for(j=i1;j<i2;j++){
      if(ss_side[j]==side1){
        nod1=1; nod2=2;
      }else if(ss_side[j]==side2){
        nod1=2; nod2=3;
      }else if(ss_side[j]==side3){
        nod1=3; nod2=4;
      }else if(ss_side[j]==side4){
        nod1=4; nod2=1;
      }else{
        fprintf(stderr,"ERROR: wrong side ID: %d for boundary!\n",ss_side[j]);
        exit(-1);
      }
			
      e=ss_elmt[j];
      if(strstr(ss_name[i],"absorbing")!=NULL){ 
        fprintf(outf_bnd,"%d %d %d %d %d\n",e,2,blk_enode[e-1][nod1-1],
        blk_enode[e-1][nod2-1],ss_side[j]); /* This is only for edge */
      }else{
        fprintf(outf_bnd,"%d %d %d %d\n",e,2,blk_enode[e-1][nod1-1],
        blk_enode[e-1][nod2-1]); /* This is only for edge */
      }
    }
    fclose(outf_bnd);
    i1+=ss_nelmt[i];
  }
}else{
  printf("WARNING: no boundary files written!\n");
}

if(TEST_JACOBIAN==1){
/* Check for negative jacobian */
lag4 = (float **)malloc(3 * sizeof(float *));
for(i = 0; i < 3; i++)lag4[i] = (float *)malloc(4 * sizeof(float));

s[0]=-1.0; t[0]=-1.0;
s[1]= 1.0; t[1]=-1.0;
s[2]= 1.0; t[2]= 1.0;
s[3]=-1.0; t[3]= 1.0;
ncount=0;
for(i=0; i<nelmt; i++){
  n1=blk_enode[i][0]; n2=blk_enode[i][1];
  n3=blk_enode[i][2]; n4=blk_enode[i][3];
  isdone=0;
  for(j=0; j<4; j++){
    shape(s[j],t[j],lag4);
    printf("%f %f %f %f\n",lag4[0][0],lag4[0][1],lag4[0][2],lag4[0][3]);
    printf("%f %f %f %f\n",lag4[1][0],lag4[1][1],lag4[1][2],lag4[1][3]);
    printf("%f %f %f %f\n",lag4[2][0],lag4[2][1],lag4[2][2],lag4[2][3]);
		
    /* Anticlockwise mapping */			
    x[0]=coord[0][n1]; x[1]=coord[0][n2];
    x[2]=coord[0][n3]; x[3]=coord[0][n4];
    if(ndim==2){
      z[0]=coord[1][n1]; z[1]=coord[1][n2];
      z[2]=coord[1][n3]; z[3]=coord[1][n4];
    }else{
      z[0]=coord[2][n1]; z[1]=coord[2][n2];
      z[2]=coord[2][n3]; z[3]=coord[2][n4];
    }
			
    dx_ds=0.0; dx_dt=0.0; dz_ds=0.0; dz_dt=0.0;
    for(k=0; k<4; k++){
      dx_ds+=x[k]*lag4[1][k]; dx_dt+=x[k]*lag4[2][k];
      dz_ds+=z[k]*lag4[1][k]; dz_dt+=z[k]*lag4[2][k];
    }
    detJ=dx_ds*dz_dt-dx_dt*dz_ds;
    printf("Before -element: %d j: %d %f\n",i+1,j,detJ);

    if(!isdone && detJ<=0.0){
      ncount+=1;
      isdone=1;
    }
		
    if(detJ <= 0.0){
      /* printf("Negative Jacobian: %f for element: %d\n",detJ,i+1); */
      /* Clockwise mapping */
      x[0]=coord[0][n1]; x[1]=coord[0][n4];
      x[2]=coord[0][n3]; x[3]=coord[0][n2];
      if(ndim==2){
        z[0]=coord[1][n1]; z[1]=coord[1][n2];
        z[2]=coord[1][n3]; z[3]=coord[1][n4];
      }else{
        z[0]=coord[2][n1]; z[1]=coord[2][n2];
        z[2]=coord[2][n3]; z[3]=coord[2][n4];
      }
			
      dx_ds=0.0; dx_dt=0.0; dz_ds=0.0; dz_dt=0.0;
			
      for(k=0; k<4; k++){
        dx_ds+=x[k]*lag4[1][norder[k]]; dx_dt+=x[k]*lag4[2][norder[k]];
        dz_ds+=z[k]*lag4[1][norder[k]]; dz_dt+=z[k]*lag4[2][norder[k]];
      }
      detJ=dx_ds*dz_dt-dx_dt*dz_ds;
      printf("After - element: %d j: %d %f\n",i+1,j,detJ);
      if(detJ <= 0.0){
        printf("After - element: %d j: %d %f\n",i+1,j,detJ);
      }
    }
  }
  /* exit(-1); */
}

printf("ncount:%d",ncount);
free(lag4);	
free(blk_enode);
free(ss_elmt);
free(ss_side);
}
		
printf("blocks: %d, elements: %d, nodes: %d\n",nblk,nelmt,nnode);
printf("successfully finished!\n");
printf("--------------------------------\n");
return(0);
}
/*----------------------------------------------------------------------------*/

/* shape functions */
#define quart 0.25
int shape(float s, float t, float **lag4)
{
int i,j;	
float sm,sp,tm,tp;

printf("In function: %f %f\n",s,t);
sp = s + 1.0; sm = s - 1.0;
tp = t + 1.0; tm = t - 1.0;

/*lag4 = (float **)malloc(3 * sizeof(float *));
for(i = 0; i < 3; i++)lag4[i] = (float *)malloc(4 * sizeof(float));*/

/* Shape functions */
lag4[0][0]=quart*sm*tm; lag4[0][1]=-quart*sp*tm; lag4[0][2]=quart*sp*tp; lag4[0][3]=-quart*sm*tp;

/* Derivatives with respect to s */
lag4[1][0]=quart*tm; lag4[1][1]=-quart*tm; lag4[1][2]=quart*tp; lag4[1][3]=-quart*tp;

/* Derivatives with respect to t */
lag4[2][0]=quart*sm; lag4[2][1]=-quart*sp; lag4[2][2]=quart*sp; lag4[2][3]=-quart*sm;

return 0;
}
/*----------------------------------------------------------------------------*/

/* check whether the nodes are clockwise */
#define SUMEDGE 0
#define NORMAL 1
#define APPROACH NORMAL /* SUMEDGE or NORMAL */
int isclockwise(int n, float x[n], float z[n])
{
int i,iflag;
float sum_edge;
float norm,nvec[3],v1[3],v2[3];

if(APPROACH == SUMEDGE){
  /* sum over edges approach */
  sum_edge=0.0;
  for(i=0; i<n-1; i++){
    sum_edge+=(x[i+1]-x[i])*(z[i+1]+z[i]);
  }
  sum_edge+=(x[0]-x[n-1])*(z[0]+z[n-1]);  
  if(sum_edge>0){
    iflag=1;
  }else{
    iflag=-1;
  }
}else{
  /* normal approach */
  v1[0]=  x[1]-x[0]; v1[1]=0.0; v1[2]=  z[1]-z[0];
  v2[0]=x[n-1]-x[0]; v2[1]=0.0; v2[2]=z[n-1]-z[0];

  nvec[0]=v1[1]*v2[2]-v1[2]*v2[1];
  nvec[1]=v1[2]*v2[0]-v1[0]*v2[2];
  nvec[2]=v1[0]*v2[1]-v1[1]*v2[0];
  norm=sqrt(nvec[0]*nvec[0]+nvec[1]*nvec[1]+nvec[2]*nvec[2]);
  nvec[0]/=norm;
  nvec[1]/=norm;
  nvec[2]/=norm;
  if(nvec[1]>0.0){
    iflag=1;
  }else if(nvec[1]<0.0){
    iflag=-1;
  }else{
    fprintf(stderr,"ERROR: Zero jacobian! %f %f\n",nvec[1],norm);
    fprintf(stderr,"coordX: %6.2f %6.2f %6.2f %6.2f\n",x[0],x[1],x[2],x[3]);
    fprintf(stderr,"coordZ: %6.2f %6.2f %6.2f %6.2f\n",z[0],z[1],z[2],z[3]);
    exit(-1);
  }
}
return(iflag);
}
/*----------------------------------------------------------------------------*/

