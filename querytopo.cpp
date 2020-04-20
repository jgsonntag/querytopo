/*------------------------------------------------------------------------*
 NAME:     querytopo.cpp

 PURPOSE:  

 AUTHOR:   John Gary Sonntag

 DATE:     22 October 2014
 *------------------------------------------------------------------------*/

#include "/home/sonntag/Include/mission.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


main(int argc, char *argv[])
{
  char line[85],demid[10],wpname[10],geoidid[10];
  int htrefflag;
  double lat,lon,topo,geoid;
  double querytopo(double,double,char *);
  double querygeoid(double,double,char *);
  int initquerytopo(),closequerytopo();
  int initquerygeoid(),closequerygeoid();
  FILE *fptr;

  // Check input
  if (argc != 3)
  {
    printf("Usage: querytopo <latlon filename> <height ref (1=geoid 2=ellipsoid)>\n");
    exit(0);
  }

  // Open the input file
  if ((fptr=fopen(argv[1],"r"))==NULL)
  {
    printf("Input file %s not found - exiting\n",argv[1]);
    exit(-1);
  }
  htrefflag = atoi(argv[2]);
  if (htrefflag<1||htrefflag>2)
  {
    printf("Unrecognized height reference - exiting\n");
    exit(-1);
  }

  // Loop over the input file entries
  initquerygeoid();
  initquerytopo();
  while(fgets(line,85,fptr)!=NULL)
  {

    // Parse the input and ensure longitude is within bounds
    sscanf(line,"%lf %lf %s",&lat,&lon,wpname);
    while (lon<=-180.0) lon+=360.0;
    while (lon>180.0) lon-=360.0;

    // Query the topo and geoid databases
    topo = querytopo(lat,lon,demid);
    geoid = querygeoid(lat,lon,geoidid);

    // Reference the topo heights according to request and native reference of database
    if (!strncmp(demid,"GT3\0",3))
    {
      if (htrefflag==2) topo = topo+geoid;
      if (isnan(topo)) topo = 0.0;
    }
    else if (!strncmp(demid,"BM2\0",3))
    {
      if (htrefflag==2) topo = topo+geoid;
      if (isnan(topo)) topo = 0.0;
    }
    else if (!strncmp(demid,"G90\0",3))
    {
      if (htrefflag==1) topo = topo-geoid;
      if (isnan(topo)) topo = 0.0;
    }

    // Output the result
    printf("%lf %lf %s %lf %s %lf %s\n",lat,lon,wpname,topo,demid,geoid,geoidid);

  }

  // Close the input file
  closequerygeoid();
  closequerytopo();
  fclose(fptr);

//  //  Check input
//  if ( argc != 3 )
//  {
//    printf("Usage:  querytopo <lat (deg)> <lon (deg)>\n");
//    exit(-1);
//  }
//  lat = atof(argv[1]);
//  lon = atof(argv[2]);
//  if (lat>90.0||lat<-90.0)
//  {
//    printf("Input latitude invalid - exiting\n");
//    exit(-1);
//  }
//  while (lon<=-180.0) lon+=360.0;
//  while (lon>180.0) lon-=360.0;
//
//  // Get the topography
//  topo = querytopo(lat,lon);
//
//  // Output the result
//  printf("Result is: %lf m  %lf ft\n",topo,topo*3.28);

//  initquerytopo();
//  while (gets(line)!=NULL)
//  {
//    sscanf(line,"%lf %lf",&lat,&lon);
//    if (lat>90.0||lat<-90.0)
//    {
//      printf("Input latitude invalid - exiting\n");
//      exit(-1);
//    }
//    while (lon<=-180.0) lon+=360.0;
//    while (lon>180.0) lon-=360.0;
//    topo = querytopo(lat,lon,demid);
//    printf("%lf m  %lf ft %s\n",topo,topo*3.28,demid);
//
//  }
//
//  closequerytopo();

}


#define EGM96PATH "/usr/local/share/geoid/egm96/WW15MGH.DAC\0"
#define GTOPO30PATH "/usr/local/share/dem/gtopo30/\0"
#define BEDMAP2PATH "/usr/local/share/dem/bedmap2/bedmap2_surface.flt\0"
#define GIMP90PATH "/usr/local/share/dem/gimp90/gimp90m.dem\0"
#define C 299792458.0
#define WE 7.292115147e-5
#define MU 3.986005005e14
#define PI (4.0*atan((double)(1.0)))
#define AE 6378137.0
#define FLAT (1.0/298.257223563)
#define RAD2NM (180.0*60.0/PI)
#define RAD2KM (RAD2NM*6076.1*12.0*2.54/100.0/1000.0)

FILE *fptrdem;
FILE *fptrgeoid;
int opendem; // -1 nothing open
             // 0-32 corresponding GTOPO30 tile file open
             // 33 GIMP90 file open
             // 34 Bedmap-2 file open
int opengeoid; // -1 nothing open
             // 0 EGM96 open


int initquerytopo()
{
  opendem = -1; // No DEM files are open at start
}


int closequerytopo()
{
  //printf("closing previously open DEM\n");
  if (opendem!=-1) fclose(fptrdem);
  opendem=-1;
}


int initquerygeoid()
{
  opengeoid = -1; // No geoid files are open at start
}


int closequerygeoid()
{
  //printf("closing previously open geoid\n");
  if (opengeoid!=-1) fclose(fptrgeoid);
  opengeoid=-1;
}


double querytopo(double lat, double lon, char *demid)
{
  int bm2flag,g90flag;
  double topo,xs,ys,xn,yn;
  double bm2x[5],bm2y[5],g90x[5],g90y[5];
  double querygtopo30(double, double);
  double querybedmap2(double, double);
  double querygimp90(double, double);
  bool pointinpolygon(double,double,double *,double *,int);

  // Define the BEDMAP-2 polygon boundary
  bm2x[0] =-3333500.0;  bm2y[0] = 3333500.0;
  bm2x[1] = 3333500.0;  bm2y[1] = 3333500.0;
  bm2x[2] = 3333500.0;  bm2y[2] =-3333500.0;
  bm2x[3] =-3333500.0;  bm2y[3] =-3333500.0;
  bm2x[4] =-3333500.0;  bm2y[4] = 3333500.0;

  // Define the Greenland GIMP polygon boundary
  g90x[0] =-639955.0;  g90y[0] = -655595.0;
  g90x[1] = 855845.0;  g90y[1] = -655955.0;
  g90x[2] = 855845.0;  g90y[2] =-3355595.0;
  g90x[3] =-639955.0;  g90y[3] =-3355595.0;
  g90x[4] =-639955.0;  g90y[4] = -655595.0;

  // Determine if input coords are within Bedmap-2 or GIMP DEM limits
  geod2ps(lat,lon,-71.0,0.0,1.0,AE,FLAT,&xs,&ys);
  bm2flag = pointinpolygon(xs,ys,bm2x,bm2y,5);
  geod2ps(lat,lon,70.0,315.0,1.0,AE,FLAT,&xn,&yn);
  g90flag = pointinpolygon(xn,yn,g90x,g90y,5);

  // Query BEDMAP-2 surface DEM (Antarctica) if within its bounds
  if (bm2flag)
  {
    strcpy(demid,"BM2\0");
    topo = querybedmap2(xs,ys); // answer is relative to mean sea level
  }

  // Query GIMP90m DEM (Greenland) if within its bounds
  else if (g90flag)
  {
    strcpy(demid,"G90\0");
    topo = querygimp90(xn,yn);  // answer is relative to the WGS-84 ellipsoid
  }

  // Otherwise query GTOPO30
  else
  {
    strcpy(demid,"GT3\0");
    topo = querygtopo30(lat,lon);  // answer is relative to mean sea level
  }

  // Return the result
  return(topo);

}


double querygimp90(double x,double y)
{
  short sdtemp,q11,q21,q12,q22;
  int nx,ny,m1,m2,n1,n2;
  long int offsetq11,offsetq21,offsetq12,offsetq22;
  double x0,y0,mdbl,ndbl;
  double x1,x2,y1,y2,denom,f1,f2,f3,f4,p;
  void byteswap(char *,char *,int);

  // Define size of grid
  x0 =-639955.0;
  y0 =-655595.0;
  nx = 16620;
  ny = 30000;

  // Determine row and column of surrounding grid cells
  mdbl = (x-x0)/90.0-0.5;
  if (mdbl<0.0)
  {
    m1 = 0;
    m2 = 0;
  }
  else if (mdbl>(nx-1)) 
  {
    m1 = nx-1;
    m2 = nx-1;
  }
  else
  {
    m1 = int(mdbl);
    m2 = m1+1;
  }
  ndbl = (y0-y)/90.0-0.5;
  if (ndbl<0.0)
  {
    n1 = 0;
    n2 = 0;
  }
  else if (ndbl>(ny-1)) 
  {
    n1 = ny-1;
    n2 = ny-1;
  }
  else
  {
    n1 = int(ndbl);
    n2 = n1+1;
  }
  offsetq11 = 2*(n1*nx+m1);
  offsetq21 = 2*(n1*nx+m2);
  offsetq12 = 2*(n2*nx+m1);
  offsetq22 = 2*(n2*nx+m2);
  //printf("x: %lf  y: %lf\n",x,y);
  //printf("ncols: %d  nrows: %d\n",nx,ny);
  //printf("double col# %lf  double row# %lf\n",mdbl,ndbl);
  //printf("m1: %d  m2:%d\n",m1,m2);
  //printf("n1: %d  n2: %d\n",n1,n2);
  //printf("offsetq11=%ld  offsetq21=%ld\n",offsetq11,offsetq21);
  //printf("offsetq12=%ld  offsetq22=%ld\n",offsetq12,offsetq22);

  // Open GIMP90m DEM file if not already open
  if (opendem==-1)
  {
    //printf("no DEM open, opening GIMP90 DEM file\n");
    if ((fptrdem=fopen(GIMP90PATH,"r"))==NULL) return(-9999.9);
    opendem = 33;
  }
  else if (opendem!=33)
  {
    //printf("closing previously open DEM, opening GIMP90 DEM file\n");
    fclose(fptrdem);
    if ((fptrdem=fopen(GIMP90PATH,"r"))==NULL) return(-9999.9);
    opendem = 33;    
  }

  // Read the four surrounding pixels from the DEM file
  fseek(fptrdem,offsetq11,SEEK_SET);
  fread(&q11,2,1,fptrdem);
  if (q11==-9999.0) q11 = 0.0;
  fseek(fptrdem,offsetq21,SEEK_SET);
  fread(&q21,2,1,fptrdem);
  if (q21==-9999.0) q21 = 0.0;
  fseek(fptrdem,offsetq12,SEEK_SET);
  fread(&q12,2,1,fptrdem);
  if (q12==-9999.0) q12 = 0.0;
  fseek(fptrdem,offsetq22,SEEK_SET);
  fread(&q22,2,1,fptrdem);
  if (q22==-9999.0) q22 = 0.0;
  //printf("q11 = %d  q21 = %d\n",q11,q21);
  //printf("q12 = %d  q22 = %d\n",q12,q22);

  // Compute terrain at requested lat/lon by bilinear interpolation
  x1 = x0+(m1+0.5)*90.0;
  x2 = x0+(m2+0.5)*90.0;
  y1 = y0-(n1+0.5)*90.0;
  y2 = y0-(n2+0.5)*90.0;
  //printf("x1 = %lf  x2 = %lf\n",x1,x2);
  //printf("y1 = %lf  y2 = %lf\n",y1,y2);
  denom = (x2-x1)*(y2-y1);
  f1 = ((x2-x)*(y2-y))/denom;
  f2 = ((x-x1)*(y2-y))/denom;
  f3 = ((x2-x)*(y-y1))/denom;
  f4 = ((x-x1)*(y-y1))/denom;
  p = f1*q11 + f2*q21 + f3*q12 + f4*q22;
  return(p);

}


double querybedmap2(double x,double y)
{
  int nx,ny,m1,m2,n1,n2;
  long int offsetq11,offsetq21,offsetq12,offsetq22;
  float q11,q21,q12,q22;
  double x0,y0,mdbl,ndbl,x1,x2,y1,y2,denom,f1,f2,f3,f4,p;
  void byteswap(char *,char *,int);

  // Define size of grid
  x0 =-3333500.0;
  y0 = 3333500.0;
  nx = 6667;
  ny = 6667;

  // Determine row and column of surrounding grid cells
  mdbl = (x-x0)/1000.0-0.5;
  if (mdbl<0.0)
  {
    m1 = 0;
    m2 = 0;
  }
  else if (mdbl>(nx-1)) 
  {
    m1 = nx-1;
    m2 = nx-1;
  }
  else
  {
    m1 = int(mdbl);
    m2 = m1+1;
  }
  ndbl = (y0-y)/1000.0-0.5;
  if (ndbl<0.0)
  {
    n1 = 0;
    n2 = 0;
  }
  else if (ndbl>(ny-1)) 
  {
    n1 = ny-1;
    n2 = ny-1;
  }
  else
  {
    n1 = int(ndbl);
    n2 = n1+1;
  }
  offsetq11 = 4*(n1*nx+m1);
  offsetq21 = 4*(n1*nx+m2);
  offsetq12 = 4*(n2*nx+m1);
  offsetq22 = 4*(n2*nx+m2);
  //printf("ncols: %d  nrows: %d\n",nx,ny);
  //printf("double col# %lf  double row# %lf\n",mdbl,ndbl);
  //printf("m1: %d  m2:%d\n",m1,m2);
  //printf("n1: %d  n2: %d\n",n1,n2);
  //printf("offsetq11=%ld  offsetq21=%ld\n",offsetq11,offsetq21);
  //printf("offsetq12=%ld  offsetq22=%ld\n",offsetq12,offsetq22);

  // Open BEDMAP-2 DEM file if not already open
  if (opendem==-1)
  {
    //printf("no DEM open, opening Bedmap-2 DEM file\n");
    if ((fptrdem=fopen(BEDMAP2PATH,"r"))==NULL) return(-9999.9);
    opendem = 34;
  }
  else if (opendem!=34)
  {
    //printf("closing previously open DEM, opening Bedmap-2 DEM file\n");
    fclose(fptrdem);
    if ((fptrdem=fopen(BEDMAP2PATH,"r"))==NULL) return(-9999.9);
    opendem = 34;    
  }

  // Read the four surrounding pixels from the DEM file
  fseek(fptrdem,offsetq11,SEEK_SET);
  fread(&q11,4,1,fptrdem);
  if (q11==-9999.0) q11 = 0.0;
  fseek(fptrdem,offsetq21,SEEK_SET);
  fread(&q21,4,1,fptrdem);
  if (q21==-9999.0) q21 = 0.0;
  fseek(fptrdem,offsetq12,SEEK_SET);
  fread(&q12,4,1,fptrdem);
  if (q12==-9999.0) q12 = 0.0;
  fseek(fptrdem,offsetq22,SEEK_SET);
  fread(&q22,4,1,fptrdem);
  if (q22==-9999.0) q22 = 0.0;
  //printf("q11 = %f  q21 = %f\n",q11,q21);
  //printf("q12 = %f  q22 = %f\n",q12,q22);

  // Compute terrain at requested lat/lon by bilinear interpolation
  x1 = x0+(m1+0.5)*1000.0;
  x2 = x0+(m2+0.5)*1000.0;
  y1 = y0-(n1+0.5)*1000.0;
  y2 = y0-(n2+0.5)*1000.0;
  //printf("x1 = %lf  x2 = %lf\n",x1,x2);
  //printf("y1 = %lf  y2 = %lf\n",y1,y2);
  denom = (x2-x1)*(y2-y1);
  f1 = ((x2-x)*(y2-y))/denom;
  f2 = ((x-x1)*(y2-y))/denom;
  f3 = ((x2-x)*(y-y1))/denom;
  f4 = ((x-x1)*(y-y1))/denom;
  p = f1*q11 + f2*q21 + f3*q12 + f4*q22;
  return(p);

}


double querygtopo30(double lat,double lon)
{
  bool flag;
  char tilename[35][15],filename[120];
  short sdtemp,q11,q21,q12,q22;
  int ntiles,i,nlon,nlat,n1,n2,m1,m2;
  long int offsetq11,offsetq21,offsetq12,offsetq22;
  double mdbl,ndbl,lon1,lon2,lat1,lat2;
  double denom,f1,f2,f3,f4,p;
  double tilelat[35][5],tilelon[40][5];
  bool pointinpolygon(double,double,double *,double *,int);
  void byteswap(char *,char *,int);

  // Define the GTOPO30 tile boundaries
  ntiles = 33;
  strcpy(tilename[ 0],"W180N90.DEM\0");
  tilelat[ 0][0] = 90.0;  tilelon[ 0][0] =-180.0;
  tilelat[ 0][1] = 90.0;  tilelon[ 0][1] =-140.0;
  tilelat[ 0][2] = 40.0;  tilelon[ 0][2] =-140.0;
  tilelat[ 0][3] = 40.0;  tilelon[ 0][3] =-180.0;
  tilelat[ 0][4] = 90.0;  tilelon[ 0][4] =-180.0;
  strcpy(tilename[ 1],"W140N90.DEM\0");
  tilelat[ 1][0] = 90.0;  tilelon[ 1][0] =-140.0;
  tilelat[ 1][1] = 90.0;  tilelon[ 1][1] =-100.0;
  tilelat[ 1][2] = 40.0;  tilelon[ 1][2] =-100.0;
  tilelat[ 1][3] = 40.0;  tilelon[ 1][3] =-140.0;
  tilelat[ 1][4] = 90.0;  tilelon[ 1][4] =-140.0;
  strcpy(tilename[ 2],"W100N90.DEM\0");
  tilelat[ 2][0] = 90.0;  tilelon[ 2][0] =-100.0;
  tilelat[ 2][1] = 90.0;  tilelon[ 2][1] =-060.0;
  tilelat[ 2][2] = 40.0;  tilelon[ 2][2] =-060.0;
  tilelat[ 2][3] = 40.0;  tilelon[ 2][3] =-100.0;
  tilelat[ 2][4] = 90.0;  tilelon[ 2][4] =-100.0;
  strcpy(tilename[ 3],"W060N90.DEM\0");
  tilelat[ 3][0] = 90.0;  tilelon[ 3][0] =-060.0;
  tilelat[ 3][1] = 90.0;  tilelon[ 3][1] =-020.0;
  tilelat[ 3][2] = 40.0;  tilelon[ 3][2] =-020.0;
  tilelat[ 3][3] = 40.0;  tilelon[ 3][3] =-060.0;
  tilelat[ 3][4] = 90.0;  tilelon[ 3][4] =-060.0;
  strcpy(tilename[ 4],"W020N90.DEM\0");
  tilelat[ 4][0] = 90.0;  tilelon[ 4][0] =-020.0;
  tilelat[ 4][1] = 90.0;  tilelon[ 4][1] = 020.0;
  tilelat[ 4][2] = 40.0;  tilelon[ 4][2] = 020.0;
  tilelat[ 4][3] = 40.0;  tilelon[ 4][3] =-020.0;
  tilelat[ 4][4] = 90.0;  tilelon[ 4][4] =-020.0;
  strcpy(tilename[ 5],"E020N90.DEM\0");
  tilelat[ 5][0] = 90.0;  tilelon[ 5][0] = 020.0;
  tilelat[ 5][1] = 90.0;  tilelon[ 5][1] = 060.0;
  tilelat[ 5][2] = 40.0;  tilelon[ 5][2] = 060.0;
  tilelat[ 5][3] = 40.0;  tilelon[ 5][3] = 020.0;
  tilelat[ 5][4] = 90.0;  tilelon[ 5][4] = 020.0;
  strcpy(tilename[ 6],"E060N90.DEM\0");
  tilelat[ 6][0] = 90.0;  tilelon[ 6][0] = 060.0;
  tilelat[ 6][1] = 90.0;  tilelon[ 6][1] = 100.0;
  tilelat[ 6][2] = 40.0;  tilelon[ 6][2] = 100.0;
  tilelat[ 6][3] = 40.0;  tilelon[ 6][3] = 060.0;
  tilelat[ 6][4] = 90.0;  tilelon[ 6][4] = 060.0;
  strcpy(tilename[ 7],"E100N90.DEM\0");
  tilelat[ 7][0] = 90.0;  tilelon[ 7][0] = 100.0;
  tilelat[ 7][1] = 90.0;  tilelon[ 7][1] = 140.0;
  tilelat[ 7][2] = 40.0;  tilelon[ 7][2] = 140.0;
  tilelat[ 7][3] = 40.0;  tilelon[ 7][3] = 100.0;
  tilelat[ 7][4] = 90.0;  tilelon[ 7][4] = 100.0;
  strcpy(tilename[ 8],"E140N90.DEM\0");
  tilelat[ 8][0] = 90.0;  tilelon[ 8][0] = 140.0;
  tilelat[ 8][1] = 90.0;  tilelon[ 8][1] = 180.0;
  tilelat[ 8][2] = 40.0;  tilelon[ 8][2] = 180.0;
  tilelat[ 8][3] = 40.0;  tilelon[ 8][3] = 140.0;
  tilelat[ 8][4] = 90.0;  tilelon[ 8][4] = 140.0;
  strcpy(tilename[ 9],"W180N40.DEM\0");
  tilelat[ 9][0] = 40.0;  tilelon[ 9][0] =-180.0;
  tilelat[ 9][1] = 40.0;  tilelon[ 9][1] =-140.0;
  tilelat[ 9][2] =-10.0;  tilelon[ 9][2] =-140.0;
  tilelat[ 9][3] =-10.0;  tilelon[ 9][3] =-180.0;
  tilelat[ 9][4] = 40.0;  tilelon[ 9][4] =-180.0;
  strcpy(tilename[10],"W140N40.DEM\0");
  tilelat[10][0] = 40.0;  tilelon[10][0] =-140.0;
  tilelat[10][1] = 40.0;  tilelon[10][1] =-100.0;
  tilelat[10][2] =-10.0;  tilelon[10][2] =-100.0;
  tilelat[10][3] =-10.0;  tilelon[10][3] =-140.0;
  tilelat[10][4] = 40.0;  tilelon[10][4] =-140.0;
  strcpy(tilename[11],"W100N40.DEM\0");
  tilelat[11][0] = 40.0;  tilelon[11][0] =-100.0;
  tilelat[11][1] = 40.0;  tilelon[11][1] =-060.0;
  tilelat[11][2] =-10.0;  tilelon[11][2] =-060.0;
  tilelat[11][3] =-10.0;  tilelon[11][3] =-100.0;
  tilelat[11][4] = 40.0;  tilelon[11][4] =-100.0;
  strcpy(tilename[12],"W060N40.DEM\0");
  tilelat[12][0] = 40.0;  tilelon[12][0] =-060.0;
  tilelat[12][1] = 40.0;  tilelon[12][1] =-020.0;
  tilelat[12][2] =-10.0;  tilelon[12][2] =-020.0;
  tilelat[12][3] =-10.0;  tilelon[12][3] =-060.0;
  tilelat[12][4] = 40.0;  tilelon[12][4] =-060.0;
  strcpy(tilename[13],"W020N40.DEM\0");
  tilelat[13][0] = 40.0;  tilelon[13][0] =-020.0;
  tilelat[13][1] = 40.0;  tilelon[13][1] = 020.0;
  tilelat[13][2] =-10.0;  tilelon[13][2] = 020.0;
  tilelat[13][3] =-10.0;  tilelon[13][3] =-020.0;
  tilelat[13][4] = 40.0;  tilelon[13][4] =-020.0;
  strcpy(tilename[14],"E020N40.DEM\0");
  tilelat[14][0] = 40.0;  tilelon[14][0] = 020.0;
  tilelat[14][1] = 40.0;  tilelon[14][1] = 060.0;
  tilelat[14][2] =-10.0;  tilelon[14][2] = 060.0;
  tilelat[14][3] =-10.0;  tilelon[14][3] = 020.0;
  tilelat[14][4] = 40.0;  tilelon[14][4] = 020.0;
  strcpy(tilename[15],"E060N40.DEM\0");
  tilelat[15][0] = 40.0;  tilelon[15][0] = 060.0;
  tilelat[15][1] = 40.0;  tilelon[15][1] = 100.0;
  tilelat[15][2] =-10.0;  tilelon[15][2] = 100.0;
  tilelat[15][3] =-10.0;  tilelon[15][3] = 060.0;
  tilelat[15][4] = 40.0;  tilelon[15][4] = 060.0;
  strcpy(tilename[16],"E100N40.DEM\0");
  tilelat[16][0] = 40.0;  tilelon[16][0] = 100.0;
  tilelat[16][1] = 40.0;  tilelon[16][1] = 140.0;
  tilelat[16][2] =-10.0;  tilelon[16][2] = 140.0;
  tilelat[16][3] =-10.0;  tilelon[16][3] = 100.0;
  tilelat[16][4] = 40.0;  tilelon[16][4] = 100.0;
  strcpy(tilename[17],"E140N40.DEM\0");
  tilelat[17][0] = 40.0;  tilelon[17][0] = 140.0;
  tilelat[17][1] = 40.0;  tilelon[17][1] = 180.0;
  tilelat[17][2] =-10.0;  tilelon[17][2] = 180.0;
  tilelat[17][3] =-10.0;  tilelon[17][3] = 140.0;
  tilelat[17][4] = 40.0;  tilelon[17][4] = 140.0;
  strcpy(tilename[18],"W180S10.DEM\0");
  tilelat[18][0] =-10.0;  tilelon[18][0] =-180.0;
  tilelat[18][1] =-10.0;  tilelon[18][1] =-140.0;
  tilelat[18][2] =-60.0;  tilelon[18][2] =-140.0;
  tilelat[18][3] =-60.0;  tilelon[18][3] =-180.0;
  tilelat[18][4] =-10.0;  tilelon[18][4] =-180.0;
  strcpy(tilename[19],"W140S10.DEM\0");
  tilelat[19][0] =-10.0;  tilelon[19][0] =-140.0;
  tilelat[19][1] =-10.0;  tilelon[19][1] =-100.0;
  tilelat[19][2] =-60.0;  tilelon[19][2] =-100.0;
  tilelat[19][3] =-60.0;  tilelon[19][3] =-140.0;
  tilelat[19][4] =-10.0;  tilelon[19][4] =-140.0;
  strcpy(tilename[20],"W100S10.DEM\0");
  tilelat[20][0] =-10.0;  tilelon[20][0] =-100.0;
  tilelat[20][1] =-10.0;  tilelon[20][1] =-060.0;
  tilelat[20][2] =-60.0;  tilelon[20][2] =-060.0;
  tilelat[20][3] =-60.0;  tilelon[20][3] =-100.0;
  tilelat[20][4] =-10.0;  tilelon[20][4] =-100.0;
  strcpy(tilename[21],"W060S10.DEM\0");
  tilelat[21][0] =-10.0;  tilelon[21][0] =-060.0;
  tilelat[21][1] =-10.0;  tilelon[21][1] =-020.0;
  tilelat[21][2] =-60.0;  tilelon[21][2] =-020.0;
  tilelat[21][3] =-60.0;  tilelon[21][3] =-060.0;
  tilelat[21][4] =-10.0;  tilelon[21][4] =-060.0;
  strcpy(tilename[22],"W020S10.DEM\0");
  tilelat[22][0] =-10.0;  tilelon[22][0] =-020.0;
  tilelat[22][1] =-10.0;  tilelon[22][1] = 020.0;
  tilelat[22][2] =-60.0;  tilelon[22][2] = 020.0;
  tilelat[22][3] =-60.0;  tilelon[22][3] =-020.0;
  tilelat[22][4] =-10.0;  tilelon[22][4] =-020.0;
  strcpy(tilename[23],"W020S10.DEM\0");
  tilelat[23][0] =-10.0;  tilelon[23][0] = 020.0;
  tilelat[23][1] =-10.0;  tilelon[23][1] = 060.0;
  tilelat[23][2] =-60.0;  tilelon[23][2] = 060.0;
  tilelat[23][3] =-60.0;  tilelon[23][3] = 020.0;
  tilelat[23][4] =-10.0;  tilelon[23][4] = 020.0;
  strcpy(tilename[24],"W060S10.DEM\0");
  tilelat[24][0] =-10.0;  tilelon[24][0] = 060.0;
  tilelat[24][1] =-10.0;  tilelon[24][1] = 100.0;
  tilelat[24][2] =-60.0;  tilelon[24][2] = 100.0;
  tilelat[24][3] =-60.0;  tilelon[24][3] = 060.0;
  tilelat[24][4] =-10.0;  tilelon[24][4] = 060.0;
  strcpy(tilename[25],"W100S10.DEM\0");
  tilelat[25][0] =-10.0;  tilelon[25][0] = 100.0;
  tilelat[25][1] =-10.0;  tilelon[25][1] = 140.0;
  tilelat[25][2] =-60.0;  tilelon[25][2] = 140.0;
  tilelat[25][3] =-60.0;  tilelon[25][3] = 100.0;
  tilelat[25][4] =-10.0;  tilelon[25][4] = 100.0;
  strcpy(tilename[26],"E140S10.DEM\0");
  tilelat[26][0] =-10.0;  tilelon[26][0] = 140.0;
  tilelat[26][1] =-10.0;  tilelon[26][1] = 180.0;
  tilelat[26][2] =-60.0;  tilelon[26][2] = 180.0;
  tilelat[26][3] =-60.0;  tilelon[26][3] = 140.0;
  tilelat[26][4] =-10.0;  tilelon[26][4] = 140.0;
  strcpy(tilename[27],"W180S60.DEM\0");
  tilelat[27][0] =-60.0;  tilelon[27][0] =-180.0;
  tilelat[27][1] =-60.0;  tilelon[27][1] =-120.0;
  tilelat[27][2] =-90.0;  tilelon[27][2] =-120.0;
  tilelat[27][3] =-90.0;  tilelon[27][3] =-180.0;
  tilelat[27][4] =-60.0;  tilelon[27][4] =-180.0;
  strcpy(tilename[28],"W120S60.DEM\0");
  tilelat[28][0] =-60.0;  tilelon[28][0] =-120.0;
  tilelat[28][1] =-60.0;  tilelon[28][1] =-060.0;
  tilelat[28][2] =-90.0;  tilelon[28][2] =-060.0;
  tilelat[28][3] =-90.0;  tilelon[28][3] =-120.0;
  tilelat[28][4] =-60.0;  tilelon[28][4] =-120.0;
  strcpy(tilename[29],"W060S60.DEM\0");
  tilelat[29][0] =-60.0;  tilelon[29][0] =-060.0;
  tilelat[29][1] =-60.0;  tilelon[29][1] = 000.0;
  tilelat[29][2] =-90.0;  tilelon[29][2] = 000.0;
  tilelat[29][3] =-90.0;  tilelon[29][3] =-060.0;
  tilelat[29][4] =-60.0;  tilelon[29][4] =-060.0;
  strcpy(tilename[30],"W000S60.DEM\0");
  tilelat[30][0] =-60.0;  tilelon[30][0] = 000.0;
  tilelat[30][1] =-60.0;  tilelon[30][1] = 060.0;
  tilelat[30][2] =-90.0;  tilelon[30][2] = 060.0;
  tilelat[30][3] =-90.0;  tilelon[30][3] = 000.0;
  tilelat[30][4] =-60.0;  tilelon[30][4] = 000.0;
  strcpy(tilename[31],"E060S60.DEM\0");
  tilelat[31][0] =-60.0;  tilelon[31][0] = 060.0;
  tilelat[31][1] =-60.0;  tilelon[31][1] = 120.0;
  tilelat[31][2] =-90.0;  tilelon[31][2] = 120.0;
  tilelat[31][3] =-90.0;  tilelon[31][3] = 060.0;
  tilelat[31][4] =-60.0;  tilelon[31][4] = 060.0;
  strcpy(tilename[32],"E120S60.DEM\0");
  tilelat[32][0] =-60.0;  tilelon[32][0] = 120.0;
  tilelat[32][1] =-60.0;  tilelon[32][1] = 180.0;
  tilelat[32][2] =-90.0;  tilelon[32][2] = 180.0;
  tilelat[32][3] =-90.0;  tilelon[32][3] = 120.0;
  tilelat[32][4] =-60.0;  tilelon[32][4] = 120.0;

  // Loop over all tiles
  for (i=0;i<33;i++)
  {

    // If point is within current tile, query that tile
    flag = pointinpolygon(lon,lat,tilelon[i],tilelat[i],5);
    if (flag) 
    {
      //printf("%s\n",tilename[i]);
      nlon = int(120.0*(tilelon[i][1]-tilelon[i][0]));
      nlat = int(120.0*(tilelat[i][1]-tilelat[i][2]));
      mdbl = 120.0*(lon-tilelon[i][0])-0.5;
      if (mdbl<0.0)
      {
        m1 = 0;
        m2 = 0;
      }
      else if (mdbl>(nlon-1)) 
      {
        m1 = nlon-1;
        m2 = nlon-1;
      }
      else
      {
        m1 = int(mdbl);
        m2 = m1+1;
      }
      ndbl = 120.0*(tilelat[i][0]-lat)-0.5;
      if (ndbl<0.0)
      {
        n1 = 0;
        n2 = 0;
      }
      else if (ndbl>(nlat-1)) 
      {
        n1 = nlat-1;
        n2 = nlat-1;
      }
      else
      {
        n1 = int(ndbl);
        n2 = n1+1;
      }
      offsetq11 = 2*(n1*nlon+m1);
      offsetq21 = 2*(n1*nlon+m2);
      offsetq12 = 2*(n2*nlon+m1);
      offsetq22 = 2*(n2*nlon+m2);
      //printf("ncols: %d  nrows: %d\n",nlon,nlat);
      //printf("double col# %lf  double row# %lf\n",mdbl,ndbl);
      //printf("m1: %d  m2:%d\n",m1,m2);
      //printf("n1: %d  n2: %d\n",n1,n2);
      //printf("offsetq11=%ld  offsetq21=%ld\n",offsetq11,offsetq21);
      //printf("offsetq12=%ld  offsetq22=%ld\n",offsetq12,offsetq22);

      // Open this GTOPO30 DEM tile if not already open
      strcpy(filename,GTOPO30PATH);
      strcat(filename,tilename[i]);
      if (opendem==-1)
      {
        //printf("no DEM open, opening GTOPO30 DEM tile %d\n",i);
        if ((fptrdem=fopen(filename,"r"))==NULL) return(-9999.9);
        opendem = i;
      }
      else if (opendem!=i)
      {
        //printf("closing previously open DEM, opening GTOPO30 DEM tile %d\n",i);
        fclose(fptrdem);
        if ((fptrdem=fopen(filename,"r"))==NULL) return(-9999.9);
        opendem = i;    
      }

      // Read the four surrounding pixels from the DEM file
      fseek(fptrdem,offsetq11,SEEK_SET);
      fread(&sdtemp,2,1,fptrdem);
      byteswap((char *)&sdtemp,(char *)&q11,2);
      if (q11==-9999) q11=0;
      fseek(fptrdem,offsetq21,SEEK_SET);
      fread(&sdtemp,2,1,fptrdem);
      byteswap((char *)&sdtemp,(char *)&q21,2);
      if (q21==-9999) q21=0;
      fseek(fptrdem,offsetq12,SEEK_SET);
      fread(&sdtemp,2,1,fptrdem);
      byteswap((char *)&sdtemp,(char *)&q12,2);
      if (q12==-9999) q12=0;
      fseek(fptrdem,offsetq22,SEEK_SET);
      fread(&sdtemp,2,1,fptrdem);
      byteswap((char *)&sdtemp,(char *)&q22,2);
      if (q22==-9999) q22=0;
      //printf("q11 = %d  q21 = %d\n",q11,q21);
      //printf("q12 = %d  q22 = %d\n",q12,q22);
      
      // Compute terrain at requested lat/lon by bilinear interpolation
      lon1 = tilelon[i][0]+(m1+0.5)/120.0;
      lon2 = tilelon[i][0]+(m2+0.5)/120.0;
      lat1 = tilelat[i][0]-(n1+0.5)/120.0;
      lat2 = tilelat[i][0]-(n2+0.5)/120.0;
      //printf("lon1 = %lf  lon2 = %lf\n",lon1,lon2);
      //printf("lat1 = %lf  lat2 = %lf\n",lat1,lat2);
      denom = (lon2-lon1)*(lat2-lat1);
      f1 = ((lon2-lon)*(lat2-lat))/denom;
      f2 = ((lon-lon1)*(lat2-lat))/denom;
      f3 = ((lon2-lon)*(lat-lat1))/denom;
      f4 = ((lon-lon1)*(lat-lat1))/denom;
      p = f1*q11 + f2*q21 + f3*q12 + f4*q22;
      return(p);

    }

  }

}


double querygeoid(double lat, double lon, char *geoidid)
{
  double geoid;
  double queryegm96(double, double);
  
  // EGM96 is our only available geoid currently
  if (1)
  {
    strcpy(geoidid,"E96\0");
    geoid = queryegm96(lat,lon);
  }

}


double queryegm96(double y,double x)
{
  short sdtemp,q11,q21,q12,q22;
  int nx,ny,mdbl,ndbl,m1,m2,n1,n2;
  double x0,y0,lon1,lon2,lat1,lat2,denom,f1,f2,f3,f4,p;
  long int offsetq11,offsetq21,offsetq12,offsetq22;
  void byteswap(char *,char *,int);

  // Offset longitude to between 0 and 360
  while (x<0.0) x+=360.0;

  // Define size of grid
  x0 = 0.0;
  y0 = 90.0;
  nx = 1440;
  ny = 721;

  // Determine row and column of surrounding grid cells
  mdbl = 4.0*(x-x0);
  if (mdbl<0.0)
  {
    m1 = 0;
    m2 = 0;
  }
  else if (mdbl>(nx-1)) 
  {
    m1 = nx-1;
    m2 = nx-1;
  }
  else
  {
    m1 = int(mdbl);
    m2 = m1+1;
  }
  ndbl = 4.0*(y0-y);
  if (ndbl<0.0)
  {
    n1 = 0;
    n2 = 0;
  }
  else if (ndbl>(ny-1)) 
  {
    n1 = ny-1;
    n2 = ny-1;
  }
  else
  {
    n1 = int(ndbl);
    n2 = n1+1;
  }
  offsetq11 = 2*(n1*nx+m1);
  offsetq21 = 2*(n1*nx+m2);
  offsetq12 = 2*(n2*nx+m1);
  offsetq22 = 2*(n2*nx+m2);
  //printf("ncols: %d  nrows: %d\n",nx,ny);
  //printf("double col# %lf  double row# %lf\n",mdbl,ndbl);
  //printf("m1: %d  m2:%d\n",m1,m2);
  //printf("n1: %d  n2: %d\n",n1,n2);
  //printf("offsetq11=%ld  offsetq21=%ld\n",offsetq11,offsetq21);
  //printf("offsetq12=%ld  offsetq22=%ld\n",offsetq12,offsetq22);

  // Open EGM96 geoid file if not already open
  if (opengeoid==-1)
  {
    if ((fptrgeoid=fopen(EGM96PATH,"r"))==NULL) return(-9999.9);
    opengeoid = 0;
  }

  // Read the four surrounding pixels from the geoid file
  fseek(fptrgeoid,offsetq11,SEEK_SET);
  fread(&sdtemp,2,1,fptrgeoid);
  byteswap((char *)&sdtemp,(char *)&q11,2);
  if (q11==-9999.0) q11 = 0.0;
  fseek(fptrgeoid,offsetq21,SEEK_SET);
  fread(&sdtemp,2,1,fptrgeoid);
  byteswap((char *)&sdtemp,(char *)&q21,2);
  if (q21==-9999.0) q21 = 0.0;
  fseek(fptrgeoid,offsetq12,SEEK_SET);
  fread(&sdtemp,2,1,fptrgeoid);
  byteswap((char *)&sdtemp,(char *)&q12,2);
  if (q12==-9999.0) q12 = 0.0;
  fseek(fptrgeoid,offsetq22,SEEK_SET);
  fread(&sdtemp,2,1,fptrgeoid);
  byteswap((char *)&sdtemp,(char *)&q22,2);
  if (q22==-9999.0) q22 = 0.0;
  //printf("q11 = %d  q21 = %d\n",q11,q21);
  //printf("q12 = %d  q22 = %d\n",q12,q22);

  // Compute geoid at requested lat/lon by bilinear interpolation
  lon1 = x0+m1/4.0;
  lon2 = x0+m2/4.0;
  lat1 = y0-n1/4.0;
  lat2 = y0-n2/4.0;
  //printf("lon1 = %lf  lon2 = %lf\n",lon1,lon2);
  //printf("lat1 = %lf  lat2 = %lf\n",lat1,lat2);
  denom = (lon2-lon1)*(lat2-lat1);
  f1 = ((lon2-x)*(lat2-y))/denom;
  f2 = ((x-lon1)*(lat2-y))/denom;
  f3 = ((lon2-x)*(y-lat1))/denom;
  f4 = ((x-lon1)*(y-lat1))/denom;
  p = f1*q11 + f2*q21 + f3*q12 + f4*q22;
  return(p/100.0); // heights in database are in cm

}


bool pointinpolygon(double x, double y,double xpoly[],double ypoly[],int npoly)
{
  int i,j=npoly-2;
  bool oddnodes=false;

  for (i=0; i<(npoly-1); i++) 
  {
    if (ypoly[i]<y && ypoly[j]>=y || ypoly[j]<y && ypoly[i]>=y) 
    {
      if (xpoly[i]+(y-ypoly[i])/(ypoly[j]-ypoly[i])*(xpoly[j]-xpoly[i])<x) 
      {
        oddnodes=!oddnodes; 
      }
    }
    j=i; 
  }

  return(oddnodes);

}


void byteswap(char *in,char *out,int len)
{
  int i,len2;

  len2 = len-1;
  for (i=0;i<len;i++)
  {
    out[i]=in[len2-i];
    //out[i] = in[i];
  }

}

