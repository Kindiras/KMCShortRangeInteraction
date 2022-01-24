#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#define DEBUGCHECK
//#define DEBUG  1
//#define DEBUGPIC 1
//#define PICTURE 1
#define NRUN 1 // number of runs
#define N 1000 // system size
double diffrate=1.e4;//6.0e5;// this is monomer diffusion rate
double barrier=0.1; //0.01;
#define MAXISLSIZE 100000
#define MAXISLCNT 1000000
#define MAXNUMMONOMERS  2000000
#define CMIN 1.e-5 // this is minimum coverage
#define CMAX 0.3 // this is maximum coverage (volume fraction)
#define NDATA 500// number of data points
#define NDATAp1 (NDATA+1)
//#define NCAP 5//data points for capture routine(+1)
#define NDIR 4 // number of directions for diffusion
#define NDIRM 3 // NDIR - 1
#define NDIRH 2 //NDIR/2
#define nn6 (N*N*NDIR) //number of sites * number of directions
#define nn6small (N*N*NDIR/8)
#define SEED 305 //time(NULL) //305  Changed this to a random seed to give me peace of mind- Brad

double Rate[8];
double suni(),xdouble,uni();
double cova[NDATAp1];//array of coverages at which densa and mondensa
// are filled 
double cova[NDATAp1];//array of coverages at which densa and mondensa
double dxcov,xnsq,xattach,xattachcount;
double sigma[MAXISLCNT][NCAP],flsigma[MAXISLCNT][NCAP],sigma2[MAXISLCNT][NCAP];
double avsigma[NCAP],avsigma2[NCAP],flavsigma[NCAP],avmonden,mondensbeg;
double isdi[MAXISLCNT];
double tsigma[MAXISLCNT];
double mondensa[NDATAp1];//monomer density array (as function of coverage)
double mondensqa[NDATAp1];//monomer density array (as function of coverage)
double nwa[NDATAp1];//walker density array (as function of coverage)
double densa[NDATAp1];//island density array (as function of coverage)
double dens3a[NDATAp1];//island density array (as function of coverage)
double dens7a[NDATAp1];//island density array (as function of coverage)
double sdnrun,dnrun;
double sav;
double sava[5];
double savloga[NDATAp1];//array to save log data in
double rnddf,rnddf1,rnddf2,diffq;
double wa[NDATAp1],skewa[NDATAp1];
double cth,sth;
double covlog,sumtot;

unsigned long poww2(int j);

int nw[8],ip[N],im[N];
int islsiza[MAXISLSIZE];
int list[4][MAXNUMMONOMERS];
int idisplay=0;
//int h[N][N],imark[N][N];
unsigned char *h;
unsigned char *imark;
//int ipointa[nn6], indexa[nn6];
unsigned char *indexa;
int *ipointa;
//int islsize[N][N];
int idata;
int ideposit=0;
int idiffuse=0;
int isize, islcnta[5][MAXISLCNT],isizemax,x;
int numisizea[MAXISLCNT],ismax,ismaxi;
int ct[MAXISLCNT][NCAP];
int skip;
int countflag;
int nbxa[N][6], nbya[N][6];
int nwt,ivar;
int np,npmax,npincr;
int iii,jjj,isiteflat,iflat,jflat;
int hmax;
int numWalker=0;
int pick=0;
int nsq,sum02;

FILE *wdata,*picure,*skdata,*tandata,*ggdata,*pixdata,*slopedat;
FILE *dist1,*dist2,*dist3,*dist4,*dist5;
FILE *logdist1,*logdist2,*logdist3,*logdist4,*logdist5;

void debugcheck();
void upnbhd(int i, int j);
void upnbhdspecial(int i, int j);
void add(int newindex, int isite);
void findislandsize();
void depositmonomer();
void diffuse();
void takelogdata (int ilogdata);
void takedistdata ();
void check(int i, int j, int numcl);
void search(int ii,int jj);
int icount(int i, int j, int dir);
int update(int i, int j,  int dir);
void updatefix(int i, int j, int dir);
void delete(int index, int isite);

FILE *densdata,*picdata,*capdata,*capdivdata,*capdist1,*capdist2,*capdist3,*capdist4,*capdist5;

int main()
{
  int   long jseed,icdata;
  double deprate, totalRate,dcov,covv,cov,logdcov;
  double chance=0;
  double dmsigma,msigma,density;
  double numLayers=0;
  double xtemp,xx,yy,zz,x1,sig1,sig2,sig3,x;
  double logx,logy,logz,covnext,covfactor;
  int numtotal,iloop,irun,totaldropped,ilogdata;
  int index=0;
  int counter=0;
  int i=0;
  int j=0;
  int k=0;
  int jjlog;
  double savvlog,totaldiff,sum,sumxx;
  double temperature=700.0; //temperature of the system
  double kb=8.617e-5; //boltzman in eV/K
  double v_mono=1e13/4.0; //frequency for monomer hopping / 6 since there are six directions
  double Es=0.6;//monomer diffusion energy barrier
  //  double Eb=atof(argv[1]);//0.65; //detachment barrier
  double Ee=0.6; //edge diffusion barrier 0.6 for fast diffusion and 10 for slow
  double kT=kb*temperature;
  double DoF=2.4e9;// D/F
  
  double flux=1.0;//1.0;//diffrate/DoF;  //2.5e-4;//1.0;//this is deposition rate
  double r1det=diffrate*13.0; //0.01;//this is 1-bond detach
  double r1edge=diffrate*13.0; //0.1;//this is 1-bond edge
  double r2det=diffrate*13.0*0.025; //0.0001;//this is 2-bond detach for 4.90 nm particles
  double r2edge=diffrate*13.0*0.025; //0.001;//this is 2-bond edge for 4.9 nm particles
  double r3det=diffrate*13.0*0.000625;//this is 3-bond detach for 4.0 nm particles
  dnrun = (double) NRUN;
  sdnrun = sqrt(dnrun-1.0);
  if (NRUN==1) sdnrun=1.0;
  xnsq=N*N;
  dcov=CMAX/5;
  xnsq = (double) xnsq;

//int h[N][N],imark[N][N];
//unsigned char *h;
//unsigned char *imark;
//int ipointa[nn6], indexa[nn6];
//short *indexa;
//int *ipointa;

  h = (unsigned char *)malloc(sizeof(unsigned char)*(N*N));
  imark = (unsigned char *)malloc(sizeof(unsigned char)*(N*N));
  indexa = (unsigned char *)malloc(sizeof(unsigned char)*(N*N*NDIR));
  ipointa = (int *)malloc(sizeof(int)*(N*N*NDIR));

  
  jseed=SEED;
  xdouble=suni(jseed);
  densdata=fopen("densbarrcodenormalN1kr1e8barr1p0-1r-.3ML","w");
  dist1=fopen("isdbarrcodenormalsqN1kr1e8barr1p0-1r-.06ML","w");
  dist2=fopen("isdbarrcodenormalsqN1kr1e8barr1p0-1r-.12ML","w");
  dist3=fopen("isdbarrcodenormalsqN1kr1e8barr1p0-1r-.18ML","w");
  dist4=fopen("isdbarrcodenormalsqN1kr1e8barr1p0-1r-.24ML","w");
  dist5=fopen("isdbarrcodenormalsqN1kr1e8barr1p0-1r-.3ML","w");
  picdata=fopen("picbarrcodenormalsqN1kr1e8barr1p0-1r-.3ML","w");

  nsq=N*N;

#ifdef DEBUG
  printf("about to initialize direction arrays\n");
#endif
  // now initialize direction pointers nbxa, nbya
  for(i=0;i<N;i++){
    ip[i]=i+1;
    im[i]=i-1;  }                                                     
  ip[N-1]=0;
  im[0]=N-1;
  for(i=0;i<N;i++){
    //c dir 0 is -> (E)
    nbxa[i][0]=ip[i];
    nbya[i][0]=i;
    //c
    //c dir 1 is S
    nbxa[i][1]=i;
    nbya[i][1]=im[i];
    //c
    //c dir 2 is W
    nbxa[i][2]=im[i];
    nbya[i][2]=i;
    //c
    //c dir 3 is W
    nbxa[i][3]=i;
    nbya[i][3]=ip[i];
    //c
    //c dir 4 is NW
    nbxa[i][4]=im[i];
    nbya[i][4]=ip[i];
    //c
    //c dir 5 is NE
    nbxa[i][5]=i;
    nbya[i][5]=ip[i];
  }
  Rate[0]=diffrate;
  Rate[1]=r1edge;//Rate[0]*r1edge;
  Rate[2]=r1det;//Rate[0]*r1det;
  Rate[3]=r2edge;//Rate[0]*r2edge;
  Rate[4]=r2det;//Rate[0]*r2det;
  Rate[5]=r3det;
  Rate[6]=0.0;
  Rate[7]=0.0;
  #ifdef DEBUG
   for(i=0;i<8;i++){
      printf("rate(%d)=%f\n",i,Rate[i]);
     }
 #endif
  // now initialize data arrays nwa[NDATAp1]  etc.
  for(i=0;i<NDATAp1;i++){
    nwa[i]=0; 
    mondensa[i]=0;
    mondensqa[i]=0;
    densa[i]=0; 
    dens3a[i]=0; 
    dens7a[i]=0; 
  }
  for(i=0;i<5;i++){
    for(j=0;j<MAXISLCNT;j++){
      islcnta[i][j]=0;
    }
  }
  isizemax=0;
  
  for(i=0;i<MAXISLCNT;i++) {
    numisizea[i]=0;
    tsigma[i]=0.0;
    for (j=0;j<NCAP;j++){
      ct[i][j]=0; //this is the counter for capture number
    }
  }
  
  for(i=0;i<NCAP;i++) {
    avsigma[i]=0;
    avsigma2[i]=0;
    for(j=0;j<MAXISLCNT;j++) {
      sigma[j][i]=0.0;
      sigma2[j][i]=0.0;
    }
  }
  
  for(irun=0;irun < NRUN;irun++){
    //main part of the code
    printf("irun=%d\n",irun);	  
    idata=-1;
    ilogdata=0; 
    ideposit=0;
    countflag=0;
    isizemax=0;
    skip=0;
    covnext=CMIN;
    covfactor=(log(CMAX)-log(CMIN))/(NDATA);
    covfactor=exp(covfactor);
    printf("covfactor=%g\n",covfactor);
    for (i=0;i<N;i++)    {// set lattice heights to zero 
      for (j=0;j<N;j++)	{
	h[i*N+j]='0';	}    }  
    
    //initialize indexa, ipointa, nw 
    for(i=0;i<nn6;i++){
      ipointa[i]=-1;
      indexa[i]='5';
    }

    /*
    for (i=0;i<N;i++){// initialize the lists 
      for (j=0;j<N;j++){
	upnbhd(i,j);	
      }    
    }
    */
    
    // initialize base of tree (is this necessary?)
    for(i=0;i<8;i++) {   
      nw[i]=0;//this line was commented out -brad // number of atoms of type i = 0
    }
   #ifdef DEBUG  
    for(i=0;i<8;i++){
      printf("nw[%d]=%d\n",i,nw[i]);
    }
   
    printf("deprate=%g\n",deprate);
    #endif
    xx=CMAX*N*N;
    totaldropped=(int) (xx+0.5);
    deprate=N*N;
    #ifdef DEBUG
    printf("xx=%g totaldropped=%d\n",xx,totaldropped);
    printf("about to start main loop deprate=%g diffrate=%g \n",deprate,diffrate);
   #endif 
    diffq=diffrate/NDIR;
    
    
    while (ideposit <= totaldropped+1){

#ifdef DEBUGPIC
	  for(int il=0;il<N;il++){
	    for(int jl=0;jl<N;jl++){
	    printf("%c",h[il*N+jl]);
	}
	printf("\n");
      }
#endif
#ifdef DEBUGCHECK
      debugcheck();
#endif

      //printf ("%d, %f \n",ideposit, xx);
      //      totaldiff = tree[0][0];
      double rate0=nw[0]*diffrate;
      double ratebarr=nw[1]*diffrate;
      sum02= nw[0]+nw[2];//type 2 is ballistic (fresh) monomers - only have 1 direction
      //      printf("in main sum02=%d\n",sum02);
 #ifdef DEBUG
       printf("in main nw[0]=%d nw[1]=%d nw[2]=%d sum02=%d\n",nw[0],nw[1],nw[2],sum02);
       for(int il=0;il<nw[2];il++){
	 int isite=list[2][il];
	 int dir=isite%NDIR;
	 i=(isite/NDIR)/N;
	 j=(isite/NDIR)%N;
	 printf("ballistic at i=%d j=%d dir=%d\n",i,j,dir);
       }
  #endif     
      //type 0 is regular (no barrier) monomers - have 4 directions
      //both type 0 and type 1 are assumed to have the same rate-per-direction (type 0 has 4 directions)
      //type 1 is monomer which is one-step away in any of the 4 directions from another particle (island or monomer)
      // it has barrier in any of the "sticking" directions, and no barrier in the others
      sumtot=sum02+nw[1]*barrier;//barrrier = equal reduction factor for sticking direction
      totaldiff=sumtot*diffrate;
      //      totaldiff=(nw[0]+nw[2]+nw[1]*barrier)*diffrate;
      totalRate = deprate + totaldiff;
      xtemp=totalRate*uni();
      if (xtemp <= deprate){
	cov=(float)ideposit/xnsq;
	if(cov >= covnext){ 
	  cova[ilogdata]=cov; 
	  takelogdata(ilogdata);
	  ilogdata++;
	  //    #ifdef DEBUG
	    printf("ilogdata=%d cov=%g\n",ilogdata,cov);
	  //    #endif
	  covnext=covnext*covfactor;
	}
	
	/* about to deposit */
	/* test for time to take isd data (5 intervals) */
	if (cov >= (idata+2)*dcov){ 
	  idata++;			    
	  //printf("idata=%d\n",idata);
	  takedistdata(); 
	}
	depositmonomer();
	ideposit++;
  #ifdef DEBUG
		  printf("ideposit=%d\n",ideposit);
  #endif    
      }
      else		  
	{	 
    diffuse();
	}
    }/* end of deposition loop */
  }/* loop over irun */
  
#ifdef PICTURE
  printf("now to save data\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      fprintf(picdata,"%c\n",h[i*N+j]);
    }}
#endif
  
  //    for(icdata=0;icdata<5;icdata++)      capnum(icdata);
  
  /* now output results */
  /* output some densities first */
  
 output: fprintf(densdata,"N=%d d/f=%g barrier=%g CMAX=%g\n",N,diffrate,barrier,CMAX);
        fprintf(dist1,"N=%d d/f=%g barrier=%g CMAX=%g\n",N,diffrate,barrier,CMAX);  
        fprintf(dist2,"N=%d d/f=%g barrier=%g CMAX=%g\n",N,diffrate,barrier,CMAX);
        fprintf(dist3,"N=%d d/f=%g barrier=%g CMAX=%g\n",N,diffrate,barrier,CMAX); 
        fprintf(dist4,"N=%d d/f=%g barrier=%g CMAX=%g\n",N,diffrate,barrier,CMAX);
        fprintf(dist5,"N=%d d/f=%g barrier=%g CMAX=%g\n",N,diffrate,barrier,CMAX); 
        fprintf(picdata,"N=%d d/f=%g barrier=%g CMAX=%g\n",N,diffrate,barrier,CMAX);    
  printf("***we are at output and cov=%g\n",cov);
  
  fprintf(densdata,"cov N1 N2p\n");
  for (i=0;i<NDATAp1;i++) {
    covv=cova[i];
    x1=mondensa[i]/(xnsq*NRUN); /* monomer density */
    xx=densa[i]/(xnsq*NRUN); /* island density */
    double yy=xnsq*NRUN;
    fprintf(densdata,"%g  %g       %g      \n",covv,x1,xx);
  }
  //printf("printed density data\n");
  /* now to output island-size distributions */   
  
  for(i=0;i<5;i++){
    //printf("*\n");
    cov=(i+1)*CMAX*.2;
    //printf("cov\n");
    if(i==0) fprintf(dist1,"s  Ns     s/S      f(u)    sav  \n");
    //printf("1\n");
    if(i==1) fprintf(dist2,"s  Ns     s/S      f(u)    sav  \n");
    //printf("2\n");
    if(i==2) fprintf(dist3,"s  Ns     s/S      f(u)    sav  \n");
    //printf("3\n");
    if(i==3) fprintf(dist4,"s  Ns     s/S      f(u)    sav  \n");
    //printf("4\n");
    if(i==4) fprintf(dist5,"s  Ns     s/S      f(u)    sav  \n");
    //printf("5\n");
    //N=%d diffrate=%g CMAX=%g\n",N,diffrate,CMAX);  
    
    //printf("finished ifs \n");
    xx=0; yy=0; zz=0; sav=0;
    for(j=2;j<=isizemax;j++){
      xx=xx+islcnta[i][j];
      yy=yy+j*islcnta[i][j];
    }// end of first j loop
    //printf("ran middle j loop\n");
    sav=yy/xx;
    //printf("start of j loop \n");
    for(j=1;j<=isizemax;j++)
      {
	//printf("j\n");
	xx=islcnta[i][j]/(xnsq*NRUN);
	yy=j/sav;
	zz=xx*sav*sav/cov;
	
	if(xx!=0)  
	  {
	    //printf("yy\n");
	    if(i==0) fprintf(dist1,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	    if(i==1) fprintf(dist2,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	    if(i==2) fprintf(dist3,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	    if(i==3) fprintf(dist4,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	    if(i==4) fprintf(dist5,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	  }
	
      }/* end of loop over j */
    //printf("end of j loop \n");
  }/* end of loop over i */
  //printf("end of main loop \n");
}/* end of main */


void debugcheck(){
  int i,j,isite,isitedir,dir,index,iidir,iloop,itype,nn;
  //we are going to go thru the entire lattice and for each occupied site we are 
  //going to look in all 4 directions, check for barrier, check to make sure on lists etc.
  //CHECK BELOW ADDED 1/8/2018 JGA
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(h[i*N+j]=='0'){
      isite=i*N+j;
      for(dir=0;dir<NDIR;dir++){
	isitedir=isite*NDIR+dir;
	if(indexa[isitedir]!='5'){printf("empty site has monomer!!!\n");abort();}
	continue;
      }}
  //CHECK ABOVE ADDED 1/8/2018 JGA
      isite=i*N+j;
      int iflag=0;//check for ballistic
      for(dir=0;dir<NDIR;dir++){
	isitedir=isite*NDIR+dir;
	if(indexa[isitedir]=='2'){iflag=1;iidir=dir;break;}
      }
      if(iflag==1){	//so iflag=1 (one of directions is ballistic)
	if(icount(i,j,iidir)!=2){
	  printf("problem with ballistic\n");abort();
	}
	for(iloop=1;iloop<=3;iloop++){
	  isitedir++;
	  if(indexa[isitedir]!='5'){
	    printf("another problem with ballistic\n");abort();
	  }
	  return;//done with ballistic checks
	}
      }
      else{
	for(dir=0;dir<NDIR;dir++){
	  index=icount(i,j,dir);//this will be either 0,1,2 or -1
	  if(index==5)continue;
	  isitedir=isite*NDIR+dir;
	  if((int)(indexa[isitedir]-'0')!=index){
	    printf("i=%d j=%d dir=%d indexa=%d is not equal to index=%d\n",i,j,dir,(int)(indexa[isitedir]-'0'),index);
	    abort();
	  }
	  if(ipointa[isitedir]==-1){
	    printf("error2 i=%d j=%d dir=%d indexa=%d is not equal to index=%d\n",i,j,dir,(int)(indexa[isitedir]-'0'),index);
	    abort();
	  }
	  if(ipointa[isitedir] >= nw[index]){
	    printf("error3 i=%d j=%d dir=%d indexa=%d is not equal to index=%d\n",i,j,dir,(int)(indexa[isitedir]-'0'),index);
	  }
	}//for dir=0,NDIR
      }//if iflag==0
    }}//for i,j

  //CHECK BELOW ADDED 1/8/2018 JGA
  //now let's check the lists (0,1,2)
  for(itype=0;itype<=2;itype++){
    for(nn=0;nn<nw[j];nn++){
      isitedir=list[itype][nn];
      if(ipointa[isitedir]!=nn || (int)(indexa[isitedir]-'0')!=itype){
	printf("problem with nw lists\n");abort();}
    }
  }
  //CHECK ABOVE ADDED 1/8/2018 JGA

}//end of debugcheck();


void takelogdata (int ilogdata)
{
  int i,j,k,numcl,icntr;
  numcl=0;
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      imark[i*N+j]='0';
    }}
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(h[i*N+j]>'0' && imark[i*N+j]=='0') 
	{
	  numcl++;isize=1;imark[i*N+j]='1';
	  check(i,j,numcl);// this is recursive call
	  if(isize==1)mondensa[ilogdata]=mondensa[ilogdata]+1;
	  if(isize>=2)densa[ilogdata]=densa[ilogdata]+1;
	  //	  if(isize>=3)dens3a[ilogdata]=dens3a[ilogdata]+1;
	  //	  if(isize>=7)dens7a[ilogdata]=dens7a[ilogdata]+1;
	}
    }/* end of j  =0-N loop */
  }/* end of i  =0-N loop */
  mondensqa[ilogdata]=mondensqa[ilogdata]+nw[0]*nw[0];
  nwa[ilogdata]=nwa[ilogdata]+nw[0]; // walker density
  //  numcapturea[ilogdata]=numcapturea[ilogdata]+numcapturetot;
  //  numnucla[ilogdata]=numnucla[ilogdata]+numnucltot;
}/* end of takelogdata */


void takedistdata ()
{
  int i,j,k,numcl;
  double x,y;
  for(i=0;i<MAXISLSIZE;i++) {
    numisizea[i]=0;}/* I added(jga) this here */
  
      numcl=0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      imark[i*N+j]='0';
    }}

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(h[i*N+j]>'0' && imark[i*N+j]=='0'){
	numcl++;isize=1;imark[i*N+j]='1';
	check(i,j,numcl);
	islcnta[idata][isize]=islcnta[idata][isize]+1;
	numisizea[isize]=islcnta[idata][isize];
	if(isize > isizemax)isizemax=isize;
      }
    }}/* i,j loop over system */
  printf("in takedistdata isizemax=%d idata=%d\n",isizemax,idata);
  /* I added(jga) this */ 
  x=0;
  y=0;
  for(i=2;i<=isizemax;i++){
    y=y+islcnta[idata][i];
    x=x+i*islcnta[idata][i];}
  sava[idata]=x/y;
  printf("sava[%d]=%g\n",idata,sava[idata]);
  /* I added(jga) this */ 

#ifdef FRACTAL
        fractal();//calculate radius of gyration etc. 
	printf("just called fractal\n");
#endif 
}/* end of takedistdata */



void check(int i, int j, int numcl)
{
  int inb,jnb,knb,idir;
  for(idir=0;idir<NDIR;idir++){
    inb=nbxa[i][idir];
    jnb=nbya[j][idir];
    //    if(inb<0 || inb >= N){printf("in check inb=%d\n",inb);abort();}
    //    if(jnb<0 || jnb >= N){printf("in check jnb=%d\n",jnb);abort();}
    if(h[inb*N+jnb] >'0' && imark[inb*N+jnb]=='0'){
      isize++;imark[inb*N+jnb]='1';
      check(inb,jnb,numcl);}
  }
}


void depositmonomer(void)
{
  int i,j,k,dir,isite,ndir;
  unsigned char hh;  
  
  i = N*uni();
  if(i==N)i=N-1;
  j = N*uni();
  if(j==N)j=N-1;
  //  dir=uni()*NDIR;
  isite=(i*N+j)*NDIR;//this dir 0
  //  abort();
  hh=h[i*N+j];
  if(hh=='1')
    {
      search(i,j);
      i=iflat;//iflat,jflat, and isiteflat are globals
      j=jflat;
      isite=isiteflat;  
  h[i*N+j]++;
  upnbhd(i,j);
  return;
}
  h[i*N+j]++;
#ifdef DEBUG
  printf("we are in deposit  i=%d j=%d ideposit=%d dir=%d h=%c\n",i,j,ideposit,dir,h[i*N+j]);
#endif
  //add(2,isite);//assume its ballistic (type 2) travelling in direction dir
  add(0,isite);//assume its monomer (type 0) travelling in direction 0
  isite++;
  add(0,isite);//assume its monomer (type 0) travelling in direction 1
  isite++;
  add(0,isite);//assume its monomer (type 0) travelling in direction 2
  isite++;
  add(0,isite);//assume its monomer (type 0) travelling in direction 2
  upnbhd(i,j);
  #ifdef DEBUG
  printf("indexa=%d\n",(int)(indexa[isite]-'0'));
  printf("isite in deposition is  %d\n",isite);
  #endif
 /* for (i=0;i<4;i++)
    {
      for(j=0;j<nn6small;j++)
      {
    printf("list[%d][%d] %d\n",i,j,list[i][j]);
  }
  printf("\t");
}
  int iicount=icount(i,j,dir);
  #ifdef DEBUG
  printf("in deposit icount=%d iicount\n",iicount);
 #endif 
  if(iicount==2){upnbhdspecial(i,j);

    #ifdef DEBUG
    printf("indeposit after upnbhdspecial, nw[0]=%d nw[1]=%d nw[2]=%d sum02=%d\n",nw[0],nw[1],nw[2],sum02);
   #endif*/ 

}// end of void deposit ()  

void search(int i, int j)
{
  int inb,jnb,knb,idir;

  while(h[i*N+j]!='0')
  { 
    //idir=dir;
    idir=uni()*NDIR;
    if(idir==NDIR)idir=NDIRM;
    //    i=isitenb/N;
    //    j=isitenb%N;
    i=nbxa[i][idir];
    j=nbya[j][idir];
    //    isiteflat=i*N+j;
  }
  isiteflat=((i*N+j)*NDIR+idir);
  iflat=i;
  jflat=j;
  //randomly pick new direction for search
}/* end of void search(i,j)*/



void add(int index, int isite)
     {
     list[index][nw[index]]=isite;
     indexa[isite]=(unsigned char)index+'0';
     ipointa[isite]=nw[index];
     nw[index]=nw[index]+1;
     }/* end of add(int index, int isite)*/

int icount(int i, int j, int dir) // this determines class of particle 
{
     // here we only check for monomer (bond in any direction)
  int inb, jnb, knb, idir,nbonds,dirp,dirm,iidir;
  int icount,ibond,ichk,ibonda[6],nchk,nchk1,nchk2;
  int bdira[10];
  unsigned char hx,hnew;
  //check to make sure not blocked "bonded" in direction
    hx=h[i*N+j];
    if(hx=='0')return(5);//-1 if site is not occupied !
    inb=nbxa[i][dir];
    jnb=nbya[j][dir];
    hnew=h[inb*N+jnb];
    if(hnew =='1'){return(5);}//-1 if new site is already occupied
  // now count bonds
  for(ichk=0;ichk<NDIR;ichk++){
    if(h[nbxa[i][ichk]*N+nbya[j][ichk]]=='1'){
      return(5);//-1 if existing site has any bonds
    }
  }
  // so number of bonds is 0  check if site we are moving to has any bonds (if so, then barrier)
  //printf("skipdir=%d\n",(dir+NDIRH)%NDIR);
  for(ichk=0;ichk<NDIR;ichk++){
    if(ichk==(dir+NDIRH)%NDIR)continue;
    if(h[nbxa[inb][ichk]*N+nbya[jnb][ichk]]=='1'){
      #ifdef DEBUG
      printf("ichk=%d\n",ichk);
      #endif
      return(1);//return 1 if site we are moving to already has bonds - barrier
      //      ibonda[nbonds]=ichk;
      //      nbonds=nbonds+1;
      //      if(nbonds>=4)return(-1);
    }
  }
  //SO NO BARRIER!
  if(indexa[(i*N+j)*NDIR+dir]=='2')return(2);//if no barrier then return 2 if we are already type 2
  return(0);//if no barrier, then return 0 otherwise (ordinary monomer)
}/* end of int icount (int i, int j, in   */  

void upnbhdspecial(int i, int j)// update nbhd of i,j,k
{
 
  int inb,jnb,knb,dir1,dir2,dir3,dir4,iloop,iic,idirr;

  //    for(dir2=0;dir2<NDIR;dir2++){  update(i,j,dir2);  }// update the site itself

  /*
  for(dir2= 0;dir2<NDIR;dir2++){// update the immediate neighbors all directions
    inb=nbxa[i][dir2];
    jnb=nbya[j][dir2];
    update(inb,jnb,(dir2)%NDIR);
    update(inb,jnb,(dir2+1)%NDIR);
    update(inb,jnb,(dir2+2)%NDIR);
    update(inb,jnb,(dir2+3)%NDIR);
  }
  */

  //update the 4 n.n.nbrs. in two special directions each
  dir1=0;dir2=3;dir3=2;dir4=1;
  for(iloop=0;iloop<NDIR;iloop++){
  inb=nbxa[nbxa[i][dir1]][dir2];
  jnb=nbya[nbya[j][dir1]][dir2];
  #ifdef DEBUG
  printf("in upnbhdspecial inb=%d jnb=%d,dir3=%d dir4=%d\n",inb,jnb,dir3,dir4);
  #endif
  updatefix(inb,jnb,dir3);
  updatefix(inb,jnb,dir4);
  dir1=(dir1+1)%NDIR;
  dir2=(dir2+1)%NDIR;
  dir3=(dir3+1)%NDIR;
  dir4=(dir4+1)%NDIR;
  }
  //2nd nbrs in same direction ("the cross")
  dir1=0;
  for(iloop=0;iloop<NDIR;iloop++){
  inb=nbxa[nbxa[i][dir1]][dir1];
  jnb=nbya[nbya[j][dir1]][dir1];
  #ifdef DEBUG
  printf("in upnbhdspecial - the cross - inb=%d jnb=%d,dir1=%d\n",inb,jnb,dir1);
  #endif
  idirr=(dir1+NDIRH)%NDIR;
  updatefix(inb,jnb,idirr);
  dir1++;
  }
}//end of void upnbhdspecial(i,j,dir) */

void updatefix(int inb, int jnb,int idirr){
  int iflag,isitec,iic,indexx,ic,d1,d2,d3,i,j;
  iflag=0;
  isitec=(inb*N+jnb)*NDIR;
#ifdef DEBUG
  printf("we are in updatefix\n");
#endif  
  //  printf("indexa=%d\n",indexa[isitec+idirr]);
  //  printf("indexa=%d\n",indexa[isitec+(idirr+1)%NDIR]);
  //  printf("indexa=%d\n",indexa[isitec+(idirr+2)%NDIR]);
  //  printf("indexa=%d\n",indexa[isitec+(idirr+3)%NDIR]);
  indexx=(int)(indexa[isitec+idirr]-'0');
  if(indexx!=5){
	iic=update(inb,jnb,idirr);
    if(indexx==2 && iic==1){  
    #ifdef DEBUG
      printf("we are adding\n");
    #endif
      i=inb;
      j=jnb;
    d1=(idirr+1)%NDIR;
    d2=(idirr+2)%NDIR;
    d3=(idirr+3)%NDIR;
    ic=icount(i,j,d1);
    if(ic==1)
      add(1,((i*N+j)*NDIR+d1)); //Edited 7-Jan-2018
    else
       add(0,((i*N+j)*NDIR+d1));

     ic=icount(i,j,d2);
    if(ic==1)
      add(1,((i*N+j)*NDIR+d2));
     else
      add(0,((i*N+j)*NDIR+d2));
     ic=icount(i,j,d3);
    if(ic==1)
      add(1,((i*N+j)*NDIR+d3));
      else
      add(0,((i*N+j)*NDIR+d3));

      //add(0,isitec+(idirr+1)%NDIR);
      //add(0,isitec+(idirr+2)%NDIR);
      //add(0,isitec+(idirr+3)%NDIR);
    }//end of if indexx==2 && icc==1
  }//end of if indexx!=-1
}//end of void updatefix

int update(int i, int j, int dir)  // update i,j
{
   int isite,index,newindex;
  isite=NDIR*(i*N+j)+dir;
  index=(int)(indexa[isite]-'0');
  newindex=icount(i,j,dir);
  #ifdef DEBUG
  printf("in update, newindex=%d i=%d j=%d dir=%d\n",newindex,i,j,dir);
  #endif
  if(newindex != index){
    if(index != 5){delete(index,isite);}
    if(newindex != 5){add(newindex,isite);}
  }
  return(newindex);
}/* end of update(int i, int j, int dir) */

void delete(int index, int isite)
   {
     int ipos,endsite,endpos;
     int i,j,idir;
     idir=isite%NDIR;//isite=(i*N+j)*NDIR+idir
     j=(isite/NDIR)%N;
     i=(isite/NDIR)/N;
     #ifdef DEBUG
     printf("we are deleting i=%d j=%d idir=%d\n",i,j,idir);
     #endif
     ipos=ipointa[isite];
     endpos=nw[index]-1;
     endsite=list[index][endpos];
     if(endpos != ipos)
       {
   // move end to hole in list
   list[index][ipos]=endsite;//move end to hole in list
   ipointa[endsite]=ipos;//reset pointer for endsite to new position (hole)
       }
     indexa[isite]='5'; //indicate no walker at given site
     ipointa[isite]=-1;
     #ifdef DEBUG
     printf("type subtraciton nw[%d] %d\n",index, nw[index]);
     #endif
     nw[index]=nw[index]-1;//decrement number of walker of given type
   }/* end of delete (index, isite) */


void upnbhd(int i, int j)// update nbhd of i,j,k
{
  int inb,jnb,knb,dir1,dir2,dir3,dir4,iloop,iic;

    for(dir2=0;dir2<NDIR;dir2++){  iic=update(i,j,dir2);  }// update the site itself
  
  for(dir2=0;dir2<NDIR;dir2++){// update the immediate neighbors all directions
    inb=nbxa[i][dir2];
    jnb=nbya[j][dir2];
    iic=update(inb,jnb,0);
    iic=update(inb,jnb,1);
    iic=update(inb,jnb,2);
    iic=update(inb,jnb,3);
  }
  
  upnbhdspecial(i,j);
}//end of void upnbhd(i,j,dir) 

/*
  //update the 4 n.n.nbrs. in two special directions each
  dir1=0;dir2=3;dir3=2;dir4=1;
  for(iloop=0;iloop<NDIR;iloop++){
  inb=nbxa[nbxa[i][dir1]][dir2];
  jnb=nbya[nbya[j][dir1]][dir2];
  iic=update(inb,jnb,dir3);
  iic=update(inb,jnb,dir4);
  dir1=(dir1+1)%NDIR;
  dir2=(dir2+1)%NDIR;
  dir3=(dir3+1)%NDIR;
  dir4=(dir4+1)%NDIR;
  }
  //2nd nbrs in same direction ("the cross")
  dir1=0;
  for(iloop=0;iloop<NDIR;iloop++){
  inb=nbxa[nbxa[i][dir1]][dir1];
  jnb=nbya[nbya[j][dir1]][dir1];
  iic=update(inb,jnb,(dir1+NDIRH)%NDIR);
  dir1++;
  }
}//end of void upnbhd(i,j,dir) 
*/

void diffuse()
{
     int x,y,i,j,hxy,a,b,site,isite,isite2;
     int dir,iwalk,xi,yi,itype;
     double xtemp;
  int sum,ih;
  int xnb,ynb,idir,islmark;
  int inb,jnb;
  double rate0,ratebal,ratebarr,ratetot;
  //      sum02= nw[0]+nw[2];//type 2 is ballistic (fresh) monomers - only have 1 direction
      //type 0 is regular (no barrier) monomers - have 4 directions
      //both type 0 and type 1 are assumed to have the same rate-per-direction (type 0 has 4 directions)
      //type 1 is monomer which is one-step away in any of the 4 directions from another particle (island or monomer)
      // it has barrier in any of the "sticking" directions, and no barrier in the others
  //      sumtot=sum02+nw[1]*barrier;//barrrier = equal reduction factor for sticking direction
  xtemp=uni()*sumtot;//
  itype=1;//barrier
  #ifdef DEBUG
  printf("in diffuse sumtot=%f \n",sumtot);
  printf("in diffuse sum02= %d\n",sum02);
  #endif
  if(xtemp<=nw[0])itype=0;//regular monomer

  //else if(xtemp <= sum02)itype=2;//ballistic
     iwalk = uni()*nw[itype];
     #ifdef DEBUG
      printf("iwalk=%d and xtemp=%f and nw=%d type=%d \n", iwalk,xtemp,nw[itype],itype);
     #endif
     if (iwalk==nw[itype]) iwalk = nw[itype]-1; //fixed JGA 2/26/11
     isite  = list[itype][iwalk];//  isite=NDIR*(i*N+j)+dir;
     #ifdef DEBUG
     printf("isite %d\n",isite );
      printf("iwalk %d\n",iwalk );
      #endif
     dir   = isite % NDIR;
     #ifdef DEBUG
     printf("dir=%d and N=%d are\n",dir,N );
     #endif
       isite2 = isite /NDIR;//  isite=NDIR*(i*N+j)+dir;
     x     = isite2/N;
     y     = isite2 % N;
     //     h[x][y]=h[x][y]-1;
     h[x*N+y]='0';
     #ifdef DEBUG
     printf("x=%d and y=%d \n",x,y);
     #endif
     idiffuse++;
     #ifdef DEBUG
     printf("h  is %c\n",h[x*N+y]);
     #endif
  
     xnb=nbxa[x][dir];
     ynb=nbya[y][dir];
     int isitenb=(xnb*N+ynb)*NDIR+dir;
     //     ih=h[xnb][ynb];
     //     h[xnb][ynb]=ih+1;
     h[xnb*N+ynb]++;
     #ifdef DEBUG
     printf("we are in diffuse x=%d y=%d h=%c dir=%d -> xnb=%d ynb=%d h=%c dir=%d isitenb=%d \n",x,y,h[x*N+y],dir,xnb,ynb,h[xnb*N+ynb],dir,isitenb);
     printf("indexa=%d\n",(int)(indexa[isite]-'0'));
     printf("we are moving particle of itype=%d dir=%d\n",itype,dir);
     #endif
     //     printf("we are moving particle of type %d\n",itype);

  if(itype==2){//we are moving a ballistic monomer so just manually update the list for that monomer
    indexa[isitenb]='2';
    if(icount(xnb,ynb,dir)==2){ //it remains ballistic
    //       delete(2,isite);//deleted ballistic monomer that was at the old site
    //       add(2,isitenb);//add the ballistic monomer at the new site.
    //    int ilist=list[2][ipointa[isite];
    //this replaces delete(2,x,y,dir) add(2,xnb,ynb,dir)
  #ifdef DEBUG
      printf("we are in diffuse and verified that icount==2\n");
      #endif
    list[2][iwalk]=isitenb;//point to new site in ballistic monomer list
    ipointa[isitenb]=iwalk;//fix inverse pointer for new site
    ipointa[isite]=-1;//fix inverse pointer for old site
    indexa[isitenb]='2';//fix indexa for new site
    indexa[isite]='5';//fix indexa for old site
    //we have moved the ballistic monomer to another position, now we just need to set-up/change barriers for 2nd nbr atoms of original and new sites
    upnbhdspecial(x,y);//this skips nbr and self
    upnbhdspecial(xnb,ynb);//this skips nbr and self
  //this just takes care of other monomers it might be affecting
    return;
    }//end of if icount==2
    //so it's changing to type 1 
    indexa[isitenb]='5';
  }//end of if(itype==2) 
  //so it was (or is now) a regular monomer of type 0 or type 1
  //if(icount(xnb,ynb,dir)==1 && itype==2){nw[2]=nw[2]+1;}
  upnbhd(x,y);
  upnbhd(xnb,ynb);
  //       delete(itype,isite);//deleted ballistic monomer that was at the old site
  //     add(itype,isitenb);//add the ballistic monomer at the new site.
}/* end of diffuse () */



    /* uni.c */
    /* 3D triangulation program written by Isabel Beichl */
    /* based on an algorithm designed by Isabel Beichl & Francis Sullivan */
    /*    National Institute of Standards & Technology */
    /*    Gaithersburg, MD 20899  */ 
     
    static unsigned long count = 0;
    #define BIGCOUNT 10000000    /* how many to do without re-initializing */

    static unsigned long m1 = 32767;
    static long int mb[607];
    /*=
    {
            30788, 23052, 2053, 19346, 10646, 19427, 23975,
            19049, 10949, 19693, 29746, 26748, 2796, 23890,
            29168, 31924, 16499
    };
    */
    static int mdig = 32;
    static unsigned long m2 = 256;
    static int iran = 272;
    static int jran = 606;




    double suni(jseed)
          unsigned long jseed;
    {
            long int j0, j1, k0, k1;
            double uni();

    /*      printf(" suni %ld\n", jseed);*/
            m1 = poww2(mdig-2) - 1;     /* avoid overflow if m1 is full size */
            m1 += m1;
            m1++;
       /* printf(" m1 %lu, m2 %lu, mdig %d, jseed %u\n", m1, m2, mdig, jseed); */
            m2 = poww2((int)(mdig/2));
       /* printf(" m1 %lu, m2 %lu, mdig %d, jseed %u\n", m1, m2, mdig, jseed); */
            jseed %= m1;                    /* jseed should less than m1 */
            if ((jseed & 1) == 0)           /* jseed should be odd */
                    jseed--;
            k0 = 9069 % m2;                 /* simple congruential generator */
            k1 = 9069 / m2;                 /* the fanciness avoids overflow */
            j0 = jseed % m2;
            j1 = jseed / m2;
            for (iran = 0; iran < 607; iran++)
            {       jseed = j0*k0;
                    j1 = (jseed/m2 + j0*k1 + j1*k0) % (m2/2);
                    j0 = jseed % m2;
                    mb[iran] = j0 + m2*j1;
        /* printf("%2d %10u\n", i, mb[i]); */
            }
            iran = 272;
            jran = 606;
            return uni();
    }

    double uni()
    {
            long int k;

            k = mb[iran] - mb[jran];
            if (k < 0)
                    k += m1;
       /* printf(" In UNI -- k = %ld\n",k); */
            if (++count >= BIGCOUNT)
            {       count = 0;
                    suni(k);
            }
            else
            {       mb[jran] = k;
                    if (--iran < 0)
                            iran = 606;
                    if (--jran < 0)
                            jran = 606;
            }
       /* printf("%lf\n", (double)k/(double)m1); */
       /* putchar(k%2 ?'+':'-');*/
            return ((double)k/(double)m1);
    }

unsigned long poww2(int j)
    {
            unsigned long x = 1;
       /* printf("poww2--j= %d\n",j); */

            while (j--)
                    x *= 2;
          /* printf("poww2--x= %lu\n",x); */
            return (unsigned long) x;
    }
    /* end of uni2.c */














/*
void capnum(int icdata)
{
  int i,irun;
  double mondens,denom,xncubdcovn1,msigma,dmsigma;
  double density,scaledsigma,dcov,densi;
  
  dcov=CMAX*0.2;
  xncubdcovn1 = diffq*islcnta[icdata][1]*.05*(icdata+1)*dcov/NRUN;
  density = 0.0;
  msigma = 0.0;
  printf("in capnum icdata=%d isizemax=%d\n",icdata,isizemax);
  for(i=2;i<=isizemax;i++){
  densi=islcnta[icdata][i]/xnsq;
  denom=xncubdcovn1*densi;
  density = density + densi;
  if(denom>0) {//I added(jga) brackets here 
    tsigma[i]=ct[i][icdata]/denom;
  //  msigma = msigma + tsigma[i]*densi;
  //    avsigma[icdata]=avsigma[icdata] + dmsigma;
  //  avsigma[icdata]=tsigma[i]*densi;
        sigma[i][icdata]=sigma[i][icdata]+tsigma[i];
    //    sigma2[i][icdata]=sigma2[i][icdata]+tsigma[i]*tsigma[i];
  }//I added(jga) brackets here 
  }

  //    dmsigma=msigma/density;
  //    avsigma2[icdata]=avsigma2[icdata] + (dmsigma*dmsigma);
    
    printf("leaving capnum\n");
}//end of void capnum(int icdata) 
*/
