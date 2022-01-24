#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
//#include "savehts2.h"
/*
 * GL01Hello.cpp: Test OpenGL/GLUT C/C++ Setup
 * Tested under Eclipse CDT with MinGW/Cygwin and CodeBlocks with MinGW
 * To compile with -lfreeglut -lglu32 -lopengl32
 */
//#include <windows.h>  // for MS Windows
#include "/Users/jamar2/GL/glut.h"  // GLUT, include glu.h and gl.h
//#include <GLUT/glut.h>// Header File For The GLUT Library
//#include <OpenGL/gl.h>// Header File For The OpenGL32 Library
#include "/Users/jamar2/GL/gl.h"  // GLUT, include glu.h and gl.h
//#include <OpenGL/glu.h>// Header File For The GLu32 Library
#include "/Users/jamar2/GL/glu.h"  // GLUT, include glu.h and gl.h
float rr[11]={1.0,1,0.0,0,1.0,0,0,0.3,0.5,0.5,2.0};
float GG[11]={1.0,0,0,0,0.5,0.5,0.9,0.6,0.0,0.8};
float bb[11]={1.0,1,1,1,0,1,0,0.5,0.1,0.5};
GLUquadricObj * quadObj;
void fullcg_(float *i, float *j, int *hh,float *size);
void sgisave();

#define N 500
#define wsize 500
int h[N][N];
float dx=1.0/(float)N;
float  scale=2.*0.05/(N/40);
 
/* Handler for window-repaint event. Call back when the window first appears and
   whenever the window needs to be re-painted. */
void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set background color to black and opaque
  glClear(GL_COLOR_BUFFER_BIT);         // Clear the color buffer (background)
	glColor3f(0.0f,0.0f,1.0f); //blue color
    glPointSize(1);//set point size to 10 pixels
    glEnable(GL_POINT_SMOOTH);
	glTranslatef(-1.0f, -1.0f, 0.0f);
  float xx,yy;  int iloop;	int i,j;
float size= 2.*0.025/(N/40);
// printf("about to loop over sites \n");
 int iflag=0;float xiim,xjjm;
	for (i = 0; i < N; i++){		for (j = 0; j < N; j++){
	    //	    printf("i=%d j=%d\n",i,j);
	    //			float xii = i + 0.5*j;//for triangular lattice
	    //			float xjj = 0.5*j*sqrt(3.0);//for triangular lattice
	    float xii = i;//for square lattice
	    float xjj = j;//for square lattice
			int ih = h[i][j];
			if (ih == 0)				continue;		    
			if(iflag==0){xiim=xii;xjjm=xjj;iflag=1;}
			if(xii > N)    xii = xii - N;			
			//converting to 1 unit size//
						xii = 2*xii / (N);
						xjj = 2*xjj/ (N);
			/*


   xii = xii*scale; xjj = xjj*scale;
   ih=5*ih;
   printf("about to call fullcg i=%d j=%d xii=%f xjj=%f ih=%d size=%f\n",i,j,xii,xjj,ih,size);
   //   size=1;
   fullcg_(&xii,&xjj,&ih,&size);
			*/

   			glBegin(GL_POINTS); //starts drawing of points
			glVertex2f(xii,xjj);
			glEnd();//end drawing of points
			glFlush();  // Render now
   
		}}	

	

  /*  glVertex2f(-0.5f, -0.5f);    // x, y
  glVertex2f( 0.5f, -0.5f);
  glVertex2f( 0.5f,  0.5f);
  glVertex2f(-0.5f,  0.5f);
  */
  
	
	/*for(iloop=0;iloop<2;iloop++){
    if(iloop==0)  glColor3f(1.0f, 0.0f, 0.0f); // Red
    else   glColor3f(0.0f, 0.0f, 1.0f); // Blue
    printf("in display dx=%f\n",dx);
  glBegin(GL_QUADS);              // Each set of 4 vertices form a quad
   xx=-dx+(float)(iloop)*0.5;   yy=-dx;
  //  glVertex2f(-0.5f, -0.5f);    // x, y
  glVertex2f(xx, yy);    // x, y
  xx+=2*dx;
  glVertex2f(xx, yy);    // x, y
  //  glVertex2f( 0.5f, -0.5f);
  yy=dx;
  glVertex2f(xx, yy);    // x, y
  //  glVertex2f( 0.5f,  0.5f);
  xx-=2*dx;
  glVertex2f(xx, yy);    // x, y
  //  glVertex2f(-0.5f,  0.5f);
  glEnd();
  }
 */
  //glFlush();  // Render now
		sgisave();
}


void readfile(){
	
  //  	FILE *file1 = fopen("filejga.txt", "r");
    	FILE *file1 = fopen("readfile", "r");
	//	FILE *file1 = fopen("file40", "r");
	printf("File read\n");
	int ii;
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			fscanf(file1,"%d\n",&ii);
			h[i][j]=ii;
		}}
	
	fclose(file1);
	
	
}
 
/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
  printf("dx=%f\n",dx);
	readfile();
  glutInit(&argc, argv);                 // Initialize GLUT
  printf("about to initialize window size\n");
  glutInitWindowSize(wsize, wsize);   // Set the window's initial width & height
  printf("just did initialize window size\n");
  glutInitWindowPosition(0, wsize); // Position the window's initial top-left corner
  printf("just initialized window position\n");
  glutCreateWindow("OpenGL Setup Test"); // Create a window with the given title
  printf("just created window\n");
  glutDisplayFunc(display); // Register display callback handler for window re-paint
  glutMainLoop();           // Enter the event-processing loop
  return 0;
}


#define PicSize wsize  //1024?
#define winsize wsize
float hff[winsize][winsize];
unsigned char outbuf[PicSize];
void putbyte(FILE *outf,unsigned char val)
{
  unsigned char buf[1];
  buf[0]=val;
  fwrite(buf,1,1,outf);
}

void putshort(FILE *outf, unsigned short val){
  unsigned char buf[2];
  buf[0]=(val>>8);
  buf[1]=(val>>0);
  fwrite(buf,2,1,outf);
}

static int putlong(FILE *outf,unsigned long val)
{
  unsigned char buf[4];
  
  buf[0] = (val>>24);
  buf[1] = (val>>16);
  buf[2] = (val>>8);
  buf[3] = (val>>0);
  return fwrite(buf,4,1,outf);
}

void sgisave(){
  int nn=winsize;
  int ih, i,j,x,y,z;
  FILE *fsave;
  char chuck,strB[10],strA[10]=".sgi";
  char iname[80]="slkmc-hts-";
  //  static int numb=1;
int idisplay=0;
  sprintf(strB,"%d",idisplay);
  strcat(iname,strB);
  strcat(iname,strA);
  fsave=fopen(iname,"wb");
  idisplay++;
  if(!fsave) {
    printf("sgiimage: can't open output file\n");
    fprintf(stderr,"sgiimage: can't open output file\n");
    exit(1);
  }

  putshort(fsave,474);       /* MAGIC                       */
  putbyte(fsave,0);          /* STORAGE is VERBATIM         */
  putbyte(fsave,1);          /* BPC is 1                    */
  putshort(fsave,3);         /* DIMENSION is 2 BW, 3 color             */
  putshort(fsave,PicSize);    /* XSIZE                       */
  putshort(fsave,PicSize);    /* YSIZE                       */
  /*    putshort(fsave,IXSIZE);   XSIZE                       
	putshort(fsave,IYSIZE);   YSIZE                       */
  putshort(fsave,3);        //1 BW, 3 color picture
  
      putlong(fsave,0);          /* PIXMIN is 0                 */
    putlong(fsave,255);        /* PIXMAX is 255               */
    for(i=0; i<4; i++)putbyte(fsave,0);      /* DUMMY 4 bytes       */
    /*    strcpy(iname,"No Name");*/
    /*    strcpy(iname,strPic);*/
    fwrite(iname,80,1,fsave);  /* IMAGENAME           */
    putlong(fsave,0);          /* COLORMAP is 0       */
    for(i=0; i<404; i++)putbyte(fsave,0);    /* DUMMY 404 bytes     */
        for(z=0; z<3; z++) {//z<3 for color
                      if(z==0)  glReadPixels (0,0,nn,nn,GL_RED,GL_FLOAT,*hff);
                      if(z==1)      glReadPixels (0,0,nn,nn,GL_GREEN,GL_FLOAT,*hff);
                if(z==2)  glReadPixels (0,0,nn,nn,GL_BLUE,GL_FLOAT,*hff);

  for(x=0; x<nn; x++) {
    for(y=0; y<nn; y++) {
            outbuf[y] = (int)(255.*hff[x][y]);
	    //            if(h[x][y]==0)outbuf[y]=255;else outbuf[y]=0;
    }
//        fwrite(&onechar,1,1,fsave);
    fwrite(outbuf,nn,1,fsave);
	//       onechar = (int)(255.*hg[x][y]);
	//        fwrite(&onechar,1,1,fsave);
	//       onechar = (int)(255.*hb[x][y]);
	//        fwrite(&onechar,1,1,fsave);
  }
  }
    fclose(fsave);
    printf("picture saved idisplay=%d\n",idisplay);
}// end of sgisave()
