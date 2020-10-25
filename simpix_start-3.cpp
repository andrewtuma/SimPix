// simple example of using ROOT libraries in a C++ program with graphics
// and use of TASImage class

#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"
#include "TSystem.h"


#include "assert.h"

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <time.h>
using namespace std;

#define NUM_TRIES_FACTOR 2
#define NUM_SUCCESS 5000
#define T_STEP .01

struct RGB {
	float r;
	float g;
	float b;
} ;

RGB* getRGB(int code){
	RGB* c = new RGB;
	c->b = ((code) & 0xff)/255.;
	c->g = ((code >> 8) & 0xff)/255.;
	c->r = ((code >> 16) & 0xff)/255.;
	return c;
}

long getInt(RGB* c){
	long i = 0xff;
	i *= 256;
	i += (int)(c->r * 255);
	i *= 256;
	i += (int)(c->g * 255);
	i *= 256;
	i += (int)(c->b * 255);
	return i;
}	

float distRGB(RGB* c1, RGB* c2){
	float d = 0;
	float dr, dg, db;
	dr = c1->r - c2->r;
	dg = c1->g - c2->g;
	db = c1->b - c2->b;
	return sqrt(dr*dr + dg*dg + db*db);
}

float dist(RGB* c1, RGB* c2){
	return distRGB(c1,c2);
}

void doSwap(RGB** img, int p1, int p2){
	RGB* temp = img[p1];
	img[p1] = img[p2];
	img[p2] = temp;
}

float swapDEnergy(RGB** img1, RGB**img2, int p1, int p2){
	double d1, d2;
	d1 = dist(img1[p1],img2[p1]) + dist(img1[p2],img2[p2]);
	d2 = dist(img1[p1],img2[p2]) + dist(img1[p2],img2[p1]);
	return d2 - d1;
}

int getRandPix(int size){
	return lrand48()%size;
}

int doMetropolisSwap(RGB** img1, RGB** img2, int size, float T, double& totalDist){
	int p1 = getRandPix(size);
	int p2 = getRandPix(size);

	float d = swapDEnergy(img1,img2,p1,p2);

	if(d<0 || drand48()< exp(-d/T)){
		doSwap(img2, p1,p2);
		totalDist += d;
		return 1;
	}	
	return 0;
}

bool doMetropolisStep(RGB** img1, RGB** img2, int size, float& T, double& totalDist){
	int numSuccess = 0;
	for (int i=0; i<size*NUM_TRIES_FACTOR; i++){
		numSuccess += doMetropolisSwap(img1, img2, size, T, totalDist);
		
	}
	if (numSuccess >= NUM_SUCCESS){
			cout << "totalDist = "<< totalDist << "\tnumSuccess = " << numSuccess << "\tT = "<< T << endl;
			return true;
		}
	return false;
}

void doAllMetropolis(RGB** img1, RGB** img2, int size, float T, double& totalDist){
	bool happy = true;
	double T0 = T;
	while (happy && T>0){
		happy = doMetropolisStep(img1, img2, size, T, totalDist);		
		T -= T_STEP*T0;
	}
}

int main(int argc, char **argv){

  if (argc<3) {
    cout << "Usage: simapix_start image1 image2 <output=out.png>" << endl;
    return 0; 
  }
  TString fsrc=argv[1];
  TString ftgt=argv[2];
  TString fout;
  argc>3 ? fout = argv[3] : fout="out.png";
  cout << "Reading images: source= " << fsrc << " target= " << ftgt << endl;
  cout << "Output= " << fout << endl;

  TApplication theApp("App", &argc, argv);

  // create image objects
  TASImage *src = new TASImage(fsrc.Data());
  TASImage *tgt = new TASImage(ftgt.Data());
  TASImage *out = new TASImage(*src); // start with copy of source

  // Test image geometry, exit if they are not the same dimensions
  assert ( src->GetWidth() == tgt->GetWidth() && src->GetHeight() == tgt->GetHeight() );
  cout << "Pixel Geometry: " << src->GetWidth() << " x " << src->GetHeight() << endl;
  Long_t numPix=src->GetWidth()*src->GetHeight();

  // *** The work happens here
  // access the pixels for the output image 
  // each pixel is a 32-bit word, 1 byte each for (alpha,red,green,blue)
  // don't touch alpha (bits 31:28)
  UInt_t *outPix = out->GetArgbArray();  

	/*
  // examples of pixel manipulations 
  for (int i=0;i< numPix; i++){
    //  outPix[i]&=0xff00ffff;  // turn off red
    outPix[i]&=0xffff00ff;  // turn off green
    //  outPix[i]&=0xffffff00;  // turn off blue
    //  cout << hex << outPix[i]<<endl;  // print pixel values in hex
  }
  // flip the image
  for (int i=0;i< numPix/2; i++){
    unsigned pxl=outPix[i];
    outPix[i]=outPix[numPix-i-1];
    outPix[numPix-i-1]=pxl;
  }
	*/
	
  // *************************

	UInt_t *inPix = src->GetArgbArray();
	UInt_t *tgtPix = tgt->GetArgbArray();

	srand48(time(NULL));

	//build the RGB arrays
	int width = src->GetWidth(); 
	int height = src->GetHeight();

	RGB** img1 = new RGB*[numPix]; //the target image
	RGB** img2 = new RGB*[numPix]; //the starting image, scrambled
	double totalDist = 0;
	for (int i = 0; i< numPix; i++){
			img1[i] = getRGB(tgtPix[i]);
			img2[i] = getRGB(inPix[i]);
			totalDist += dist(img1[i],img2[i]);
	}
	
	double T = 5;
	
	doAllMetropolis(img1,img2,numPix,T,totalDist);

	for (int i =0; i<numPix; i++){
		outPix[i] = getInt(img2[i]);
		//cout << getInt(img2[i]) << endl;
	}


  // print the results
  TCanvas *c1 = new TCanvas("c1", "images", 640, 480);
  c1->Divide(2,2);

  c1->cd(1);
  c1->Draw();
  src->Draw("X");
  c1->cd(2);
  tgt->Draw("X");
  c1->cd(3);
  out->Draw("X");
  c1->Print("collage.png");
  
  // save the new image
  out->WriteImage(fout.Data());

  // comment out the lines for running in batch mode
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();

}
