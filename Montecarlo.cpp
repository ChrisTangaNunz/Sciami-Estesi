#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TGraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include <bits/stdc++.h>
#include <algorithm>
#include <cmath>
#include <TAxis.h>
#include <iomanip>
#include <TBox.h>
#include <TPaveText.h>
#include <TRandomGen.h>
using namespace std;




void coincidenza_doppia12(int n, double Lu1, double La, double spessore, double Lu2, double h ){

    double S1=La*Lu1;
    double S2=La*Lu2;
    TRandomMixMax randomGenerator;
    /*
    TH1F *histogram= new TH1F("histogram","titlevero.c_str()",100, -1,2);
  
    new TCanvas;
    TH2F * ProbabilityDensity= new TH2F("Probability Density2","Density of Probability2",100,-3.14/2,3.14/2,100,-1,2);
   
    */
   
    double tentativi=0, successi=0, brutti=0;
    for(int i=0; i<n; i++){
        
        //genero un punto nel piano x1*y1
        double x1=randomGenerator.Uniform(La)-La/2;
        double y1=randomGenerator.Uniform(Lu1)-Lu1/2; 
        double  phi=randomGenerator.Uniform(),z1=0;
        double theta=randomGenerator.Uniform();
        double check_theta=randomGenerator.Uniform();
        brutti++;
        if(check_theta>pow(cos(theta),2)) continue;
        tentativi++;
        double x2=x1-h*tan(theta)*cos(phi);
        double y2=y1-h*tan(theta)*sin(phi);
        double z2=z1-h;
        if(x2<=La/2 && x2>=-La/2 && y2<=Lu2/2 && y2>=-Lu2/2){
            successi++;
        }
    }
    double efficienza_geom=successi/tentativi;
    printf("successi = %lf \n", successi);
    printf("tentativi = %lf \n", tentativi);
    printf("effgeom = %lf \n", efficienza_geom);
    //printf("efficienza fisica = %lf \n", efficienza_geom);
 

}

int Montecarlo(){

    coincidenza_doppia12(100, 52, 39.5, 1.5, 52.2, 16.3); //n iterazioni, Lunghezza scintillatore 1, Larghezza scintillatore 2, spessore scintillatore, lunghezza scintillatore 2, altezza (distanza tra i due scinitllatori)

    return 0; 

}
