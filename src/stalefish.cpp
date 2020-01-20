#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include "tools.h"
using namespace cv;
using namespace std;


// *** GLOBALS ***

vector<frameData> F;
int I;
Mat img;

// *** DISPLAY ***
void onmouse(int event, int x, int y, int flags, void* param) {
    Point pt = Point(x,y);
    if(event==CV_EVENT_LBUTTONDOWN){ F[I].P.push_back(pt); }
    img = F[I].frame.clone();
    circle(img,pt,5,Scalar(0,0,255),1);
    for(int i=1;i<F[I].P.size();i++){
        line(img,F[I].P[i-1],F[I].P[i],Scalar(255,0,0),1);
    }
    for(int i=0;i<F[I].P.size();i++){
        circle(img,F[I].P[i],5,Scalar(255,0,0),1);
    }
    if (F[I].P.size()){
        line(img,F[I].P[F[I].P.size()-1],pt,Scalar(255,0,0),1);
    }
    for(int i=1;i<F[I].fitted.size();i++){
        line(img,F[I].fitted[i-1],F[I].fitted[i],Scalar(0,255,0),1);
    }
    line(img,F[I].axis[0],F[I].axis[1],Scalar(0,0,255),1);
    for(int i=0;i<F[I].origins.size();i++){
        line(img,F[I].origins[i],F[I].tangents[i],Scalar(0,255,255),1);
    }
    stringstream ss;
    ss<<"Frame: "<<I+1<<"/"<<F.size()<<", "<<"Poly order: "<<F[I].polyOrder<<", Bins: "<<F[I].nBins;
    putText(img,ss.str(),Point(30,30), FONT_HERSHEY_SIMPLEX, 0.8, Scalar(0,0,0), 1, CV_AA);
    imshow("StaleFish", img);
}

// *** MAIN PROGRAM ***

int main(int argc, char** argv){

    double lenA = 10.;
    double lenB = 50.;
    I = 0;
    for (int i=1;i<argc;i++){
        char* imageName = argv[i];
        Mat frame = imread(imageName, IMREAD_COLOR );
        if( frame.empty()){ cout <<  "Could not open or find the image" << std::endl; return -1;}
        F.push_back(frameData(frame.clone()));
    }
    namedWindow("StaleFish",1);
    setMouseCallback("StaleFish", onmouse, &img);

    // *** MAIN LOOP ***

    while(1) {
        onmouse(CV_EVENT_MOUSEMOVE,0,0,0,NULL);
        char k = waitKey(0);
        switch(k){
            case('x'):{
                F[I].updateFit();
                F[I].refreshBoxes(lenA,lenB);
                F[I].getBoxMeans();
            } break;
            case('c'):{ F[I].removeLastPoint(); } break;
            case ('n'):{ I++; I%=F.size();} break;
            case ('p'):{
                F[I].polyOrder ++;
                F[I].polyOrder %= 10;
            } break;
            case ('b'):{
                F[I].nBins ++;
                F[I].nBins %= 100;
            } break;
            case ('w'):{ F[I].printMeans(); } break;
            case ('q'):{return 0;} break;
        }
    }
}
