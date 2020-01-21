#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include "tools.h"
#include "FrameData.h"
using namespace cv;
using namespace std;

//! Singleton pattern class to hold framedata
class DataManager
{
private:
    //! Private constructor/destructor
    DataManager();
    ~DataManager();

    //! A pointer returned to the single instance of this class
    static DataManager* pInstance;

    //! The frame data which is being managed
    vector<FrameData> frames;

    //! The index into frames which is the current image
    int I = 0;

    //! The current image
    Mat img;

public:
    static DataManager* instance (void) {
        if (DataManager::pInstance == 0) {
            DataManager::pInstance = new DataManager;
        }
        return DataManager::pInstance;
    }

    void addFrame (FrameData& fr) {
        this->frames.push_back (fr);
    }

    vector<FrameData>* getFrames (void) {
        return &(this->frames);
    }

    unsigned int getNumFrames (void) {
        return this->frames.size();
    }

    FrameData* getCurrentFrame (void) {
        if (!this->frames.empty()) {
            return &(this->frames[I]);
        }
        return (FrameData*)0;
    }

    int getI (void) {
        return I;
    }

    Mat* getImg (void) {
        return &(img);
    }

    void nextFrame (void) {
        ++I;
        if (I == this->frames.size()) {
            // Ooops, off the end, restep
            --I;
        }
    }

    void previousFrame (void) {
        --I;
        if (I < 0) {
            // Ooops, off the end, restep
            ++I;
        }
    }

    //! Clone the current frame into Mat img
    void cloneFrame (void) {
        this->img = this->frames[I].frame.clone();
    }
};

// Globally init datamanager instance pointer to null
DataManager* DataManager::pInstance = 0;

// *** DISPLAY ***
void onmouse(int event, int x, int y, int flags, void* param) {
    cout << "onmouse called" << endl;
    Point pt = Point(x,y);
    if (event==CV_EVENT_LBUTTONDOWN) {
        DataManager::instance()->getCurrentFrame()->P.push_back(pt);
    }
    DataManager::instance()->cloneFrame();

    circle(img,pt,5,Scalar(0,0,255),1);
    for(int i=1;i<DataManager::instance()->getCurrentFrame()->P.size();i++){
        line(img,DataManager::instance()->getCurrentFrame()->P[i-1],DataManager::instance()->getCurrentFrame()->P[i],Scalar(255,0,0),1);
    }
    for(int i=0;i<DataManager::instance()->getCurrentFrame()->P.size();i++){
        circle(img,DataManager::instance()->getCurrentFrame()->P[i],5,Scalar(255,0,0),1);
    }
    cout << "fors mid" << endl;
    if (DataManager::instance()->getCurrentFrame()->P.size()){
        line(img,DataManager::instance()->getCurrentFrame()->P[DataManager::instance()->getCurrentFrame()->P.size()-1],pt,Scalar(255,0,0),1);
    }
    cout << "fors mid2" << endl;
    for(int i=1;i<DataManager::instance()->getCurrentFrame()->fitted.size();i++){
        line(img,DataManager::instance()->getCurrentFrame()->fitted[i-1],DataManager::instance()->getCurrentFrame()->fitted[i],Scalar(0,255,0),1);
    }
    cout << "fors mid3" << endl;
    line(img,DataManager::instance()->getCurrentFrame()->axis[0],DataManager::instance()->getCurrentFrame()->axis[1],Scalar(0,0,255),1);
    cout << "fors mid4" << endl;
    for(int i=0;i<DataManager::instance()->getCurrentFrame()->origins.size();i++){
        line(img,DataManager::instance()->getCurrentFrame()->origins[i],DataManager::instance()->getCurrentFrame()->tangents[i],Scalar(0,255,255),1);
    }
    cout << "fors done" << endl;
    stringstream ss;
    ss<<"Frame: "<<DataManager::instance()->getI()<<"/"<<DataManager::instance()->getNumFrames()<<", "<<"Poly order: "<<DataManager::instance()->getCurrentFrame()->polyOrder<<", Bins: "<<DataManager::instance()->getCurrentFrame()->nBins;
    putText(img,ss.str(),Point(30,30), FONT_HERSHEY_SIMPLEX, 0.8, Scalar(0,0,0), 1, CV_AA);
    cout << "imshow()..." << endl;
    imshow("StaleFish", img);
}

// *** MAIN PROGRAM ***

int main(int argc, char** argv){

    cout << "Start" << endl;
    F.clear();
    double lenA = 10.;
    double lenB = 50.;
    I = 0;
    for (int i=1;i<argc;i++){
        char* imageName = argv[i];
        cout << "imread " << imageName << endl;
        Mat frame = imread(imageName, IMREAD_COLOR );
        if (frame.empty()) {
            cout <<  "Could not open or find the image" << endl;
            return -1;
        }
        cout << "push_back on F" << endl;
        F.push_back (FrameData (frame.clone()));
    }

    // Fails to draw the window
    namedWindow("StaleFish", 1);
    //img = F[0].frame.clone();
    //imshow("Stalefish", img);

    // Causes "no opengl in opencv" error
    //updateWindow("StaleFish");

    cout << "Set mouse callback..." << endl;
    setMouseCallback("StaleFish", onmouse, &img);
    cout << "Done setting mouse callback" << endl;

    // *** MAIN LOOP ***

    while (1) {
        onmouse(CV_EVENT_MOUSEMOVE,0,0,0,NULL);
        char k = waitKey(0);
        switch(k){
            case('x'):{
                DataManager::instance()->getCurrentFrame()->updateFit();
                DataManager::instance()->getCurrentFrame()->refreshBoxes(lenA,lenB);
                DataManager::instance()->getCurrentFrame()->getBoxMeans();
            } break;
            case('c'):{ DataManager::instance()->getCurrentFrame()->removeLastPoint(); } break;
            case ('n'):{ I++; I%=F.size();} break;
            case ('p'):{
                DataManager::instance()->getCurrentFrame()->polyOrder ++;
                DataManager::instance()->getCurrentFrame()->polyOrder %= 10;
            } break;
            case ('b'):{
                DataManager::instance()->getCurrentFrame()->nBins ++;
                DataManager::instance()->getCurrentFrame()->nBins %= 100;
            } break;
            case ('w'):{ DataManager::instance()->getCurrentFrame()->printMeans(); } break;
            case ('q'):{return 0;} break;
        }
    }
}
