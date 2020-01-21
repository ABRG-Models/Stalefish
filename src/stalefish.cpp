#include <stdio.h>
#include <unistd.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include "tools.h"
#include "FrameData.h"
using namespace cv;
using namespace std;

//! Singleton pattern data manager class to hold framedata
class DM
{
private:
    //! Private constructor/destructor
    DM() {};
    ~DM() {};
    //! A pointer returned to the single instance of this class
    static DM* pInstance;
    //! The frame data which is being managed
    vector<FrameData> vFrameData;
    //! The index into frames which is the current image
    int I = 0;
    //! The current image
    Mat img;

public:
    //! The instance public function. Short on purpose
    static DM* i (void) {
        if (DM::pInstance == 0) {
            DM::pInstance = new DM;
        }
        return DM::pInstance;
    }
    //! Add a frame to vFrameData
    void addFrame (Mat& frameImg) {
        this->vFrameData.push_back (FrameData (frameImg.clone()));
    }
    //! Return the size of vFrameData
    unsigned int getNumFrames (void) {
        return this->vFrameData.size();
    }
    //! get current frame. Short name on purpose.
    FrameData* gcf (void) {
        if (!this->vFrameData.empty()) {
            return &(this->vFrameData[I]);
        }
        return (FrameData*)0;
    }
    //! Get the current frame number, counting from 1 like a human.
    int getFrameNum (void) {
        return 1+I;
    }
    // Get a pointer to the persistent Mat img member attribute
    Mat* getImg (void) {
        return &(img);
    }
    //! Make the next frame current (or cycle back to the first)
    void nextFrame (void) {
        ++I %= this->vFrameData.size();
    }
    //! Clone the current frame into Mat img
    void cloneFrame (void) {
        this->img = this->vFrameData[I].frame.clone();
    }
};

// Globally init datamanager instance pointer to null
DM* DM::pInstance = 0;

// *** DISPLAY ***
void onmouse (int event, int x, int y, int flags, void* param) {

    Point pt = Point(x,y);
    if (event == CV_EVENT_LBUTTONDOWN) {
        DM::i()->gcf()->P.push_back(pt);
    }
    DM::i()->cloneFrame();

    // Make copies of pointers to neaten up the code, below
    Mat* pImg = DM::i()->getImg();
    FrameData* cf = DM::i()->gcf();

    circle (*pImg, pt, 5, Scalar(0,0,255), 1);

    for (int i=1; i<cf->P.size(); i++) {
        line (*pImg, cf->P[i-1], cf->P[i], Scalar(255,0,0), 1);
    }
    for (int i=0; i<cf->P.size(); i++){
        circle (*pImg, cf->P[i], 5, Scalar(255,0,0), 1);
    }

    if (cf->P.size()) {
        line (*pImg, cf->P[cf->P.size()-1], pt, Scalar(255,0,0),1);
    }

    for(int i=1; i<cf->fitted.size(); i++) {
        line (*pImg, cf->fitted[i-1], cf->fitted[i], Scalar(0,255,0), 1);
    }

    line (*pImg, cf->axis[0], cf->axis[1], Scalar(0,0,255), 1);

    for (int i=0; i<cf->origins.size(); i++) {
        line (*pImg, cf->origins[i], cf->tangents[i], Scalar(0,255,255), 1);
    }

    stringstream ss;
    ss << "Frame: " << DM::i()->getFrameNum() << "/" << DM::i()->getNumFrames()
       << ", Poly order: " << cf->polyOrder
       << ", Bins: " << cf->nBins;
    putText (*pImg, ss.str(), Point(30,30), FONT_HERSHEY_SIMPLEX, 0.8, Scalar(0,0,0), 1, CV_AA);

    imshow ("StaleFish", *pImg);
}

// *** MAIN PROGRAM ***

int main(int argc, char** argv){

    if (argc < 2) {
        cout << "Please supply at least one image filename" << endl;
        return 1;
    }

    double lenA = 10.;
    double lenB = 50.;
    for (int i=1;i<argc;i++){
        char* imageName = argv[i];
        cout << "imread " << imageName << endl;
        Mat frame = imread(imageName, IMREAD_COLOR );
        if (frame.empty()) {
            cout <<  "Could not open or find the image" << endl;
            return -1;
        }
        DM::i()->addFrame (frame);
    }
    namedWindow("StaleFish", WINDOW_AUTOSIZE);
    setMouseCallback("StaleFish", onmouse, DM::i()->getImg());

    // *** MAIN LOOP ***
    while (1) {
        onmouse(CV_EVENT_MOUSEMOVE,0,0,0,NULL);
        char k = waitKey(0);
        switch(k) {
        case('x'):
        {
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes(lenA,lenB);
            DM::i()->gcf()->getBoxMeans();
            break;
        }
        case('c'):
        {
            DM::i()->gcf()->removeLastPoint();
            break;
        }
        case ('n'):
        {
            DM::i()->nextFrame();
            break;
        }
        case ('p'):
        {
            DM::i()->gcf()->polyOrder ++;
            DM::i()->gcf()->polyOrder %= 10;
            break;
        }
        case ('b'):{
            DM::i()->gcf()->nBins ++;
            DM::i()->gcf()->nBins %= 100;
            break;
        }
        case ('w'):
        {
            DM::i()->gcf()->printMeans();
            break;
        }
        case ('q'):
        {
            return 0;
            break;
        }
        default:
            break;
        }

        // Allow graphing system to catch up
        usleep (1000);
    }
}
