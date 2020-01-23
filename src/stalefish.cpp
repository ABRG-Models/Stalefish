// The Data Manager class
#include "DM.h"
// The image-and-associated-fit-parameters class
#include "FrameData.h"
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using namespace cv;
#include <iostream>
using std::cout;
using std::endl;
#include <sstream>
using std::stringstream;
#include <unistd.h>

//! Global mouse callback. Could perhaps be a member of DM
void onmouse (int event, int x, int y, int flags, void* param)
{
    Point pt = Point(x,y);
    if (event == CV_EVENT_LBUTTONDOWN) {
        DM::i()->gcf()->P.push_back (pt);
    }
    DM::i()->cloneFrame();

    // Make copies of pointers to neaten up the code, below
    Mat* pImg = DM::i()->getImg();
    FrameData* cf = DM::i()->gcf();

    circle (*pImg, pt, 5, Scalar(0,0,255), 1);

    for (size_t i=1; i<cf->P.size(); i++) {
        line (*pImg, cf->P[i-1], cf->P[i], Scalar(255,0,0), 1);
    }
    for (size_t i=0; i<cf->P.size(); i++){
        circle (*pImg, cf->P[i], 5, Scalar(255,0,0), 1);
    }

    if (cf->P.size()) {
        line (*pImg, cf->P[cf->P.size()-1], pt, Scalar(255,0,0),1);
    }

    for(size_t i=1; i<cf->fitted.size(); i++) {
        line (*pImg, cf->fitted[i-1], cf->fitted[i], Scalar(0,255,0), 1);
    }

    line (*pImg, cf->axis[0], cf->axis[1], Scalar(0,0,255), 1);

    for (size_t i=0; i<cf->origins.size(); i++) {
        line (*pImg, cf->origins[i], cf->tangents[i], Scalar(0,255,255), 1);
    }

    stringstream ss;
    ss << "Frame: " << DM::i()->getFrameNum() << "/" << DM::i()->getNumFrames()
       << " " << DM::i()->gcf()->getFitInfo();
    putText (*pImg, ss.str(), Point(30,30), FONT_HERSHEY_SIMPLEX, 0.8, Scalar(0,0,0), 1, CV_AA);

    imshow ("StaleFish", *pImg);
}

//! Main entry point
int main (int argc, char** argv)
{
    if (argc < 2) {
        cout << "Please supply at least one image filename" << endl;
        return 1;
    }

    double lenA = 10.;
    double lenB = 50.;
    for (int i=1; i<argc; i++){
        char* imageName = argv[i];
        cout << "imread " << imageName << endl;
        Mat frame = imread(imageName, IMREAD_COLOR);
        if (frame.empty()) {
            cout <<  "Could not open or find the image" << endl;
            return -1;
        }
        DM::i()->addFrame (frame);
    }
    namedWindow ("StaleFish", WINDOW_AUTOSIZE);
    setMouseCallback ("StaleFish", onmouse, DM::i()->getImg());

    // *** MAIN LOOP ***
    while (1) {
        onmouse (CV_EVENT_MOUSEMOVE, 0, 0, 0, NULL);
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
        case ('m'):{
            // Change 'mode'
            DM::i()->gcf()->toggleCurveType();
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
