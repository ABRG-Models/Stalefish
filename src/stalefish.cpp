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

// OpenCV functions mostly expect colours in Blue-Green-Red order
#define SF_BLUE Scalar(255,0,0,10)
#define SF_GREEN Scalar(0,255,0,10)
#define SF_RED Scalar(0,0,255,10)
#define SF_YELLOW Scalar(0,255,255,10)
#define SF_BLACK Scalar(0,0,0)
#define SF_WHITE Scalar(255,255,255)
#define SF_C1 Scalar(238,121,159) // mediumpurple2
#define SF_C2 Scalar(238,58,178) // darkorchid2

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

    // red circle under the cursor
    circle (*pImg, pt, 5, SF_RED, 2);

    // First the lines in the preceding PP point-sets:
    for (size_t j=0; j<cf->PP.size(); j++) {
        Scalar linecol = j%2 ? SF_RED : SF_BLUE;
        for (size_t i=0; i<cf->PP[j].size(); i++) {
            circle (*pImg, cf->PP[j][i], 5, linecol, 1);
            if (i) { line (*pImg, cf->PP[j][i-1], cf->PP[j][i], linecol, 2); }
        }
    }

#if 1
    // Add the control points in similar colours
    list<BezCurve<double>> theCurves = cf->bcp.curves;
    size_t j = 0;
    for (auto curv : theCurves) {
        Scalar linecol = j%2 ? SF_RED : SF_BLUE;
        vector<pair<double,double>> ctrls = curv.getControls();
        for (size_t cc = 0; cc<ctrls.size(); ++cc) {
            Point p1(ctrls[cc].first, ctrls[cc].second);
            circle (*pImg, p1, 5, linecol, -1);
#if 0
            if (cc==0 || cc==ctrls.size()-1) {
                circle (*pImg, p1, 2, SF_BLACK, -1);
            } else {
                circle (*pImg, p1, 2, SF_WHITE, -1);
            }
#endif
        }
        Point ps(ctrls[0].first, ctrls[0].second);
        Point pe(ctrls[1].first, ctrls[1].second);
        line (*pImg, ps, pe, SF_GREEN, 1);
        Point ps2(ctrls[ctrls.size()-2].first, ctrls[ctrls.size()-2].second);
        Point pe2(ctrls[ctrls.size()-1].first, ctrls[ctrls.size()-1].second);
        line (*pImg, ps2, pe2, SF_GREEN, 1);

        j++;
    }
#endif

    // Then the current point set:
    for (size_t i=0; i<cf->P.size(); i++) {
        circle (*pImg, cf->P[i], 5, SF_C1, 1);
        if (i) { line (*pImg, cf->P[i-1], cf->P[i], SF_C1, 2); }
    }

    // blue
    if (cf->P.size()) {
        line (*pImg, cf->P[cf->P.size()-1], pt, SF_C1, 1);
    }

    // green. This is the fit.
    if (cf->showFit == true) {
        for(size_t i=1; i<cf->fitted.size(); i++) {
            line (*pImg, cf->fitted[i-1], cf->fitted[i], SF_GREEN, 1);
        }

        // Axis line, relevant for polynomial fit
        line (*pImg, cf->axis[0], cf->axis[1], SF_RED, 1);

        // yellow. pointsInner to pointsOuter
        for (size_t i=0; i<cf->pointsInner.size(); i++) {
            line (*pImg, cf->pointsInner[i], cf->pointsOuter[i], SF_YELLOW, 1);
        }
    }

    stringstream ss;
    ss << "Frame: " << DM::i()->getFrameNum() << "/" << DM::i()->getNumFrames()
       << " " << DM::i()->gcf()->getFitInfo();
    putText (*pImg, ss.str(), Point(30,30), FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, CV_AA);

    imshow ("StaleFish", *pImg);
}

//! Main entry point
int main (int argc, char** argv)
{
    if (argc < 2) {
        cout << "Please supply at least one image filename" << endl;
        return 1;
    }

    double lenA = 0.;
    double lenB = 40.;
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
            DM::i()->gcf()->setShowFit(true);
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes(lenA,lenB);
            DM::i()->gcf()->getBoxMeans();
            break;
        }
        case('c'):
        {
            DM::i()->gcf()->removeLastPoint();
            //DM::i()->gcf()->updateFit();
            //DM::i()->gcf()->refreshBoxes(lenA,lenB);
            //DM::i()->gcf()->getBoxMeans();
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
        case ('s'):
        {
            // show/unshow the fit and boxes
            DM::i()->gcf()->toggleShowFit();
            break;
        }
        case ('b'):
        {
            DM::i()->gcf()->incBins();
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes(lenA,lenB);
            DM::i()->gcf()->getBoxMeans();
            break;
        }
        case ('B'):
        {
            DM::i()->gcf()->incBins(10);
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes(lenA,lenB);
            DM::i()->gcf()->getBoxMeans();
            break;
        }
        case ('m'):
        {
            // Change 'mode'
            DM::i()->gcf()->toggleCurveType();
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes(lenA,lenB);
            DM::i()->gcf()->getBoxMeans();
            break;
        }
        case ('w'):
        {
            DM::i()->gcf()->printMeans();
            break;
        }
        case (' '):
        {
            DM::i()->gcf()->nextCurve();
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
