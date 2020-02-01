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

//! Global mouse/drawing callback. Could perhaps be a member of DM?
static void onmouse (int event, int x, int y, int flags, void* param)
{
    Point pt = Point(x,y);
    if (x==-1 && y==-1) {
        pt = Point(DM::i()->x, DM::i()->y);
    } else {
        DM::i()->x = x;
        DM::i()->y = y;
    }
    if (event == CV_EVENT_LBUTTONDOWN) {
        DM::i()->gcf()->P.push_back (pt);
        DM::i()->gcf()->setShowUsers(true);
    }
    DM::i()->cloneFrame();

    // Make copies of pointers to neaten up the code, below
    Mat* pImg = DM::i()->getImg();
    FrameData* cf = DM::i()->gcf();

    // red circle under the cursor
    circle (*pImg, pt, 5, SF_RED, 1);

    if (cf->flags.test(ShowUsers) == true) {
        // First the lines in the preceding PP point-sets:
        for (size_t j=0; j<cf->PP.size(); j++) {
            Scalar linecol = j%2 ? SF_RED : SF_BLUE;
            for (size_t i=0; i<cf->PP[j].size(); i++) {
                circle (*pImg, cf->PP[j][i], 5, linecol, -1);
                if (i) { line (*pImg, cf->PP[j][i-1], cf->PP[j][i], linecol, 2); }
            }
        }
    }

    if (cf->flags.test(ShowCtrls)) {
        // Add the control points in similar colours
        list<BezCurve<double>> theCurves = cf->bcp.curves;
        size_t j = 0;
        for (auto curv : theCurves) {
            Scalar linecol = j%2 ? SF_RED : SF_BLUE;
            vector<pair<double,double>> ctrls = curv.getControls();
            for (size_t cc = 0; cc<ctrls.size(); ++cc) {
                Point p1(ctrls[cc].first, ctrls[cc].second);
                circle (*pImg, p1, 5, linecol, -1);
            }
            Point ps(ctrls[0].first, ctrls[0].second);
            Point pe(ctrls[1].first, ctrls[1].second);
            line (*pImg, ps, pe, SF_GREEN, 1);
            Point ps2(ctrls[ctrls.size()-2].first, ctrls[ctrls.size()-2].second);
            Point pe2(ctrls[ctrls.size()-1].first, ctrls[ctrls.size()-1].second);
            line (*pImg, ps2, pe2, SF_GREEN, 1);

            j++;
        }
    }

    if (cf->flags.test(ShowUsers) == true) {
        // Then the current point set:
        if (cf->PP.empty() || (!cf->PP.empty() && cf->P.size() > 1)) {
            for (size_t i=0; i<cf->P.size(); i++) {
                circle (*pImg, cf->P[i], 5, SF_BLACK, -1);
                if (i) { line (*pImg, cf->P[i-1], cf->P[i], SF_BLACK, 2); }
            }
        }
#if 1
        // line to cursor position?
        if ((cf->PP.empty() && cf->P.size() > 0)
            || (!cf->PP.empty() && cf->P.size() > 1)) {
            line (*pImg, cf->P[cf->P.size()-1], pt, SF_BLACK, 1);
        }
#endif
    }

    // green. This is the fit.
    if (cf->flags.test(ShowFits) == true) {
        for (size_t i=1; i<cf->fitted.size(); i++) {
            line (*pImg, cf->fitted[i-1], cf->fitted[i], SF_GREEN, 1);
        }
        // Axis line, relevant for polynomial fit
        if (cf->ct == CurveType::Poly) {
            line (*pImg, cf->axis[0], cf->axis[1], SF_RED, 1);
        }
    }

    if (cf->flags.test(ShowBoxes) == true) {
        // yellow. pointsInner to pointsOuter
        for (size_t i=0; i<cf->pointsInner.size(); i++) {
            line (*pImg, cf->pointsInner[i], cf->pointsOuter[i], SF_YELLOW, 1);
        }
    }

    stringstream ss;
    ss << "Frame: " << DM::i()->getFrameNum() << "/" << DM::i()->getNumFrames()
       << " " << DM::i()->gcf()->getFitInfo();
    putText (*pImg, ss.str(), Point(30,30), FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, CV_AA);

    imshow (DM::i()->winName, *pImg);
}

//! On any trackbar changing, refresh the boxes
static void ontrackbar_boxes (int val, void*)
{
    FrameData* cf = DM::i()->gcf();
    cf->binA = DM::i()->binA;
    cf->binB = DM::i()->binB;
    cf->setShowBoxes (true);
    cf->refreshBoxes (-cf->binA, cf->binB);
    onmouse (CV_EVENT_MOUSEMOVE, -1, -1, 0, NULL);
}

static void ontrackbar_nbins (int val, void*)
{
    FrameData* cf = DM::i()->gcf();
    cf->nBinsTarg = DM::i()->nBinsTarg;
    if (cf->nBinsTarg < 2) {
        cf->nBinsTarg = 2;
    }
    cf->setBins (cf->nBinsTarg);
    cf->setShowBoxes (true);
    cf->updateFit();
    cf->refreshBoxes (-cf->binA, cf->binB);
    onmouse (CV_EVENT_MOUSEMOVE, -1, -1, 0, NULL);
}

static void createTrackbars (void)
{
    // Set up trackbars. Have to do this for each frame
    string tbBinA = "Box A";
    string tbBinB = "Box B";
    string tbNBins = "Num bins";
    createTrackbar (tbBinA, DM::i()->winName, &DM::i()->binA, 200, ontrackbar_boxes);
    setTrackbarPos (tbBinA, DM::i()->winName, DM::i()->binA);
    createTrackbar (tbBinB, DM::i()->winName, &DM::i()->binB, 200, ontrackbar_boxes);
    setTrackbarPos (tbBinB, DM::i()->winName, DM::i()->binB);
    createTrackbar (tbNBins, DM::i()->winName, &DM::i()->nBinsTarg, 200, ontrackbar_nbins);
    setTrackbarPos (tbNBins, DM::i()->winName, DM::i()->nBinsTarg);
}

//! Main entry point
int main (int argc, char** argv)
{
    if (argc < 2) {
        cout << "Please supply at least one image filename" << endl;
        return 1;
    }

    for (int i=1; i<argc; i++){
        char* imageName = argv[i];
        cout << "imread " << imageName << endl;
        Mat frame = imread (imageName, IMREAD_COLOR);
        if (frame.empty()) {
            cout <<  "Could not open or find the image" << endl;
            return -1;
        }
        DM::i()->addFrame (frame);
    }
    namedWindow (DM::i()->winName, WINDOW_AUTOSIZE);
    // Make sure there's an image in DM to start with
    DM::i()->cloneFrame();
    setMouseCallback (DM::i()->winName, onmouse, DM::i()->getImg());
    createTrackbars();
    // Init current frame with binA, binB and nBinsTarg:
    DM::i()->gcf()->binA = DM::i()->binA;
    DM::i()->gcf()->binB = DM::i()->binB;
    DM::i()->gcf()->nBinsTarg = DM::i()->nBinsTarg;

    // *** MAIN LOOP ***
    while (1) {
        onmouse (CV_EVENT_MOUSEMOVE, -1, -1, 0, NULL);
        char k = waitKey(0);
        switch(k) {
        // 1 to 4 - select what is shown
        case ('1'):
        {
            DM::i()->gcf()->toggleShowCtrls();
            break;
        }
        case ('2'):
        {
            DM::i()->gcf()->toggleShowUsers();
            break;
        }
        case ('3'):
        {
            DM::i()->gcf()->toggleShowFits();
            break;
        }
        case ('4'):
        {
            DM::i()->gcf()->toggleShowBoxes();
            break;
        }
        // Perform a fit
        case('x'):
        {
            DM::i()->gcf()->setShowFits (true);
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes (-DM::i()->gcf()->binA, DM::i()->gcf()->binB);
            break;
        }
        // Remove last point
        case('c'):
        {
            DM::i()->gcf()->removeLastPoint();
            break;
        }
        // Go to the next frame
        case ('n'):
        {
            DM::i()->nextFrame();
            break;
        }
        // Polyfit mode only - increment order
        case ('p'):
        {
            if (DM::i()->gcf()->ct == CurveType::Bezier) {
                break;
            }
            DM::i()->gcf()->polyOrder++;
            DM::i()->gcf()->polyOrder %= 10;
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes (-DM::i()->gcf()->binA, DM::i()->gcf()->binB);
            break;
        }
        // Increment number of bins
        case ('b'):
        {
            DM::i()->gcf()->incBins();
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes (-DM::i()->gcf()->binA, DM::i()->gcf()->binB);
            break;
        }
        // Increment number of bins faster
        case ('B'):
        {
            DM::i()->gcf()->incBins(10);
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes (-DM::i()->gcf()->binA, DM::i()->gcf()->binB);
            break;
        }
        case ('m'):
        {
            // Change fitting mode
            DM::i()->gcf()->toggleCurveType();
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes (-DM::i()->gcf()->binA, DM::i()->gcf()->binB);
            break;
        }
        case ('w'):
        {
            DM::i()->gcf()->getBoxMeans();
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
