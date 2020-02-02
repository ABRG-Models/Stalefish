// The Data Manager/drawing class
#include "DM.h"
// The image-and-associated-fit-parameters class
#include "FrameData.h"
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using cv::namedWindow;
using cv::setMouseCallback;
using cv::waitKey;
using cv::WINDOW_AUTOSIZE;
using cv::imread;
using cv::IMREAD_COLOR;
#include <iostream>
using std::cout;
using std::endl;
#include <sstream>
using std::stringstream;
#include <unistd.h>

//! Main entry point
int main (int argc, char** argv)
{
    if (argc < 2) {
        cout << "Please supply at least one image filename" << endl;
        return 1;
    }

    // FIXME: Use morph::Config to read in the files to read, and pertient information
    // such as the slice thickness, anything about how to read the image data, where to
    // save the logs, etc.

    for (int i=1; i<argc; i++){
        char* imageName = argv[i];
        cout << "imread " << imageName << endl;
        Mat frame = imread (imageName, IMREAD_COLOR);
        if (frame.empty()) {
            cout <<  "Could not open or find the image" << endl;
            return -1;
        }
        DM::i()->addFrame (frame, imageName);
    }
    namedWindow (DM::i()->winName, WINDOW_AUTOSIZE);
    // Make sure there's an image in DM to start with
    DM::i()->cloneFrame();
    setMouseCallback (DM::i()->winName, DM::onmouse, DM::i()->getImg());
    DM::createTrackbars();
    // Init current frame with binA, binB and nBinsTarg:
    DM::i()->gcf()->binA = DM::i()->binA;
    DM::i()->gcf()->binB = DM::i()->binB;
    DM::i()->gcf()->nBinsTarg = DM::i()->nBinsTarg;

    // *** MAIN LOOP ***
    while (1) {
        DM::onmouse (CV_EVENT_MOUSEMOVE, -1, -1, 0, NULL);
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
            DM::i()->writeFrames();
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
