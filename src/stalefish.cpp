// The Data Manager/drawing class
#include "DM.h"
// The image-and-associated-fit-parameters class
#include "FrameData.h"
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using cv::waitKey;
#include <iostream>
using std::cout;
using std::endl;
#include <sstream>
using std::stringstream;
#include <unistd.h>
#include <morph/Config.h>
using morph::Config;

//! Main entry point
int main (int argc, char** argv)
{
    if (argc < 2) {
        cout << "Please supply path to the json conf file" << endl;
        return 1;
    }
    string paramsfile (argv[1]);

    DM::i()->setup (paramsfile);

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
