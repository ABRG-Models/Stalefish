/*
 * Stalefish - main program
 * Copyright (C) 2022 Sebastian S. James, Stuart P. Wilson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

// The Data Manager/drawing class
#include "DM.h"
// The image-and-associated-fit-parameters class
#include "FrameData.h"
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

//! Main entry point
int main (int argc, char** argv)
{
    std::string paramsfile = "";
    if (argc < 2) {
        // If Mac, then show a window with this message saying 'using this as a cmd line app'
        std::cout << "Please supply path to the json conf file. Showing help/example screen." << std::endl;
        DM::i()->appmode = AppMode::NoFile;
        DM::i()->argv0 = std::string(argv[0]);
    } else {
        std::string pfile(argv[1]);
        paramsfile = pfile;
    }

    // In case we need to do something special when upgrading data formats
    if (argc > 2) {
        if (std::string(argv[2]) == std::string("old")) {
            DM::i()->readOldFormat = true;
        }
    }

    DM::i()->setup (paramsfile);

    while (1) {
#if 0
        // For some reason, cv::WND_PROP_VISIBLE is *always* -1, so this code doesn't work.
        std::cout << "cv::WND_PROP_VISIBLE for \"" << DM::i()->winName << "\" = "
                  << cv::getWindowProperty(DM::i()->winName, cv::WND_PROP_VISIBLE) << std::endl;
        if (DM::i()->firstCall == false && cv::getWindowProperty(DM::i()->winName, cv::WND_PROP_VISIBLE) == -1.0) {
            // Time to exit
            std::cout << "Would be time to exit now..." << std::endl;
        } // else carry on
# define WAIT_TIME 30 // ms
#else
# define WAIT_TIME 0 // wait infinitely
#endif
        DM::onmouse (cv::EVENT_MOUSEMOVE, -1, -1, 0, NULL);
        char k = cv::waitKey (WAIT_TIME);
        switch(k) {
        // 1 to 4 - select what is shown
        case ('1'):
        {
            DM::i()->toggleShowCtrls();
            DM::i()->gcf()->setShowCtrls (DM::i()->flags.test(AppShowCtrls));
            break;
        }
        case ('2'):
        {
            DM::i()->toggleShowUsers();
            DM::i()->gcf()->setShowUsers (DM::i()->flags.test(AppShowUsers));
            break;
        }
        case ('3'):
        {
            DM::i()->toggleShowFits();
            DM::i()->gcf()->setShowFits (DM::i()->flags.test(AppShowFits));
            break;
        }
        case ('4'):
        {
            DM::i()->toggleShowBoxes();
            DM::i()->gcf()->setShowBoxes (DM::i()->flags.test(AppShowBoxes));
            break;
        }
        case ('5'):
        {
            DM::i()->toggleShowText();
            break;
        }
        case ('8'):
        {
            DM::i()->moveFrameBack();
            break;
        }
        case ('9'):
        {
            DM::i()->moveFrameNext();
            break;
        }
        // Apply bin values across all frames
        case ('B'):
        {
            DM::i()->updateAllBins();
            break;
        }
        // Clear all Curves
        case ('C'):
        {
            if (DM::i()->clearAllPending == true) {
                DM::i()->clearAllCurves();
                DM::i()->clearAllPending = false;
            } else {
                DM::i()->clearAllPending = true;
            }
            break;
        }
        // Cancel a pending 'clear all curves' or 'export file'. char(27) is the Esc key.
        case (char(27)):
        {
            DM::i()->clearAllPending = false;
            DM::i()->exportPending = false;
            DM::i()->importPending = false;
            DM::i()->exitPending = false;
            break;
        }
        // Perform all fits
        case ('F'):
        {
            DM::i()->updateAllFits();
            break;
        }
        // Perform a fit
        case ('f'):
        {
            DM::i()->setShowFits (true);
            DM::i()->gcf()->setShowFits (true);
            DM::i()->setShowBoxes (true);
            DM::i()->gcf()->setShowBoxes (true);
            DM::i()->gcf()->updateFit();
            DM::i()->gcf()->refreshBoxes (-DM::i()->gcf()->binA, DM::i()->gcf()->binB);
            break;
        }
        // Remove last point or freehand region (or first point if in 'work at the start' mode)
        case ('c'):
        {
            DM::i()->removeLastThing();
            break;
        }
        // Go to the next frame
        case ('n'):
        {
            DM::i()->nextFrame();
            break;
        }
        // Go back to the previous frame
        case ('b'):
        {
            DM::i()->previousFrame();
            break;
        }
        case ('l'):
        {
            if (DM::i()->importPending == true) {
                DM::i()->importLandmarks();
                DM::i()->importPending = false;
            } else {
                DM::i()->importPending = true;
            }
            break;
        }
        case ('i'):
        {
            if (DM::i()->importPending == true) {
                DM::i()->importCurves();
                DM::i()->importPending = false;
            } else {
                DM::i()->importPending = true;
            }
            break;
        }
        case ('j'):
        {
            if (DM::i()->importPending == true) {
                DM::i()->importFreehand();
                DM::i()->importPending = false;
            } else {
                DM::i()->importPending = true;
            }
            break;
        }
        case ('k'):
        {
            if (DM::i()->exportPending == true) {
                DM::i()->exportUserpoints();
                DM::i()->exportPending = false;
            } else {
                DM::i()->exportPending = true;
            }
            break;
        }
        case ('['):
        {
            if (DM::i()->exportPending == true) {
                DM::i()->exportInputModePointsCurrentFrame();
                DM::i()->exportPending = false;
            } else {
                DM::i()->exportPending = true;
            }
            break;
        }
        case ('p'):
        {
            if (DM::i()->exportPending == true) {
                // Depending on InputMode, export landmarks, curves or freehand
                DM::i()->exportInputModePoints();
                DM::i()->exportPending = false;
            } else {
                DM::i()->exportPending = true;
            }
            break;
        }
        case ('q'):
        {
            // Toggle hiding of the image
            DM::i()->hideImage = DM::i()->hideImage ? false : true;
            break;
        }
        case ('m'):
        {
            // Mirror
            DM::i()->gcf()->mirror();
            DM::onmouse (cv::EVENT_MOUSEMOVE, -1, -1, 0, NULL);
            break;
        }
        case ('o'):
        {
            // Change input mode
            DM::i()->cycleInputMode(); // cycle globally
            break;
        }
        case ('P'):
        {
            // Duplicate objects from previous frame (mode dependent)
            DM::i()->duplicateFromPrev();
            break;
        }
        case ('N'):
        {
            // Duplicate objects from next frame (mode dependent)
            DM::i()->duplicateFromNext();
            break;
        }
        case ('s'):
        {
            // Toggle to 's'tart of curve - add/remove points there, instead
            DM::i()->toggleStartEnd();
            break;
        }
        case ('h'):
        {
            DM::i()->toggleShowHelp();
            if (!DM::i()->flags.test(AppShowText)) {
                for (auto ht : DM::i()->helptxt) { std::cout << ht << std::endl; }
            }
            DM::onmouse (cv::EVENT_MOUSEMOVE, -1, -1, 0, NULL);
            break;
        }
        case ('w'):
        {
            DM::i()->writeFrames();
            break;
        }
        case ('W'):
        {
            // Write out 2D map to separate file (only)
            DM::i()->writeMap();
            break;
        }
        case (' '):
        {
            DM::i()->gcf()->nextCurve();
            break;
        }
        case ('r'):
        {
            DM::i()->toggleBlurWindow();
            break;
        }
        case ('E'):
        {
            DM::i()->toggleOffsWindow();
            break;
        }
        case ('e'):
        {
            DM::i()->showExampleProject();
            break;
        }
        case ('.'):
        {
            // Output current cursor position
            DM::i()->coutCursorpos();
            break;
        }
        case ('x'):
        {
            // Warn user, to avoid inadvertent data loss.
            if (DM::i()->exitPending == true) {
                DM::i()->exitPending = false;
                return 0;
            } else {
                DM::i()->exitPending = true;
            }
            break;
        }
        default:
            break;
        }

        // Allow graphing system to catch up
        usleep (1000);
    }
}
