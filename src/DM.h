#pragma once

#include <vector>
using std::vector;
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using cv::Mat;
using cv::FONT_HERSHEY_SIMPLEX;
using cv::createTrackbar;
using cv::setTrackbarPos;
#include <morph/HdfData.h>
using morph::HdfData;
#include "FrameData.h"

// OpenCV functions mostly expect colours in Blue-Green-Red order
#define SF_BLUE Scalar(255,0,0,10)
#define SF_GREEN Scalar(0,255,0,10)
#define SF_RED Scalar(0,0,255,10)
#define SF_YELLOW Scalar(0,255,255,10)
#define SF_BLACK Scalar(0,0,0)
#define SF_WHITE Scalar(255,255,255)
#define SF_C1 Scalar(238,121,159) // mediumpurple2
#define SF_C2 Scalar(238,58,178) // darkorchid2

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
    //! The thickness of each brain slice in mm
    float thickness = 0.2f;

public:
    //! The instance public function. Short on purpose
    static DM* i (void) {
        if (DM::pInstance == 0) {
            DM::pInstance = new DM;
            DM::i()->init();
        }
        return DM::pInstance;
    }
    void init (void) {
        this->vFrameData.clear();
    }
    //! Add a frame to vFrameData
    void addFrame (Mat& frameImg, const string& frameImgFilename) {
        FrameData fd(frameImg);
        fd.filename = frameImgFilename;
        // Increment layer info. Best might be to use JSON info for layer positions as
        // they are unlikely always to be in perfect increments.
        if (this->vFrameData.empty()) {
            fd.layer_x = 0.0f;
            fd.idx = 0;
        } else {
            fd.layer_x = this->vFrameData.back().layer_x + DM::i()->thickness;
            fd.idx = this->vFrameData.back().idx + 1;
        }
        this->vFrameData.push_back (fd);
    }
    //! Return the size of vFrameData
    unsigned int getNumFrames (void) const {
        return this->vFrameData.size();
    }
    //! get current frame. Short name on purpose.
    FrameData* gcf (void) {
        if (!this->vFrameData.empty()) {
            return &(this->vFrameData[this->I]);
        }
        return (FrameData*)0;
    }
    //! Get the current frame number, counting from 1 like a human.
    int getFrameNum (void) const {
        return 1+this->I;
    }
    // Get a pointer to the persistent Mat img member attribute
    Mat* getImg (void) {
        return &(this->img);
    }
    //! Make the next frame current (or cycle back to the first)
    void nextFrame (void) {
        ++this->I %= this->vFrameData.size();
        DM::i()->gcf()->binA = this->binA;
        DM::i()->gcf()->binB = this->binB;
        DM::i()->gcf()->nBinsTarg = this->nBinsTarg;
    }
    //! Clone the current frame into Mat img
    void cloneFrame (void) {
        this->img = this->vFrameData[this->I].frame.clone();
    }
    //! Write frames to HdfData
    void writeFrames (void) {
        HdfData d(this->logname);
        for (auto f : this->vFrameData) {
            f.write (d);
        }
    }
    //! The application window name
    const string winName = "StaleFish";
    //! Saved/last cursor position
    int x = 0;
    int y = 0;
    //! Target number of bins; used by bins slider. Apply this to the framedata
    int nBinsTarg = 100;
    //! The bin lengths, set with a slider.
    int binA = 0;
    int binB = 40;
    //! Filename for writing
    string logname = "./stalefish.log";

    /*!
     * UI methods. Could probably un-static these.
     */
    static void onmouse (int event, int x, int y, int flags, void* param) {
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
    static void ontrackbar_boxes (int val, void*) {
        FrameData* cf = DM::i()->gcf();
        cf->binA = DM::i()->binA;
        cf->binB = DM::i()->binB;
        cf->setShowBoxes (true);
        cf->refreshBoxes (-cf->binA, cf->binB);
        DM::onmouse (CV_EVENT_MOUSEMOVE, -1, -1, 0, NULL);
    }

    //! On trackbar change, refresh the size of the boxes
    static void ontrackbar_nbins (int val, void*) {
        FrameData* cf = DM::i()->gcf();
        cf->nBinsTarg = DM::i()->nBinsTarg;
        if (cf->nBinsTarg < 2) {
            cf->nBinsTarg = 2;
        }
        cf->setBins (cf->nBinsTarg);
        cf->setShowBoxes (true);
        cf->updateFit();
        cf->refreshBoxes (-cf->binA, cf->binB);
        DM::onmouse (CV_EVENT_MOUSEMOVE, -1, -1, 0, NULL);
    }

    static void createTrackbars (void) {
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
};

//! Globally initialise DM instance pointer to NULL
DM* DM::pInstance = 0;
