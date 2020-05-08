#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <morph/HdfData.h>
#include <morph/Config.h>
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

//! OpenCV sliders can't be negative, so we have an offset.
#define BIN_A_OFFSET 200

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
    std::vector<FrameData> vFrameData;
    //! The index into vFrameData which is the current image
    int I = 0;
    //! The current image
    cv::Mat img;
    //! The thickness of each brain slice in mm
    float thickness = 0.05f;
    //! The application configuration
    morph::Config conf;
    //! Should the help text be shown?
    bool showHelp = false;
    //! Colour space parameters
    //@{
    std::string colourmodel = "monochrome";
    std::array<float, 9> colour_rot; // Colour space rotation
    std::array<float, 3> colour_trans; // Colour space pre-translation
    std::array<float, 2> ellip_axes; // red-green ellipse for "elliptical tube of expressing colours
    float luminosity_factor; // The slope of the linear luminosity vs signal fit.
    float luminosity_cutoff; // at what luminosity does the signal cut off to zero?
    //@}

    // Called by next/previousFrame. Take binA, binB from the frame and change the
    // sliders. Update the fit and refresh boxes.
    void refreshFrame (void) {
        this->binA = this->gcf()->binA+BIN_A_OFFSET;
        this->binB = this->gcf()->binB;
        this->nBinsTarg = this->gcf()->nBinsTarg;
        DM::updateTrackbars();
        this->gcf()->updateFit();
        this->gcf()->refreshBoxes (-(this->binA-BIN_A_OFFSET), this->binB);
    }

public:

    //! The instance public function. Short on purpose
    static DM* i (void) {
        if (DM::pInstance == 0) {
            DM::pInstance = new DM;
            DM::i()->init();
        }
        return DM::pInstance;
    }

    //! Initialize by clearing out vFrameData.
    void init (void) {
        this->vFrameData.clear();
    }

    /*!
     * Add frame @frameImg to vFrameData, setting the metadata attributes
     * @frameImgFilename (The filename for the image), @slice_x (position in the x
     * dimension) and @ppm (pixels per mm; the scale).
     */
    void addFrame (cv::Mat& frameImg, const std::string& frameImgFilename, const float& slice_x) {
        std::cout << "********** DM::addFrame ***********" << std::endl;
        FrameData fd(frameImg);
        fd.filename = frameImgFilename;
        fd.setParentStack (&this->vFrameData);
        // Increment layer index. Best might be to use JSON info for layer positions as
        // they are unlikely always to be in perfect increments.
        if (this->vFrameData.empty()) {
            fd.idx = 0;
            std::cout << "First frame; not setting previous" << std::endl;
        } else {
            fd.idx = this->vFrameData.back().idx + 1;
            // Set pointer to previous so slices can be aligned during FrameData::write or updateFit
            fd.setPrevious (this->vFrameData.back().idx);
            std::cout << "Subsequent frame; setting previous to " << (&this->vFrameData[this->vFrameData.back().idx]) << std::endl;
        }
        fd.layer_x = slice_x;
        fd.pixels_per_mm = (double)this->pixels_per_mm;
        fd.thickness = this->thickness;
        if (this->colourmodel == "allen") {
            fd.cmodel = ColourModel::AllenDevMouse;
        } else {
            fd.cmodel = ColourModel::Greyscale;
        }
        fd.colour_rot = this->colour_rot;
        fd.colour_trans = this->colour_trans;
        fd.ellip_axes = this->ellip_axes;
        fd.luminosity_factor = this->luminosity_factor;
        fd.luminosity_cutoff = this->luminosity_cutoff;

        std::cout << "Before read, binA=" << fd.binA << std::endl;
        std::cout << "             binB=" << fd.binB << std::endl;

        // Read, opportunistically
        try {
            morph::HdfData d(this->datafile, true); // true for read
            fd.read (d, this->readOldFormat);
            if (fd.flags.test (Mirrored)) {
                fd.mirror_image_only();
            }
            std::cout << "DM::addFrame: Calling FrameData::updateFit()" << std::endl;
            fd.updateFit();
        } catch (...) {
            // No problem, just carry on
        }
        std::cout << "After read,  binA=" << fd.binA << std::endl;
        std::cout << "             binB=" << fd.binB << std::endl;

        this->vFrameData.push_back (fd);
    }

    //! Return the size of vFrameData
    unsigned int getNumFrames (void) const {
        return this->vFrameData.size();
    }

    //! Copy the current frame's bin parameters (binA, binB, nBinsTarg) to all the other frames.
    void updateAllBins (void) {
        int nfr = DM::i()->getNumFrames();
        int idx = DM::i()->gcf()->idx;
        for (int f = 0; f < nfr; ++f) {
            if (idx == f) {
                continue;
            }
            this->vFrameData[f].binA = this->vFrameData[idx].binA;
            this->vFrameData[f].binB = this->vFrameData[idx].binB;
            this->vFrameData[f].nBinsTarg = this->vFrameData[idx].nBinsTarg;
            this->vFrameData[f].updateFit();
            this->vFrameData[f].refreshBoxes (-this->vFrameData[f].binA, this->vFrameData[f].binB);
        }
    }

    //! Update all fits - i.e. for every frame in the stack
    void updateAllFits (void) {
        int nfr = DM::i()->getNumFrames();
        for (int f = 0; f < nfr; ++f) {
            this->vFrameData[f].setShowFits (true);
            this->vFrameData[f].setShowBoxes (true);
            this->vFrameData[f].updateFit();
            this->vFrameData[f].refreshBoxes (-this->vFrameData[f].binA, this->vFrameData[f].binB);
        }
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
    cv::Mat* getImg (void) {
        return &(this->img);
    }

    //! Make the next frame current (or cycle back to the first)
    void nextFrame (void) {
        ++this->I %= this->vFrameData.size();
        this->refreshFrame();
    }

    //! Back up a frame
    void previousFrame (void) {
        this->I = --this->I < 0 ? this->vFrameData.size()-1 : this->I;
        this->refreshFrame();
    }

    //! Clone the current frame into Mat img
    void cloneFrame (void) {
        this->img = this->vFrameData[this->I].frame.clone();
    }

    //! Write frames to HdfData
    void writeFrames (void) {
        morph::HdfData d(this->datafile);
        for (auto f : this->vFrameData) {
            f.getBoxMeans();
            f.write (d);
            // Also build up an "overall" data store of the bins
        }
        int nf = this->vFrameData.size();
        d.add_val("/nframes", nf);
        std::cout << "writeFrames complete: All frames written to HDF5" << std::endl;
    }

    //! Toogle showHelp
    void toggleHelp (void) {
        this->showHelp = !this->showHelp;
    }

    //! The application window name
    const std::string winName = "StaleFish";
    //! Saved/last cursor position
    int x = 0;
    int y = 0;
    //! Target number of bins; used by bins slider. Apply this to the framedata
    int nBinsTarg = 100;
    //! The bin lengths, set with a slider.
    int binA = 0+BIN_A_OFFSET; // 200 means the slider is in the middle
    int binB = 40;
    //! Filename for writing
    std::string datafile = "unset.h5";
    //! How many pixels in the image is 1mm?
    float pixels_per_mm = 100.0f;

    //! Set to true to read in an old data format project, which is then written out in new format.
    bool readOldFormat = false;

    //! Application setup
    void setup (const std::string& paramsfile) {

        // Set the HDF5 data file path based on the .json file path
        std::string::size_type jsonpos = paramsfile.find(".json");
        if (jsonpos == std::string::npos) {
            this->datafile = paramsfile + ".h5";
        } else {
            this->datafile = paramsfile.substr (0,jsonpos) + ".h5";
        }

        this->conf.init (paramsfile);
        if (!this->conf.ready) {
            std::cerr << "Error setting up JSON config: " << this->conf.emsg << ", exiting." << std::endl;
            exit (1);
        }

        // Set the scale from JSON, too
        this->pixels_per_mm = conf.getFloat ("pixels_per_mm", 100.0f);
        this->thickness = conf.getFloat ("thickness", 0.05f);

        // The colour space information, if relevant (Allen ISH images)
        // colourmodel - if exists, a string
        this->colourmodel = conf.getString ("colourmodel", "monochrome");
        // colour_rot - array<float, 9>
        const Json::Value cr = conf.getArray ("colour_rot");
        for (unsigned int i = 0; i < cr.size(); ++i) {
            this->colour_rot[i] = cr[i].asFloat();
        }
        // colour_trans - array<float, 3>
        const Json::Value ct = conf.getArray ("colour_trans");
        for (unsigned int i = 0; i < ct.size(); ++i) {
            this->colour_trans[i] = ct[i].asFloat();
        }
        // ellip_axes array<float, 2>
        const Json::Value ea = conf.getArray ("ellip_axes");
        this->ellip_axes[0] = ea[0].asFloat();
        this->ellip_axes[1] = ea[1].asFloat();
        // luminosity linear fit parameters
        this->luminosity_cutoff = conf.getFloat ("luminosity_cutoff", 255.0);
        this->luminosity_factor = conf.getFloat ("luminosity_factor", -0.00392); // -1/255

        // Loop over slices, creating a FrameData object for each.
        const Json::Value slices = conf.getArray ("slices");
        for (unsigned int i = 0; i < slices.size(); ++i) {
            Json::Value slice = slices[i];
            std::string fn = slice.get ("filename", "unknown").asString();
            float slice_x = slice.get ("x", 0.0).asFloat();

            std::cout << "imread " << fn << std::endl;
            cv::Mat frame = cv::imread (fn.c_str(), cv::IMREAD_COLOR);
            if (frame.empty()) {
                std::cout <<  "Could not open or find the image '" << fn << "', exiting." << std::endl;
                exit (1);
            }

            // scaling routine //
            float scaleFactor = conf.getFloat("scaleFactor", 1.0f); // pull scale factor from config json

            if (scaleFactor != 1.0f) {
                std::cout << "rescaling frame to scaleFactor: " << scaleFactor << std::endl;

                cv::Size scaledSize = cv::Size(round(frame.cols * scaleFactor), round(frame.rows * scaleFactor));
                cv::Mat scaledFrame = cv::Mat(scaledSize, frame.type());
                cv::resize (frame, scaledFrame, scaledSize, scaleFactor, scaleFactor, cv::INTER_LINEAR);

                frame.release(); // free original frame since we have resized it

                this->addFrame(scaledFrame, fn, slice_x);
            } else {
                // if we are at the default scale factor do not do anything
                this->addFrame(frame, fn, slice_x);
            }
        }

        //cv::namedWindow (this->winName, cv::WINDOW_AUTOSIZE);
        cv::namedWindow (this->winName, cv::WINDOW_NORMAL|cv::WINDOW_FREERATIO);
        // Make sure there's an image in DM to start with
        this->cloneFrame();
        cv::setMouseCallback (this->winName, DM::onmouse, this->getImg());
        // Init current frame with binA, binB and nBinsTarg taken from the current
        // frame... But... is that information stored? Yes, it is.
        this->binA = this->gcf()->binA+BIN_A_OFFSET;
        this->binB = this->gcf()->binB;
        this->nBinsTarg = this->gcf()->nBinsTarg;
        DM::createTrackbars();
        //this->gcf()->updateFit();
        this->gcf()->refreshBoxes (-(this->binA-BIN_A_OFFSET), this->binB);
    }

    //! In Bezier or Polynomial modes, draw curves and users points
    void draw_curves (const cv::Point& pt) {

        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        FrameData* cf = _this->gcf();

        // red circle under the cursor
        circle (*pImg, pt, 5, SF_RED, 1);

        if (cf->flags.test(ShowUsers) == true) {
            // First the lines in the preceding PP point-sets:
            for (size_t j=0; j<cf->PP.size(); j++) {
                Scalar linecol = j%2 ? SF_RED : SF_BLUE;
                for (size_t ii=0; ii<cf->PP[j].size(); ii++) {
                    circle (*pImg, cf->PP[j][ii], 5, linecol, -1);
                    if (ii) { line (*pImg, cf->PP[j][ii-1], cf->PP[j][ii], linecol, 2); }
                }
            }
        }

        if (cf->flags.test(ShowCtrls)) {
            // Add the control points in similar colours
            std::list<morph::BezCurve<double>> theCurves = cf->bcp.curves;
            size_t j = 0;
            for (auto curv : theCurves) {
                cv::Scalar linecol = j%2 ? SF_RED : SF_BLUE;
                std::vector<std::pair<double,double>> ctrls = curv.getControls();
                for (size_t cc = 0; cc<ctrls.size(); ++cc) {
                    cv::Point p1(ctrls[cc].first, ctrls[cc].second);
                    cv::circle (*pImg, p1, 5, linecol, -1);
                }
                cv::Point ps(ctrls[0].first, ctrls[0].second);
                cv::Point pe(ctrls[1].first, ctrls[1].second);
                line (*pImg, ps, pe, SF_GREEN, 1);
                cv::Point ps2(ctrls[ctrls.size()-2].first, ctrls[ctrls.size()-2].second);
                cv::Point pe2(ctrls[ctrls.size()-1].first, ctrls[ctrls.size()-1].second);
                line (*pImg, ps2, pe2, SF_GREEN, 1);

                j++;
            }
        }

        if (cf->flags.test(ShowUsers) == true) {
            // Then draw the current point set:
            if (cf->PP.empty() || (!cf->PP.empty() && cf->P.size() > 1)) {
                for (size_t ii=0; ii<cf->P.size(); ii++) {
                    circle (*pImg, cf->P[ii], 5, SF_GREEN, -1);
                    if (ii) { line (*pImg, cf->P[ii-1], cf->P[ii], SF_GREEN, 1); }
                }
            }
            // also draw a thin line to the cursor position
            if ((cf->PP.empty() && cf->P.size() > 0)
                || (!cf->PP.empty() && cf->P.size() > 1)) {
                line (*pImg, cf->P[cf->P.size()-1], pt, SF_GREEN, 1);
            }
        }

        // This is the fit line
        if (cf->flags.test(ShowFits) == true) {
            for (size_t ii=1; ii<cf->fitted.size(); ii++) {
                line (*pImg, cf->fitted[ii-1], cf->fitted[ii], SF_GREEN, 2);
            }
            // Axis line, relevant for polynomial fits only
            if (cf->ct == CurveType::Poly) {
                line (*pImg, cf->axis[0], cf->axis[1], SF_RED, 1);
            }
        }

        if (cf->flags.test(ShowBoxes) == true) {
            // The bins; pointsInner to pointsOuter
            for (size_t ii=0; ii<cf->pointsInner.size(); ii++) {
                line (*pImg, cf->pointsInner[ii], cf->pointsOuter[ii], SF_YELLOW, 1);
            }
        }
    }

    //! Draw freehand loops when in CurveType::Freehand mode
    void draw_freehand (const cv::Point& pt) {

        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        FrameData* cf = _this->gcf();

        for (size_t j=0; j<cf->FPE.size(); j++) {
            for (size_t ii=0; ii<cf->FPE[j].size(); ii++) {
                rectangle (*pImg, cf->FPE[j][ii], cf->FPE[j][ii], SF_RED, 1);
            }
        }

        // Then draw the current point set:
        for (size_t ii=0; ii<cf->FP.size(); ii++) {
            rectangle (*pImg, cf->FP[ii], cf->FP[ii], SF_GREEN, 1);
        }
    }

    /*!
     * UI methods. Could possibly un-static these.
     */
    static void onmouse (int event, int x, int y, int flags, void* param) {

        // Make copies of pointers to neaten up the code, below
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        FrameData* cf = _this->gcf();

        // What's the cv::Point under the mouse pointer?
        cv::Point pt = cv::Point(x,y);
        if (x==-1 && y==-1) {
            pt = cv::Point(_this->x, _this->y);
        } else {
            _this->x = x;
            _this->y = y;
        }

        // If the button is down, add it
        if (event == cv::EVENT_LBUTTONDOWN) {
            if (cf->ct == CurveType::Bezier || cf->ct == CurveType::Poly) {
                cf->P.push_back (pt);
                cf->setShowUsers(true);
            } else if (cf->ct == CurveType::Freehand) {
                //_this->lbutton_down = true;
                cf->addToFP (pt);
            }
        } else if (event == cv::EVENT_MOUSEMOVE
                   && (flags & cv::EVENT_FLAG_LBUTTON) == cv::EVENT_FLAG_LBUTTON) {
            // Now button is down, want to add any pixel that the mouse moves over
            cf->addToFP (pt);
        }

        _this->cloneFrame();

        // Code for drawing stuff when we're in a curve-fitting mode
        if (cf->ct == CurveType::Bezier || cf->ct == CurveType::Poly) {
            _this->draw_curves (pt);
        } else if (cf->ct == CurveType::Freehand) {
            _this->draw_freehand (pt);
        } else {
            std::cerr << "WARNING: unknown curve type in current frame (this should not occur)" << std::endl;
        }

        std::stringstream ss;
        int xh = 30;
        ss << "Frame: " << _this->getFrameNum() << "/" << _this->getNumFrames()
           << " " << cf->getFitInfo() << ". 'h' to toggle help.";
        putText (*pImg, ss.str(), cv::Point(xh,30), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);

        int yh = 90;
        int yinc = 40;
        if (_this->showHelp) {
            putText (*pImg, std::string("Use the sliders to control the bin parameters"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("1:   Toggle Bezier controls"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("2:   Toggle user points"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("3:   Toggle the fit line"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("4:   Toggle the bins"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("Spc: Next curve"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("c:   Cancel last point"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("f:   Update the fit"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("F:   Update ALL fits"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("B:   Copy current bin params to all frames"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            std::stringstream hh;
            hh << "w:   Save to file: " << _this->datafile;
            putText (*pImg, hh.str(),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("o:   Fit mode (Bezier, polynomial or freehand)"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("p:   In polynomial mode, change order"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("n:   Next frame"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("b:   Back to previous frame"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("m:   Mirror this frame"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("x:   Exit the program"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
        }

        imshow (_this->winName, *pImg);
    }

    //! On any trackbar changing, refresh the boxes
    static void ontrackbar_boxes (int val, void*) {
        DM* _this = DM::i();
        FrameData* cf = _this->gcf();
        cf->binA = _this->binA-BIN_A_OFFSET;
        cf->binB = _this->binB;
        cf->setShowBoxes (true);
        cf->refreshBoxes (-cf->binA, cf->binB);
        DM::onmouse (cv::EVENT_MOUSEMOVE, -1, -1, 0, NULL);
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
        DM::onmouse (cv::EVENT_MOUSEMOVE, -1, -1, 0, NULL);
    }

    static void createTrackbars (void) {
        // Set up trackbars. Have to do this for each frame
        std::string tbBinA = "Box A";
        std::string tbBinB = "Box B";
        std::string tbNBins = "Num bins";
        DM* _this = DM::i();
        std::cout << "createTrackbars: _this->binA=" << _this->binA << std::endl;
        cv::createTrackbar (tbBinA, _this->winName, &_this->binA, 400, ontrackbar_boxes);
        cv::setTrackbarPos (tbBinA, _this->winName, _this->binA);
        std::cout << "createTrackbars: _this->binB=" << _this->binB << std::endl;
        cv::createTrackbar (tbBinB, _this->winName, &_this->binB, 200, ontrackbar_boxes);
        cv::setTrackbarPos (tbBinB, _this->winName, _this->binB);
        cv::createTrackbar (tbNBins, _this->winName, &_this->nBinsTarg, 200, ontrackbar_nbins);
        cv::setTrackbarPos (tbNBins, _this->winName, _this->nBinsTarg);
    }

    static void updateTrackbars (void) {
        std::string tbBinA = "Box A";
        std::string tbBinB = "Box B";
        std::string tbNBins = "Num bins";
        DM* _this = DM::i();
        cv::setTrackbarPos (tbBinA, _this->winName, _this->binA);
        cv::setTrackbarPos (tbBinB, _this->winName, _this->binB);
        cv::setTrackbarPos (tbNBins, _this->winName, _this->nBinsTarg);
    }
};

//! Globally initialise DM instance pointer to NULL
DM* DM::pInstance = 0;
