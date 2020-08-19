#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <bitset>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <morph/HdfData.h>
#include <morph/Config.h>
#include "FrameData.h"

// OpenCV functions mostly expect colours in Blue-Green-Red order
#define SF_BLUE     cv::Scalar(255,0,0,10)
#define SF_GREEN    cv::Scalar(0,255,0,10)
#define SF_RED      cv::Scalar(0,0,255,10)
#define SF_YELLOW   cv::Scalar(0,255,255,10)
#define SF_BLACK    cv::Scalar(0,0,0)
#define SF_WHITE    cv::Scalar(255,255,255)
#define SF_C1       cv::Scalar(238,121,159) // mediumpurple2
#define SF_C2       cv::Scalar(238,58,178) // darkorchid2
#define SF_TRANS    cv::Scalar(255,0,0,100)

//! OpenCV sliders can't be negative, so we have an offset.
#define BIN_A_OFFSET 200

//! Flags used by the application for application-wide state
enum AppFlag {
    AppShowBoxes, // Show the yellow boxes?
    AppShowUsers, // Show the user points?
    AppShowCtrls, // Show the ctrl points of the fits?
    AppShowFits,  // Show the fits?
    AppShowHelp
};

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
    //! The signal image
    cv::Mat sImg;
    //! The thickness of each brain slice in mm
    float thickness = 0.05f;
    //! The application configuration
    morph::Config conf;

    //! Colour space parameters
    std::string colourmodel = "monochrome";
    std::array<float, 9> colour_rot; // Colour space rotation
    std::array<float, 3> colour_trans; // Colour space pre-translation
    std::array<float, 2> ellip_axes; // red-green ellipse for "elliptical tube of expressing colours
    float luminosity_factor; // The slope of the linear luminosity vs signal fit.
    float luminosity_cutoff; // at what luminosity does the signal cut off to zero?

    // Called by next/previousFrame. Take binA, binB from the frame and change the
    // sliders. Update the fit and refresh boxes. Update the view of boxes/fit line/control points
    void refreshFrame (void)
    {
        this->binA = this->gcf()->binA+BIN_A_OFFSET;
        this->binB = this->gcf()->binB;
        this->nBinsTarg = this->gcf()->getBins();
        DM::updateTrackbars();
        this->gcf()->ct = this->input_mode;
        this->gcf()->updateFit();
        this->gcf()->refreshBoxes (-(this->binA-BIN_A_OFFSET), this->binB);
        this->gcf()->setShowCtrls (this->flags.test(AppShowCtrls));
        this->gcf()->setShowUsers (this->flags.test(AppShowUsers));
        this->gcf()->setShowFits (this->flags.test(AppShowFits));
        this->gcf()->setShowBoxes (this->flags.test(AppShowBoxes));
    }

public:
    //! A bit set containing flags to track application state
    std::bitset<8> flags;
    //! What's the global input mode?
    InputMode input_mode = InputMode::Bezier;

    //! The instance public function. Uses the very short name 'i' to keep code tidy.
    static DM* i (void)
    {
        if (DM::pInstance == 0) {
            DM::pInstance = new DM;
            DM::i()->init();
        }
        return DM::pInstance;
    }

    //! Initialize by clearing out vFrameData.
    void init (void) { this->vFrameData.clear(); }

    /*!
     * Add frame \a frameImg to vFrameData, setting the metadata attributes
     * \a frameImgFilename (The filename for the image), \a slice_x (position in the x
     * dimension) and \a ppm (pixels per mm; the scale).
     */
    void addFrame (cv::Mat& frameImg, const std::string& frameImgFilename, const float& slice_x)
    {
        FrameData fd(frameImg, this->bgBlurScreenProportion, this->bgBlurSubtractionOffset);
        fd.ct = this->default_mode;
        fd.filename = frameImgFilename;
        fd.setParentStack (&this->vFrameData);
        // Increment layer index. Best might be to use JSON info for layer positions as
        // they are unlikely always to be in perfect increments.
        if (this->vFrameData.empty()) {
            // First frame; not setting previous
            fd.idx = 0;
        } else {
            // Subsequent frame; setting previous
            fd.idx = this->vFrameData.back().idx + 1;
            // Set pointer to previous so slices can be aligned during FrameData::write or updateFit
            fd.setPrevious (this->vFrameData.back().idx);
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

        // Read, opportunistically
        try {
            morph::HdfData d(this->datafile, true); // true for read
            fd.read (d, this->readOldFormat);
            if (fd.flags.test (Mirrored)) {
                fd.mirror_image_only();
            }
            // Update the Bezier fit so that the boxes can be drawn
            fd.updateFit();
            // Update landmarks, so they exist in the vFrameData?

        } catch (...) {
            // No problem, just carry on
        }

        this->vFrameData.push_back (fd);
    }

    //! Return the size of vFrameData
    unsigned int getNumFrames (void) const { return this->vFrameData.size(); }

    //! Toggle between curve fitting, freehand loop drawing or alignment mark (landmark) input.
    void cycleInputMode()
    {
        if (this->input_mode == InputMode::Landmark) {
            this->input_mode = InputMode::Bezier;
        } else if (this->input_mode == InputMode::Bezier) {
            this->input_mode = InputMode::Freehand;
        } else if (this->input_mode == InputMode::Freehand) {
            this->input_mode = InputMode::Landmark;
        } else {
            // Shouldn't get here...
            this->input_mode = InputMode::Bezier;
        }
    }

    //! Copy the current frame's bin parameters (binA, binB, nBinsTarg) to all the other frames.
    void updateAllBins (void)
    {
        int nfr = DM::i()->getNumFrames();
        int idx = DM::i()->gcf()->idx;
        std::cout << nfr << " frames; current is " << idx << "\n";
        for (int f = 0; f < nfr; ++f) {
            if (idx == f) {
                continue;
            }
            this->vFrameData[f].binA = this->vFrameData[idx].binA;
            this->vFrameData[f].binB = this->vFrameData[idx].binB;
            std::cout << "current frame's nBinsTarg = " << this->vFrameData[idx].getBins() << "\n";
            std::cout << "current DM::nBinsTarg = " << this->nBinsTarg << "\n";
            this->vFrameData[f].setBins(this->vFrameData[idx].getBins());
            this->vFrameData[f].updateFit();
            this->vFrameData[f].refreshBoxes (-this->vFrameData[f].binA, this->vFrameData[f].binB);
        }
        std::cout << "End of updateAllBins()\n";
    }

    //! Update all fits - i.e. for every frame in the stack
    void updateAllFits (void)
    {
        int nfr = DM::i()->getNumFrames();
        for (int f = 0; f < nfr; ++f) {
            this->setShowFits (true);
            this->vFrameData[f].setShowFits (true);
            this->setShowBoxes (true);
            this->vFrameData[f].setShowBoxes (true);
            this->vFrameData[f].updateFit();
            this->vFrameData[f].refreshBoxes (-this->vFrameData[f].binA, this->vFrameData[f].binB);
        }
    }

    //! Call before write to ensure boxes are all created from the current fits.
    void refreshAllBoxes (void)
    {
        int nfr = DM::i()->getNumFrames();
        for (int f = 0; f < nfr; ++f) {
            this->vFrameData[f].refreshBoxes (-this->vFrameData[f].binA, this->vFrameData[f].binB);
        }
    }

    //! get current frame. Short name on purpose.
    FrameData* gcf (void)
    {
        if (!this->vFrameData.empty()) {
            return &(this->vFrameData[this->I]);
        }
        return (FrameData*)0;
    }

    //! Get the current frame number, counting from 1 like a human.
    int getFrameNum (void) const { return 1+this->I; }

    //! Get a pointer to the persistent Mat img member attribute
    cv::Mat* getImg (void) { return &(this->img); }
    //! Signal image
    cv::Mat* getSImg (void) { return &(this->sImg); }

    //! Make the next frame current (or cycle back to the first)
    void nextFrame (void)
    {
        ++this->I %= this->vFrameData.size();
        this->refreshFrame();
    }

    //! Back up a frame
    void previousFrame (void)
    {
        this->I = --this->I < 0 ? this->vFrameData.size()-1 : this->I;
        this->refreshFrame();
    }

    //! Clone the current frame into Mat img
    void cloneFrame (void) {
        this->img = this->vFrameData[this->I].frame.clone();
        FrameData* cf = this->gcf();
        if (cf) {
            this->sImg = cf->frame_signalU.clone(); // Don't get an alpha channel with frame_signal, unlike frame_bgoffU
        }
    }

    //! Write frames to HdfData
    void writeFrames (void)
    {
        // Call updateAllFits() before writing only to ensure that all the boxes have
        // been refreshed. Seems these are not read out of the .h5 file. Bit of a hack, this.
        this->refreshAllBoxes();

        // Before writing, apply the slice alignment algorithms
        for (auto& f : this->vFrameData) {
            f.updateAlignments();
        }

        morph::HdfData d(this->datafile);
        for (auto f : this->vFrameData) {
            f.write (d);
            // Also build up an "overall" data store of the bins
        }
        int nf = this->vFrameData.size();
        d.add_val("/nframes", nf);
        std::cout << "writeFrames complete: All frames written to HDF5" << std::endl;
    }

    //! Toogle showHelp
    void toggleShowHelp() { this->flags[AppShowHelp] = this->flags.test(AppShowHelp) ? false : true; }
    void setShowHelp (bool t) { this->flags[AppShowHelp] = t; }

    void toggleShowBoxes() { this->flags[AppShowBoxes] = this->flags.test(AppShowBoxes) ? false : true; }
    void setShowBoxes (bool t) { this->flags[AppShowBoxes] = t; }

    void toggleShowFits() { this->flags[AppShowFits] = this->flags.test(AppShowFits) ? false : true; }
    void setShowFits (bool t) { this->flags[AppShowFits] = t; }

    void toggleShowUsers() { this->flags[AppShowUsers] = this->flags.test(AppShowUsers) ? false : true; }
    void setShowUsers (bool t) { this->flags[AppShowUsers] = t; }

    void toggleShowCtrls() { this->flags[AppShowCtrls] = this->flags.test(AppShowCtrls) ? false : true; }
    void setShowCtrls (bool t) { this->flags[AppShowCtrls] = t; }

    //! The application window name
    const std::string winName = "StaleFish";
    //! The Gaussian blur window
    std::string blurWin = "";
    //! If true, display the window with the Gaussian blur
    bool showBlurWin = false;
    //! The offset signal window (this has the blurred background subtracted)
    std::string offsWin = "";
    //! If true, display the window with the offset signal
    bool showOffsWin = true;
    //! Toggle for the blur window
    void toggleBlurWindow() { this->showBlurWin = this->showBlurWin ? false : true; }
    //! Toggle for the offset (signal) window
    void toggleOffsWindow() { this->showOffsWin = this->showOffsWin ? false : true; }
    //! Are we on the first window drawing call?
    bool firstCall = true;
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

    //! Set true to read in old data format, to be written out in new format.
    bool readOldFormat = false;

    //! Which drawing mode should the application start in? Bezier by default. JSON can
    //! be used to modify.
    InputMode default_mode = InputMode::Bezier;

    //! The sigma for the Gaussian used to blur the image to get the overall background
    //! luminance is the framewidth in pixels multiplied by this number.
    double bgBlurScreenProportion = 0.1667;

    //! A subtraction offset used when subtracting blurred background signal from image
    float bgBlurSubtractionOffset = 255.0;

    //! Application setup
    void setup (const std::string& paramsfile)
    {
        // Set the HDF5 data file path based on the .json file path
        std::string::size_type jsonpos = paramsfile.find(".json");
        if (jsonpos == std::string::npos) {
            this->datafile = paramsfile + ".h5";
        } else {
            this->datafile = paramsfile.substr (0,jsonpos) + ".h5";
        }

        this->conf.init (paramsfile);
        if (!this->conf.ready) {
            std::cerr << "Error setting up JSON config: "
                      << this->conf.emsg << ", exiting." << std::endl;
            exit (1);
        }

        // Set the scale from JSON, too
        this->pixels_per_mm = conf.getFloat ("pixels_per_mm", 100.0f);
        this->thickness = conf.getFloat ("thickness", 0.05f);

        // Set parameters for background offsetting.
        this->bgBlurScreenProportion = conf.getDouble ("bg_blur_screen_proportion", 0.1667);
        this->bgBlurSubtractionOffset = conf.getDouble ("bg_blur_subtraction_offset", 255.0f);

        this->default_mode = InputMode::Bezier;

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
        this->luminosity_cutoff = conf.getFloat ("luminosity_cutoff", 255.0f);
        this->luminosity_factor = conf.getFloat ("luminosity_factor", -0.00392f); // -1/255

        // Loop over slices, creating a FrameData object for each.
        const Json::Value slices = conf.getArray ("slices");
        for (unsigned int i = 0; i < slices.size(); ++i) {
            Json::Value slice = slices[i];
            std::string fn = slice.get ("filename", "unknown").asString();
            float slice_x = slice.get ("x", 0.0).asFloat();

            std::cout << "imread " << fn << std::endl;
            cv::Mat frame = cv::imread (fn.c_str(), cv::IMREAD_COLOR);
            if (frame.empty()) {
                std::cout <<  "Could not open or find the image '"
                          << fn << "', exiting." << std::endl;
                exit (1);
            }

            // scaling routine (pull scale factor from config json)
            float scaleFactor = conf.getFloat("scaleFactor", 1.0f);

            if (scaleFactor != 1.0f) {
                std::cout << "rescaling frame to scaleFactor: " << scaleFactor << std::endl;

                cv::Size scaledSize = cv::Size(std::round(frame.cols * scaleFactor),
                                               std::round(frame.rows * scaleFactor));
                cv::Mat scaledFrame = cv::Mat(scaledSize, frame.type());
                cv::resize (frame, scaledFrame, scaledSize,
                            scaleFactor, scaleFactor, cv::INTER_LINEAR);

                frame.release(); // free original frame since we have resized it

                this->addFrame(scaledFrame, fn, slice_x);
            } else {
                // if we are at the default scale factor do not do anything
                this->addFrame(frame, fn, slice_x);
            }
        }

        cv::namedWindow (this->winName, cv::WINDOW_NORMAL|cv::WINDOW_FREERATIO);
        // Make sure there's an image in DM to start with
        this->cloneFrame();
        cv::setMouseCallback (this->winName, DM::onmouse, this->getImg());
        // Init current frame with binA, binB and nBinsTarg taken from the current
        // frame... But... is that information stored? Yes, it is.
        this->binA = this->gcf()->binA+BIN_A_OFFSET;
        this->binB = this->gcf()->binB;
        this->nBinsTarg = this->gcf()->getBins();

        this->setShowUsers (this->gcf()->getShowUsers());
        this->setShowFits (this->gcf()->getShowFits());
        this->setShowBoxes (this->gcf()->getShowBoxes());
        this->setShowCtrls (this->gcf()->getShowCtrls());

        DM::createTrackbars();
        this->gcf()->refreshBoxes (-(this->binA-BIN_A_OFFSET), this->binB);
    }

    //! In Bezier or Polynomial modes, draw curves and users points
    void draw_curves (const cv::Point& pt)
    {
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        FrameData* cf = _this->gcf();
        // Signal image comes straight out of the FrameData
        cv::Mat* sImg = _this->getSImg();

        // red circle under the cursor
        if (cf->ct == InputMode::Bezier) {
            circle (*pImg, pt, 5, SF_RED, 1);
        }

        if (cf->flags.test(ShowUsers) == true) {
            // First the lines in the preceding PP point-sets:
            for (size_t j=0; j<cf->PP.size(); j++) {
                cv::Scalar linecol = j%2 ? SF_RED : SF_BLUE;
                for (size_t ii=0; ii<cf->PP[j].size(); ii++) {
                    circle (*pImg, cf->PP[j][ii], 5, linecol, -1);
                    if (ii) { line (*pImg, cf->PP[j][ii-1], cf->PP[j][ii], linecol, 2, cv::LINE_AA); }
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
                line (*pImg, ps, pe, SF_GREEN, 1, cv::LINE_AA);
                cv::Point ps2(ctrls[ctrls.size()-2].first, ctrls[ctrls.size()-2].second);
                cv::Point pe2(ctrls[ctrls.size()-1].first, ctrls[ctrls.size()-1].second);
                line (*pImg, ps2, pe2, SF_GREEN, 1, cv::LINE_AA);

                j++;
            }
        }

        // This is the set of green user points that will be the next Bezier curve section
        if (cf->ct == InputMode::Bezier && cf->flags.test(ShowUsers) == true) {
            // Then draw the current point set:
            if (cf->PP.empty() || (!cf->PP.empty() && cf->P.size() > 1)) {
                for (size_t ii=0; ii<cf->P.size(); ii++) {
                    circle (*pImg, cf->P[ii], 5, SF_GREEN, -1);
                    if (ii) { line (*pImg, cf->P[ii-1], cf->P[ii], SF_GREEN, 1, cv::LINE_AA); }
                }
            }
            // also draw a thin line to the cursor position
            if ((cf->PP.empty() && cf->P.size() > 0)
                || (!cf->PP.empty() && cf->P.size() > 1)) {
                line (*pImg, cf->P[cf->P.size()-1], pt, SF_GREEN, 1, cv::LINE_AA);
            }
        }

        // This is the fit line
        if (cf->flags.test(ShowFits) == true) {
            for (size_t ii=1; ii<cf->fitted.size(); ii++) {
                line (*pImg, cf->fitted[ii-1], cf->fitted[ii], SF_GREEN, 2, cv::LINE_AA);
                // line (*sImg, cf->fitted[ii-1], cf->fitted[ii], SF_BLACK, 2, cv::LINE_AA);
            }
        }

        if (cf->flags.test(ShowBoxes) == true) {
            // The bins; pointsInner to pointsOuter
            for (size_t ii=0; ii<cf->pointsInner.size(); ii++) {
                line (*pImg, cf->pointsInner[ii], cf->pointsOuter[ii], SF_YELLOW, 1, cv::LINE_AA);
                line (*sImg, cf->pointsInner[ii], cf->pointsOuter[ii], SF_BLACK, 1, cv::LINE_AA);

                if (ii > 0) {
                    line (*pImg, cf->pointsInner[ii-1], cf->pointsInner[ii], SF_YELLOW, 1, cv::LINE_AA);
                    line (*pImg, cf->pointsOuter[ii-1], cf->pointsOuter[ii], SF_YELLOW, 1, cv::LINE_AA);

                    line (*sImg, cf->pointsInner[ii-1], cf->pointsInner[ii], SF_BLACK, 1, cv::LINE_AA);
                    line (*sImg, cf->pointsOuter[ii-1], cf->pointsOuter[ii], SF_BLACK, 1, cv::LINE_AA);
                }
            }
        }
    }

    //! On \a _pImg, draw the region specified in \a vp, using \a colour. Also add text
    //! showing the mean signal/luminance (\a themean)
    void draw_region (const std::vector<cv::Point>& vp, cv::Mat* _pImg, const float& themean, const cv::Scalar& colour)
    {
        double alpha = 0.3;
        int xmean = 0;
        int ymean = 0;
        for (size_t ii=0; ii<vp.size(); ii++) {
            // This fills area with transparent blue. First get region from *_pImg
            cv::Point cur = vp[ii];
            cv::Mat roi = (*_pImg)(cv::Rect(cur.x, cur.y, 1, 1));
            // Compute mean x,y? for text position?
            xmean += cur.x;
            ymean += cur.y;
            // Create a colour
            cv::Mat color(roi.size(), CV_8UC3, colour);
            cv::addWeighted (color, alpha, roi, 1.0 - alpha , 0.0, roi, CV_8UC3);
        }
        xmean /= vp.size();
        ymean /= vp.size();
        // Add text for FL_signal_means
        std::stringstream flm;
        flm << themean;
        cv::Point tpt(xmean, ymean); // Could use extents_FL here.
        putText (*_pImg, flm.str(), tpt, cv::FONT_HERSHEY_SIMPLEX, 0.5, SF_BLACK, 1, cv::LINE_AA);
    }

    //! Draw freehand loops when in InputMode::Freehand mode
    void draw_freehand (const cv::Point& pt)
    {
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        cv::Mat* sImg = _this->getSImg();
        FrameData* cf = _this->gcf();

        // blue? circle under the cursor
        if (cf->ct == InputMode::Freehand) {
            circle (*pImg, pt, 3, SF_BLUE, 1);
        }

        // Draw the existing regions
        for (size_t j=0; j<cf->FLE.size(); j++) {
            if (!cf->FLE[j].empty()) {
                float themean = cf->FL_signal_means.size() > j ? cf->FL_signal_means[j] : 0.0f;
                float thepixelmean = cf->FL_pixel_means.size() > j ? cf->FL_pixel_means[j] : 0;
                draw_region (cf->FLE[j], pImg, thepixelmean, SF_BLUE);
                draw_region (cf->FLE[j], sImg, themean, SF_BLACK);
            }
        }

        // Draw the current point set in green:
        for (size_t ii=0; ii<cf->FL.size(); ii++) {
            rectangle (*pImg, cf->FL[ii], cf->FL[ii], SF_GREEN, 1);
            rectangle (*sImg, cf->FL[ii], cf->FL[ii], SF_GREEN, 1);
        }

#ifdef DEBUG_INSIDE_OUTSIDE_BOUNDARY
        for (size_t ii=0; ii<cf->inside_FL.size(); ii++) {
            rectangle (*pImg, cf->inside_FL[ii], cf->inside_FL[ii], SF_RED, 1);
        }
        for (size_t ii=0; ii<cf->outside_FL.size(); ii++) {
            rectangle (*pImg, cf->outside_FL[ii], cf->outside_FL[ii], SF_BLUE, 1);
        }
#endif

#ifdef DEBUG_EXTENTS
        line (*pImg, cf->extents_FL[0], cf->extents_FL[0]+cv::Point(0,3), SF_BLACK, 1);
        line (*pImg, cf->extents_FL[0], cf->extents_FL[0]+cv::Point(3,0), SF_BLACK, 1);
        line (*pImg, cf->extents_FL[1], cf->extents_FL[1]+cv::Point(0,-3), SF_BLACK, 1);
        line (*pImg, cf->extents_FL[1], cf->extents_FL[1]+cv::Point(-3,0), SF_BLACK, 1);
#endif
    }

    //! Input mode for drawing the numbered alignment marks
    void draw_landmarks (const cv::Point& pt)
    {
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        FrameData* cf = _this->gcf();

        // circle under the cursor
        if (cf->ct == InputMode::Landmark) {
            circle (*pImg, pt, 7, SF_BLACK, 1);
        }

        // Draw circles for the landmarks, with a number next to each one.
        cv::Point toffset(8,5); // a text offset
        bool lm_ok = cf->landmarkCheck(); // True if all landmarks are present and correct
        for (size_t ii=0; ii<cf->LM.size(); ii++) {
            circle (*pImg, cf->LM[ii], 5, (lm_ok ? SF_BLACK : SF_RED), -1);

            std::stringstream ss;
            ss << (1+ii);
            cv::Point tpt(cf->LM[ii]);
            putText (*pImg, ss.str(), tpt+toffset, cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
        }
    }

    //! Actions to take on a mouse user-interface event
    static void onmouse (int event, int x, int y, int flags, void* param)
    {
        // Make copies of pointers to neaten up the code, below
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        cv::Mat* sImg = _this->getSImg();
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
            if (cf->ct == InputMode::Bezier) {
                cf->P.push_back (pt);
                cf->setShowUsers(true);
            } else if (cf->ct == InputMode::Freehand) {
                cf->addToFL (pt);
            } else if (cf->ct == InputMode::Landmark) {
                cf->LM.push_back (pt);
            }
        } else if (event == cv::EVENT_LBUTTONUP) {
            cf->loopFinished = false;

        } else if (event == cv::EVENT_MOUSEMOVE
                   && (flags & cv::EVENT_FLAG_LBUTTON) == cv::EVENT_FLAG_LBUTTON
                   && cf->ct == InputMode::Freehand
                   && cf->loopFinished == false) {
            // Now button is down, want to add any pixel that the mouse moves over
            cf->addToFL (pt);
        }

        _this->cloneFrame();

        // Code for drawing stuff when we're in a curve-fitting mode. Draw both curves
        // *and* freehand loops on the screen, always. However, we will need to know
        // what *input* mode we're in, to draw curve points or freehand loops.
        _this->draw_curves (pt);
        _this->draw_freehand (pt);
        _this->draw_landmarks (pt);

        std::stringstream ss;
        int xh = 30;
        ss << "Frame: " << _this->getFrameNum() << "/" << _this->getNumFrames()
           << " " << cf->getFitInfo() << ". 'h' to toggle help.";
        ss << " Range: " << cf->frame_maxmin.second << "," << cf->frame_maxmin.first;
        ss.precision(3);
        ss << " (bm:"<< cf->blurmeanU << ")";
        // Float font size based on image size (frame width)?
        float fwidth = (float)cf->frame.cols;
        float fontsz = fwidth / 1727.0f;
        putText (*pImg, ss.str(), cv::Point(xh,30), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);

        std::stringstream ss2;
        ss2.precision(3);
        ss2 << " Signal range: " << cf->frame_signal_maxmin.second << "," << cf->frame_signal_maxmin.first
            << " (using blur offset: " << _this->bgBlurSubtractionOffset << ")";
        putText (*sImg, ss2.str(), cv::Point(xh,30), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_WHITE, 1, cv::LINE_AA);

        int yh = 90;
        int yinc = 40;
        if (_this->flags.test(AppShowHelp)) {
            putText (*pImg, std::string("Use the sliders to control the bin parameters"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("1:   Toggle Bezier controls"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("2:   Toggle user points"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("3:   Toggle the fit line"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("4:   Toggle the bins"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("Spc: Next curve"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("c:   Cancel last point/freehand region"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("f:   Update the fit"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("F:   Update fit (all frames)"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("B:   Copy current bin params to all frames"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            std::stringstream hh;
            hh << "w:   Save to file: " << _this->datafile;
            putText (*pImg, hh.str(),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("o:   Draw mode (Curve/freehand/landmark)"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("n:   Next frame"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("b:   Back to previous frame"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("m:   Mirror this frame"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("r:   Toggle blur window"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("e:   Toggle signal window"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("x:   Exit the program"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
        }

        // Always show the main window
        imshow (_this->winName, *pImg);
        _this->firstCall = false;

        // Optionally show blurry window...
        if (_this->showBlurWin == true) {
            if (_this->blurWin == "") {
                _this->blurWin = "blurWin";
                cv::namedWindow (_this->blurWin, cv::WINDOW_NORMAL|cv::WINDOW_FREERATIO);
                cv::setWindowTitle (_this->blurWin, "Gaussian blur");
                imshow (_this->blurWin, *cf->getBlur());
            } // else nothing to do, window is already showing (unless user closed it!)
        } else {
            if (_this->blurWin == "blurWin") {
                cv::destroyWindow (_this->blurWin);
                _this->blurWin = "";
            } // else nothing to do
        }

        // ...and the offset window
        if (_this->showOffsWin == true) {
            //if (_this->offsWin == "") { // Always re-draw offset window, as it has items on itx
            _this->offsWin = "offsWin";
            cv::namedWindow (_this->offsWin, cv::WINDOW_NORMAL|cv::WINDOW_FREERATIO);
            cv::setWindowTitle (_this->offsWin, "mRNA signal");
            imshow (_this->offsWin, *sImg);
            //}
        } else {
            if (_this->offsWin == "offsWin") {
                cv::destroyWindow (_this->offsWin);
                _this->offsWin = "";
            }
        }
    }

    //! On any trackbar changing, refresh the boxes
    static void ontrackbar_boxes (int val, void*)
    {
        DM* _this = DM::i();
        FrameData* cf = _this->gcf();
        cf->binA = _this->binA-BIN_A_OFFSET;
        cf->binB = _this->binB;
        cf->setShowBoxes (true);
        cf->refreshBoxes (-cf->binA, cf->binB);
        DM::onmouse (cv::EVENT_MOUSEMOVE, -1, -1, 0, NULL);
    }

    //! On trackbar change, refresh the size of the boxes
    static void ontrackbar_nbins (int val, void*)
    {
        FrameData* cf = DM::i()->gcf();
        unsigned int nbt = DM::i()->nBinsTarg;
        if (nbt < 2) { nbt = 2; }
        cf->setBins (nbt);
        cf->setShowBoxes (true);
        cf->updateFit();
        cf->refreshBoxes (-cf->binA, cf->binB);
        DM::onmouse (cv::EVENT_MOUSEMOVE, -1, -1, 0, NULL);
    }

    static void createTrackbars (void)
    {
        // Set up trackbars. Have to do this for each frame
        std::string tbBinA = "Box A";
        std::string tbBinB = "Box B";
        std::string tbNBins = "Num bins";
        DM* _this = DM::i();
        //std::cout << "createTrackbars: _this->binA=" << _this->binA << std::endl;
        cv::createTrackbar (tbBinA, _this->winName, &_this->binA, 400, ontrackbar_boxes);
        cv::setTrackbarPos (tbBinA, _this->winName, _this->binA);
        //std::cout << "createTrackbars: _this->binB=" << _this->binB << std::endl;
        cv::createTrackbar (tbBinB, _this->winName, &_this->binB, 200, ontrackbar_boxes);
        cv::setTrackbarPos (tbBinB, _this->winName, _this->binB);
        cv::createTrackbar (tbNBins, _this->winName, &_this->nBinsTarg, 200, ontrackbar_nbins);
        cv::setTrackbarPos (tbNBins, _this->winName, _this->nBinsTarg);
    }

    static void updateTrackbars (void)
    {
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
