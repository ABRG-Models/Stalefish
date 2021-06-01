#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <bitset>
#include <fstream>
#include <limits>
#include <utility>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <morph/HdfData.h>
#include <morph/Config.h>
#include <morph/Random.h>
#include <morph/tools.h>
#include "FrameData.h"

// OpenCV functions mostly expect colours in Blue-Green-Red order
#define SF_BLUE     cv::Scalar(255,0,0,10)
#define SF_BLUEISH  cv::Scalar(238,110,67)
#define SF_GREEN    cv::Scalar(0,255,0,10)
#define SF_RED      cv::Scalar(0,0,255,10)
#define SF_YELLOW   cv::Scalar(0,255,255,10)
#define SF_ORANGE   cv::Scalar(25,136,249,10)
#define SF_PURPLE   cv::Scalar(238,121,159,50)
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

    //! The angle which marks a point about half way around the slices, with respect to the origin axes
    double mapAlignAngle = 1.57;

    //! Colour space parameters
    std::string colourmodel = "monochrome";

    //! Allen Developing Mouse Brain Atlas colour mapping parameters
    AllenColourParams acparams;

    // Called by next/previousFrame. Take binA, binB from the frame and change the
    // sliders. Update the fit and refresh boxes. Update the view of boxes/fit line/control points
    void refreshFrame()
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
    //! Set true if we're in 'clear all pending' mode
    bool clearAllPending = false;
    //! Set true if we're in 'export pending' mode
    bool exportPending = false;
    bool importPending = false;

    //! The instance public function. Uses the very short name 'i' to keep code tidy.
    static DM* i()
    {
        if (DM::pInstance == 0) {
            DM::pInstance = new DM;
            DM::i()->init();
        }
        return DM::pInstance;
    }

    //! Initialize by clearing out vFrameData.
    void init() { this->vFrameData.clear(); }

    //! Set up, and add the FramdData object \a fd
    void addFrame (FrameData& fd, const std::string& frameImgFilename, const float& slice_x)
    {
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
        fd.pixels_per_mm = (double)this->pixels_per_mm * this->scaleFactor;
        fd.thickness = this->thickness;
        fd.savePerPixelData = this->savePerPixel;
        fd.saveAutoAlignData = this->saveAutoAlign;
        fd.saveLMAlignData = this->saveLMAlign;
        fd.saveFrameToH5 = this->saveFrameImage;

        fd.rotateLandmarkOne = this->rotateLandmarkOne;
        fd.rotateButAlignLandmarkTwoPlus = this->rotateButAlignLandmarkTwoPlus;

        // Read, opportunistically
        try {
            morph::HdfData d(this->datafile, true); // true for read
            fd.read (d, this->readOldFormat);
            if (fd.flags.test (Mirrored)) {
                fd.mirror_image_only();
            }
            // Update the Bezier fit so that the boxes can be drawn
            fd.updateFit();

        } catch (const std::exception& e) {
            // No problem, just carry on
            std::cout << "Caught: " << e.what() << std::endl;
        }

        // So that freehand loops show up their pixel/signal values in the UI:
        fd.computeFreehandMeans();

        this->vFrameData.push_back (fd);
    }

    /*!
     * Add frame to vFrameData, setting the metadata attributes
     * \a frameImgFilename (The filename for the image), \a slice_x (position in the x
     * dimension)
     */
    void addFrame (const std::string& frameImgFilename, const float& slice_x)
    {
        std::cout << __FUNCTION__ << "(1) called\n";
        // Create an empty FrameData, with no image data as yet.
        FrameData fd(this->bgBlurScreenProportion, this->bgBlurSubtractionOffset,
                     (this->colourmodel == "allen" ? ColourModel::AllenDevMouse : ColourModel::Greyscale),
                     this->acparams);
        this->addFrame (fd, frameImgFilename, slice_x);
    }

    /*!
     * Add frame \a frameImg to vFrameData, setting the metadata attributes
     * \a frameImgFilename (The filename for the image), \a slice_x (position in the x
     * dimension)
     */
    void addFrame (cv::Mat& frameImg, const std::string& frameImgFilename, const float& slice_x)
    {
        std::cout << __FUNCTION__ << "(2) called\n";
        FrameData fd(frameImg, this->bgBlurScreenProportion, this->bgBlurSubtractionOffset,
                     (this->colourmodel == "allen" ? ColourModel::AllenDevMouse : ColourModel::Greyscale),
                     this->acparams);
        this->addFrame (fd, frameImgFilename, slice_x);
    }

    //! Return the size of vFrameData
    unsigned int getNumFrames() const { return this->vFrameData.size(); }

    //! After changing DM::input_mode, apply this to the current frame
    void updateInputMode()
    {
        this->gcf()->ct = this->input_mode;
        this->gcf()->updateFit();
        this->gcf()->refreshBoxes (-this->gcf()->binA, this->gcf()->binB);
    }

    //! Toggle between curve fitting, freehand loop drawing or alignment mark (landmark) input.
    void cycleInputMode()
    {
        if (this->input_mode == InputMode::Bezier || this->input_mode == InputMode::ReverseBezier) {
            this->input_mode = InputMode::Freehand;
        } else if (this->input_mode == InputMode::Freehand) {
            this->input_mode = InputMode::Landmark;
        } else if (this->input_mode == InputMode::Landmark) {
            this->input_mode = InputMode::GlobalLandmark;
        } else if (this->input_mode == InputMode::GlobalLandmark) {
            this->input_mode = InputMode::Circlemark;
        } else if (this->input_mode == InputMode::Circlemark) {
            this->input_mode = InputMode::Axismark;
        } else if (this->input_mode == InputMode::Axismark) {
            this->input_mode = InputMode::Bezier;
        } else {
            // Shouldn't get here...
            this->input_mode = InputMode::Bezier;
        }
        this->updateInputMode();
    }

    //! If on curve mode, toggle between "insert at end" and "insert at start" mode
    void toggleStartEnd()
    {
        if (this->input_mode == InputMode::ReverseBezier) {
            this->input_mode = InputMode::Bezier;
            this->updateInputMode();
        } else if (this->input_mode == InputMode::Bezier) {
            this->input_mode = InputMode::ReverseBezier;
            this->updateInputMode();
        } // else do nothing
    }

    void clearAllCurves()
    {
        std::cout << "Clearing all user-supplied curve points...\n";
        int nfr = DM::i()->getNumFrames();
        for (int f = 0; f < nfr; ++f) {
            this->vFrameData[f].removeAllPoints();
            this->vFrameData[f].refreshBoxes (-this->vFrameData[f].binA, this->vFrameData[f].binB);
        }
    }

    //! Copy the current frame's bin parameters (binA, binB, nBinsTarg) to all the other frames.
    void updateAllBins()
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
    void updateAllFits()
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
    void refreshAllBoxes()
    {
        int nfr = DM::i()->getNumFrames();
        for (int f = 0; f < nfr; ++f) {
            this->vFrameData[f].refreshBoxes (-this->vFrameData[f].binA, this->vFrameData[f].binB);
        }
    }

    //! get current frame. Short name on purpose.
    FrameData* gcf()
    {
        if (!this->vFrameData.empty()) {
            return &(this->vFrameData[this->I]);
        }
        return (FrameData*)0;
    }

    //! Get the current frame number, counting from 1 like a human.
    int getFrameNum() const { return 1+this->I; }

    //! Get a pointer to the persistent Mat img member attribute
    cv::Mat* getImg() { return &(this->img); }
    //! Signal image
    cv::Mat* getSImg() { return &(this->sImg); }

    //! Make the next frame current (or cycle back to the first)
    void nextFrame()
    {
        ++this->I %= this->vFrameData.size();
        this->refreshFrame();
    }

    //! Back up a frame
    void previousFrame()
    {
        this->I = --this->I < 0 ? this->vFrameData.size()-1 : this->I;
        this->refreshFrame();
    }

    //! Clone the current frame into Mat img
    void cloneFrame() {
        this->img = this->vFrameData[this->I].frame.clone();
        FrameData* cf = this->gcf();
        if (cf) {
            this->sImg = cf->frame_signal.clone(); // Don't get an alpha channel with frame_signal
        }
    }

    //! Update the curve fits, re-collect signal data and perform alignment, to ready the project for writing out.
    void writePrep()
    {
        // Call updateAllFits() before writing only to ensure that all the boxes have
        // been refreshed. Seems these are not read out of the .h5 file. Bit of a hack, this.
        this->refreshAllBoxes();

        // Before writing, apply the slice alignment algorithms
        for (auto& f : this->vFrameData) { f.updateAlignments(); }

#ifdef ALIGN_BY_DISTANCE // Using the minimum distance in the y-z plane
        // For the middle slice, find the location of the 'surface box at zero
        // degrees', then for each slice, find the location of the surface box which is
        // closest to the middle slice one.
        size_t middleslice = this->vFrameData.size()/2;
        this->vFrameData[middleslice].setMiddle (this->mapAlignAngle);
        for (auto& f : this->vFrameData) {
            f.setMiddle (this->vFrameData[middleslice]);
        }
#else
        // Align all by angle. Problematic if the origin x axis does not run all the way
        // through the 'centre' of the brain - when, for a given mapAlignAngle, there
        // are >1 brain surfaces available to choose between.
        cv::Point2d am_start_lm;
        cv::Point2d am_start_aa;
        float x_start;
        bool got_start = false;
        cv::Point2d am_end_lm;
        cv::Point2d am_end_aa;
        float x_end;
        bool got_end = false;
        for (auto& f : this->vFrameData) {
            // Work just with one axis mark for now
            if (!f.AM_lmaligned.empty() && got_start && !got_end) {
                std::cout << "Got END axismark on slice " << (f.idx+1) << std::endl;
                am_end_lm = f.AM_lmaligned[0];
                am_end_aa = f.AM_autoaligned[0];
                x_end = f.layer_x;
                got_end = true;
            }
            if (!f.AM_scaled.empty() && !got_start) {
                std::cout << "Got START axismark on slice " << (f.idx+1) << std::endl;
                am_start_lm = f.AM_lmaligned[0];
                am_start_aa = f.AM_autoaligned[0];
                x_start = f.layer_x;
                got_start = true;
            }
            // Just use the first two axis marks we come upon, ignoring the rest.
            if (got_start == true && got_end == true) { break; }
        }

        // If we have alignment marks, then compute them for each frame
        if (got_start == true && got_end == true) {
            std::cout << "Have alignment marks; compute AM_origins from start / end: " << am_start_lm << "/" << am_end_lm << "\n";
            double x1 = (double)x_start;
            double x2 = (double)x_end;
            for (auto& f : this->vFrameData) {

                // First compute for the landmark aligned coordinate system
                double x = (double)f.layer_x;
                double y1 = am_start_lm.x;
                double y2 = am_end_lm.x;
                double z1 = am_start_lm.y;
                double z2 = am_end_lm.y;
                // Compute components of the slope, my and mz.
                double my = (y2-y1)/(x2-x1);
                double mz = (z2-z1)/(x2-x1);
                // Find the offsets, c
                double cy = y1 - my*x1;
                double cz = z1 - mz*x1;
                // Compute y and z for this slice
                double y = my * x + cy;
                double z = mz * x + cz;

                //std::cout << "(lm) y = " << y << ", z = " << z << std::endl;
                cv::Point2d origin_new(y,z);
                //std::cout << "alignmark origin (lm): " << origin_new << std::endl;
                if (f.AM_origins_lmaligned.empty()) {
                    f.AM_origins_lmaligned.push_back (origin_new);
                } else {
                    f.AM_origins_lmaligned[0] = origin_new;
                }
                // Now compute for the autoaligned coordinate system
                y1 = am_start_aa.x;
                y2 = am_end_aa.x;
                z1 = am_start_aa.y;
                z2 = am_end_aa.y;
                // Compute components of the slope, my and mz.
                my = (y2-y1)/(x2-x1);
                mz = (z2-z1)/(x2-x1);
                // Find the offsets, c
                cy = y1 - my*x1;
                cz = z1 - mz*x1;
                // Compute y and z for this slice
                y = my * x + cy;
                z = mz * x + cz;
                //std::cout << "(aa) y = " << y << ", z = " << z << std::endl;
                origin_new.x = y;
                origin_new.y = z;
                //std::cout << "alignmark origin (aa): " << origin_new << std::endl;
                if (f.AM_origins_autoaligned.empty()) {
                    f.AM_origins_autoaligned.push_back (origin_new);
                } else {
                    f.AM_origins_autoaligned[0] = origin_new;
                }
            }
        }

        if (got_start == true && got_end == true) {
            // We have alignment marks; define and then align around that axis
            for (auto& f : this->vFrameData) {
                f.setMiddle (this->mapAlignAngle, f.AM_origins_lmaligned[0], f.AM_origins_autoaligned[0]);
            }
        } else {
            // just align around the origin:
            for (auto& f : this->vFrameData) { f.setMiddle (this->mapAlignAngle); }
        }
#endif
    }

    //! Write only the frames to a separate data file
    void writeMap()
    {
        // First writeFrames anyway
        this->writeFrames();

        std::string mapfile = "map.h5";
        {
            // Get these things from the datafile: For each layer: Layer001/class/layer_x
            morph::HdfData d_in(datafile, true);
            morph::HdfData d_out(mapfile);

            std::vector<float> map_x;
            std::vector<float> map_y_lmalign;
            std::vector<float> map_y_autoalign;
            std::vector<float> map_means;

            d_in.read_contained_vals ("/map/x", map_x);
            d_in.read_contained_vals ("/map/y_lmalign", map_y_lmalign);
            d_in.read_contained_vals ("/map/y_autoalign", map_y_autoalign);
            d_in.read_contained_vals ("/map/means", map_means);

            d_out.add_contained_vals ("/map/x", map_x);
            d_out.add_contained_vals ("/map/y_lmalign", map_y_lmalign);
            d_out.add_contained_vals ("/map/y_autoalign", map_y_autoalign);
            d_out.add_contained_vals ("/map/means", map_means);
        }

        std::cout << "writeMap complete: 2D map written to HDF5" << std::endl;
    }

    //! Write frames to HdfData
    void writeFrames()
    {
        if (this->noFiles == true) { return; }

        this->writePrep();

        morph::HdfData d(this->datafile);

        // Pass in mapAlignAngle for generating the angle maps
        for (auto f : this->vFrameData) { f.write (d, this->mapAlignAngle); }

        std::cout << "Exporting globallandmarks which has size " << this->globalLandmarks.size() << std::endl;
        // /globallandmarks is an index. See this->globalLandmarks. This gives the index
        // of the frame and within that frame the element of FrameData::GLM for each
        // global landmark. Order is that in which user added them.
        d.add_contained_vals ("/global_landmarks", this->globalLandmarks);

        // For each frame also collect 2d Map data. That's /class/layer_x, /signal/postproc/boxes/means, lmalign/flattened/sbox_linear_distance
        std::vector<float> map_x;
        std::vector<float> map_y_lmalign;
        std::vector<float> map_y_autoalign;
        std::vector<float> map_means;
        for (auto f : this->vFrameData) {
            std::string frameName;
            {
                std::stringstream ss;
                ss << "/Frame";
                ss.width(3);
                ss.fill('0');
                ss << (1+f.idx);
                frameName = ss.str();
            }

            std::string pth = frameName + "/signal/postproc/boxes/means";
            std::vector<float> _means;
            d.read_contained_vals(pth.c_str(), _means);

            pth = frameName + "/lmalign/flattened/sbox_linear_distance";
            std::vector<float> _lmalign_lin_distances;
            d.read_contained_vals(pth.c_str(), _lmalign_lin_distances);

            pth = frameName + "/autoalign/flattened/sbox_linear_distance";
            std::vector<float> _autoalign_lin_distances;
            d.read_contained_vals(pth.c_str(), _autoalign_lin_distances);

            pth = frameName + "/class/layer_x";
            float _x;
            d.read_val(pth.c_str(), _x);
            std::vector<float> _xx (_means.size(), _x);

            map_x.insert(map_x.end(), _xx.begin(), _xx.end());
            map_y_lmalign.insert(map_y_lmalign.end(), _lmalign_lin_distances.begin(), _lmalign_lin_distances.end());
            map_y_autoalign.insert(map_y_autoalign.end(), _autoalign_lin_distances.begin(), _autoalign_lin_distances.end());
            map_means.insert(map_means.end(), _means.begin(), _means.end());
        }
        // Now write map/x etc into the HdfData.
        d.add_contained_vals ("/map/x", map_x);
        d.add_contained_vals ("/map/y_lmalign", map_y_lmalign);
        d.add_contained_vals ("/map/y_autoalign", map_y_autoalign);
        d.add_contained_vals ("/map/means", map_means);

        int nf = this->vFrameData.size();
        d.add_val("/nframes", nf);

        // Save the JSON into /config
        d.add_string("/config", this->conf.str());

        std::cout << "writeFrames complete: All frames written to HDF5" << std::endl;
    }

    //! A hardcoded tmp file location to save exported points for import into another project
    const std::string lm_exportfile = "/tmp/landmarks.h5";
    const std::string cp_exportfile = "/tmp/curves.h5";
    const std::string fh_exportfile = "/tmp/freehand.h5";

    //! Depending on the current InputMode, export landmarks, curves OR freehand loops.
    void exportInputModePoints()
    {
        FrameData* cf = DM::i()->gcf();
        if (cf->ct == InputMode::Bezier || cf->ct == InputMode::ReverseBezier) {
            this->exportCurves();
        } else if (cf->ct == InputMode::Freehand) {
            this->exportFreehand();
        } else if (cf->ct == InputMode::Landmark || cf->ct == InputMode::GlobalLandmark
                   || cf->ct == InputMode::Circlemark || cf->ct == InputMode::Axismark) {
            this->exportLandmarks();
        } else {
            std::cerr << "Unknown mode for export\n";
        }
    }

    void exportLandmarks()
    {
        if (this->noFiles == true) { return; }
        this->refreshAllBoxes();
        for (auto& f : this->vFrameData) { f.updateAlignments(); }
        int nf = this->vFrameData.size();
        morph::HdfData d(lm_exportfile);
        for (auto f : this->vFrameData) { f.exportLandmarks (d); }
        std::cout << "Exporting globallandmarks... which has size " << this->globalLandmarks.size() << std::endl;
        d.add_contained_vals ("/global_landmarks", this->globalLandmarks);
        d.add_val("/nframes", nf);
        std::cout << "Exported landmarks to " << lm_exportfile << std::endl;
    }

    void exportFreehand()
    {
        if (this->noFiles == true) { return; }
        this->refreshAllBoxes();
        for (auto& f : this->vFrameData) { f.updateAlignments(); }
        int nf = this->vFrameData.size();
        morph::HdfData d(fh_exportfile);
        for (auto f : this->vFrameData) { f.exportFreehand (d); }
        d.add_val("/nframes", nf);
        std::cout << "Exported freehand loops to " << fh_exportfile << std::endl;
    }

    void exportCurves()
    {
        if (this->noFiles == true) { return; }
        this->refreshAllBoxes();
        for (auto& f : this->vFrameData) { f.updateAlignments(); }
        int nf = this->vFrameData.size();
        morph::HdfData d(cp_exportfile);
        for (auto f : this->vFrameData) { f.exportCurves (d); }
        d.add_val("/nframes", nf);
        std::cout << "Exported curves to " << cp_exportfile << std::endl;
    }

    //! Export all user-supplied point information to files
    void exportUserpoints()
    {
        if (this->noFiles == true) { return; }
        this->refreshAllBoxes();
        for (auto& f : this->vFrameData) { f.updateAlignments(); }
        int nf = this->vFrameData.size();
        {
            morph::HdfData d(lm_exportfile);
            for (auto f : this->vFrameData) { f.exportLandmarks (d); }
            d.add_val("/nframes", nf);
            std::cout << "Exported landmarks to " << lm_exportfile << std::endl;
        }
        {
            morph::HdfData d(cp_exportfile);
            for (auto f : this->vFrameData) { f.exportCurves (d); }
            d.add_val("/nframes", nf);
            std::cout << "Exported curves to " << cp_exportfile << std::endl;
        }
        {
            morph::HdfData d(fh_exportfile);
            for (auto f : this->vFrameData) { f.exportFreehand (d); }
            d.add_val("/nframes", nf);
            std::cout << "Exported freehand loops to " << fh_exportfile << std::endl;
        }
    }

    //! Import landmark information from a file
    void importLandmarks()
    {
        if (this->noFiles == true) { return; }
        try {
            morph::HdfData d(lm_exportfile, true);
            for (auto& f : this->vFrameData) { f.importLandmarks (d); }
            // Now load, if possible, the index data about the global landmarks.
        } catch (...) {
            std::cout << "Failed to read " << lm_exportfile << std::endl;
        }
    }

    //! Import "curve points" from a file
    void importCurves()
    {
        try {
            morph::HdfData d(cp_exportfile, true);
            for (auto& f : this->vFrameData) { f.importCurves (d); }
        } catch (...) {
            std::cout << "Failed to read " << cp_exportfile << std::endl;
        }
    }

    //! Import freehand loops from a file
    void importFreehand()
    {
        if (this->noFiles == true) { return; }
        try {
            morph::HdfData d(fh_exportfile, true);
            for (auto& f : this->vFrameData) { f.importFreehand (d); }
        } catch (...) {
            std::cout << "Failed to read " << fh_exportfile << std::endl;
        }
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
    std::string winName = "StaleFish";
    //! The Gaussian blur window
    std::string blurWin = "";
    //! If true, display the window with the Gaussian blur
    bool showBlurWin = false;
    //! The offset signal window (this has the blurred background subtracted)
    std::string offsWin = "";
    //! If true, display the window with the offset signal
    bool showOffsWin = false;
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
    //! What is the scaling factor to scale the images before saving them in FrameData objects?
    float scaleFactor = 1.0f;

    //! Set true if program was invoked with no .json or .h5 file to open.
    bool noFiles = false;

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

    // Flags to control what data gets saved.
    //! Save data for every single freakin' pixel in every single sample box (instead of means and SDs)
    bool savePerPixel = false;
    //! Save auto (i.e. curve based) alignment location data
    bool saveAutoAlign = true;
    //! Save landmark alignment location data
    bool saveLMAlign = true;
    //! If true, then tell FrameData objects to save a copy of their image data (FrameData::frame)
    bool saveFrameImage = true;
    //! If true, and there are >1 landmark per slice, apply the "rotate slices about
    //! landmark 1" alignment procedure anyway. This rotational alignment is applied by
    //! default if there is ONLY 1 landmark per slice.
    bool rotateLandmarkOne = false;
    //! If true, in "rotate about landmark 1 mode" align the other landmarks, instead of the curves.
    bool rotateButAlignLandmarkTwoPlus = false;

    //! Global landmarks need to be identified in order, so that they can be matched
    //! with a corresponding set of global landmarks from another data set. The pair in
    //! the vector holds as .first, the frame number of that global landmark (counting
    //! from 1) and as .second, the index within FrameData::GLM (counting from 0).
    //! UPDATE: Working with Allen data, which is saggital, means I have to have another
    //! way to order the global landmarks. I can't rely on the slice ordering - with
    //! transverse slices ordered from rostral to caudal - to give me the order of
    //! global landmarks, so I need to preserve the *order in which the user adds the
    //! global landmarks*. So, the meaning of the pair remains unchanged, but the order
    //! of the vector is now important and has to retain the order in which the user
    //! entered the marks.
    std::vector<std::pair<unsigned int, unsigned int>> globalLandmarks;

    //! Adapted from Seb's futil library. Create unique file name for temporary file
    std::string generateRandomFilename (const std::string& prefixPath, const unsigned int numChars = 0)
    {
        std::string rtn(prefixPath);
        unsigned int nc = 8;
        if (numChars > 0) { nc = numChars; }
        morph::RandString rs((size_t)nc, morph::CharGroup::HexLowerCase);
        rtn.append (rs.get());
        return rtn;
    }

    //! Special, minimal setup to show a message to the user when they try to open the app on Mac with no filename
    void noFileSetup()
    {
        std::cout << __FUNCTION__ << " called\n";
        std::string fn("No file provided");
        cv::Mat emptyFrame = cv::Mat(cv::Size(1024,768), CV_8UC3, SF_WHITE);
        this->addFrame (emptyFrame, fn, 0);
        cv::namedWindow (this->winName, cv::WINDOW_NORMAL|cv::WINDOW_FREERATIO);
        // Make sure there's an image in DM to start with
        this->cloneFrame();
        cv::setMouseCallback (this->winName, DM::onmouse, this->getImg());
        int xh = 30;
        int yh = 90;
        int yinc = 40;
        float fontsz = 0.9f;
        cv::Mat* pImg = this->getImg();
        putText (*pImg, std::string("Welcome to Stalefish!"),
                 cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
        fontsz = 0.8f * fontsz;
        yh += 60;
        putText (*pImg, std::string("This application has a very simple user interface and has to be run"),
                 cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
        yh += yinc;
        putText (*pImg, std::string("from the command line. You must specify the JSON or HDF5 file that"),
                 cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
        yh += yinc;
        putText (*pImg, std::string("the program will read. Typically, you'll open a terminal and run:"),
                 cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
        yh += 50;
        putText (*pImg, std::string("  /path/to/stalefish myfile.json"),
                 cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
        yh += 50;
        putText (*pImg, std::string("where myfile.json has been set up as described in the documentation."),
                 cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
        yh += 60;
        putText (*pImg, std::string("Please press 'x' to exit"),
                 cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
    }

    /*!
     * Application setup
     *
     * \param paramsfile Either a path to a .json config file or a path to an HDF5
     * file. If HDF5, then the saved JSON config is searched for within the HDF as
     * '/config'
     */
    void setup (const std::string& paramsfile)
    {
        if (this->noFiles == true) { return this->noFileSetup(); }

        std::string jsonfile (paramsfile);
        bool hdf_for_reading = false;
        // Set the HDF5 data file path based on the .json file path
        std::string::size_type jsonpos = paramsfile.find(".json");
        if (jsonpos == std::string::npos) {
            // If no '.json' in the filename, then test if there's '.h5'. FIXME - should
            // really ID the file content using libmagic, rather than relying on file suffices.
            std::string::size_type h5pos = paramsfile.find(".h5");
            if (h5pos == std::string::npos) {
                std::cerr << paramsfile << " seems to be neither json nor h5; exiting." << std::endl;
                exit (1);
            } else {
                this->datafile = paramsfile;

                // Now get config from h5. Make a temporary file. Fill it with JSON
                // config. Put the path to the tmp file in jsonfile.
                std::string jsoncontent("");
                try {
                    hdf_for_reading = true;
                    morph::HdfData d(this->datafile, true); // true for read
                    d.read_string ("/config", jsoncontent);
                    // What about on Apple. Tmp location there?
                    // Strip directories off paramsfile for tmp file
                    std::string only_paramsfile(paramsfile);
                    morph::Tools::stripUnixPath (only_paramsfile);
                    std::string prefix = "/tmp/json_" + only_paramsfile + "_";
                    jsonfile = "";
                    jsonfile = this->generateRandomFilename (prefix);
                    jsonfile.append(".json");
                    std::ofstream jf;
                    jf.open (jsonfile.c_str(), std::ios::out|std::ios::trunc);
                    if (!jf.is_open()) {
                        std::stringstream ee;
                        ee << "Failed to open temporary file '" << jsonfile << "' for writing.";
                        throw std::runtime_error (ee.str());
                    }
                    jf << jsoncontent;
                    jf.close();

                } catch (const std::exception& e) {
                    std::cerr << "Could not read a configuration from that h5 file (is it too old?): " << e.what() << std::endl;
                    exit (1);
                }
            }

        } else {
            this->datafile = paramsfile.substr (0,jsonpos) + ".h5";
        }

        // Put the data file in the title bar, as that's useful to see with multiple Stalefishies
        this->winName += " " + this->datafile;

        this->conf.init (jsonfile);
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

        this->mapAlignAngle = conf.getDouble ("map_align_angle", 1.57);

        this->default_mode = InputMode::Bezier;

        // The colour space information, if relevant (Allen ISH images)
        // colourmodel - if exists, a string
        this->colourmodel = conf.getString ("colourmodel", "monochrome");
        // colour_rot - array<float, 9>
        const Json::Value cr = conf.getArray ("colour_rot");
        for (unsigned int ii = 0; ii < cr.size(); ++ii) {
            this->acparams.colour_rot[ii] = cr[ii].asFloat();
        }
        // colour_trans - array<float, 3>
        const Json::Value ct = conf.getArray ("colour_trans");
        for (unsigned int ii = 0; ii < ct.size(); ++ii) {
            this->acparams.colour_trans[ii] = ct[ii].asFloat();
        }
        // ellip_axes array<float, 2>
        const Json::Value ea = conf.getArray ("ellip_axes");
        this->acparams.ellip_axes[0] = ea[0].asFloat();
        this->acparams.ellip_axes[1] = ea[1].asFloat();
        // luminosity linear fit parameters
        this->acparams.luminosity_cutoff = conf.getFloat ("luminosity_cutoff", 255.0f);
        this->acparams.luminosity_factor = conf.getFloat ("luminosity_factor", -0.00392f); // -1/255
        // Data-saving parameters - what to save to the HDF5 file
        this->savePerPixel = conf.getBool ("save_per_pixel_data", false);
        this->saveAutoAlign = conf.getBool ("save_auto_align_data", true);
        this->saveLMAlign = conf.getBool ("save_landmark_align_data", true);
        this->saveFrameImage = conf.getBool ("save_frame_image", true);
        // A scaling applied to the original image before it is saved in FrameData::frame
        this->scaleFactor = conf.getFloat ("scaleFactor", 1.0f);

        this->rotateLandmarkOne = conf.getBool ("rotate_landmark_one", false);
        this->rotateButAlignLandmarkTwoPlus = conf.getBool ("rotate_align_landmarks", false);

        // Set false, try to open frames from json; if that fails, then set this
        // true. At that point, we'll try to add frames using in-hdf5 saved data for the
        // frame.
        bool fallbackToInternalFrames = false;

        // Loop over slices, creating a FrameData object for each. BUT if image not
        // findable, allow system to fall back to the image saved in the frame.
        const Json::Value slices = conf.getArray ("slices");
        for (unsigned int si = 0; si < slices.size(); ++si) {
            Json::Value slice = slices[si];
            std::string fn = slice.get ("filename", "unknown").asString();
            float slice_x = slice.get ("x", 0.0).asFloat();

            std::cout << "imread " << fn << std::endl;
            cv::Mat frame = cv::imread (fn.c_str(), cv::IMREAD_COLOR);
            if (frame.empty()) {
                std::cout << "Could not open or find the image '"
                          << fn << "'. Will attempt fall-back to image data stored in the h5 file." << std::endl;
                fallbackToInternalFrames = true;
            }

            if (fallbackToInternalFrames == true) {
                // Try loading the frame exclusively from file
                this->addFrame (fn, slice_x);
                // Reset the bool, incase we can read the next one.
                fallbackToInternalFrames = false;

            } else {
                if (std::abs(this->scaleFactor - 1.0f) > std::numeric_limits<float>::epsilon()) {
                    std::cout << "rescaling frame to scaleFactor: " << scaleFactor << std::endl;

                    cv::Size scaledSize = cv::Size(std::round(frame.cols * scaleFactor),
                                                   std::round(frame.rows * scaleFactor));
                    cv::Mat scaledFrame = cv::Mat(scaledSize, frame.type());
                    cv::resize (frame, scaledFrame, scaledSize,
                                scaleFactor, scaleFactor, cv::INTER_LINEAR);

                    frame.release(); // free original frame since we have resized it

                    this->addFrame (scaledFrame, fn, slice_x);
                } else {
                    // if we are at the default scale factor do not do anything
                    this->addFrame (frame, fn, slice_x);
                }
            }
        }

        // Having read frames, read in any global information, which is currently just globalLandmarks
        if (hdf_for_reading == true) {
            morph::HdfData d(this->datafile, true); // true for read
            d.read_contained_vals ("/global_landmarks", this->globalLandmarks);
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

        // green circle under the cursor indicates curve mode
        if (cf->ct == InputMode::Bezier || cf->ct == InputMode::ReverseBezier) {
            circle (*pImg, pt, 5, SF_GREEN, 1);
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
        if ((cf->ct == InputMode::Bezier || cf->ct == InputMode::ReverseBezier)
            && cf->flags.test(ShowUsers) == true) {
            // draw the "candidate" point set (for adding to the end of the curve):
            if (cf->PP.empty() || (!cf->PP.empty() && cf->P.size() > 1)) {
                for (size_t ii=0; ii<cf->P.size(); ii++) {
                    circle (*pImg, cf->P[ii], 5, SF_GREEN, -1);
                    if (ii) { line (*pImg, cf->P[ii-1], cf->P[ii], SF_GREEN, 1, cv::LINE_AA); }
                }
            }
            // draw the "candidate" point set (for adding to the *start* of the curve):
            if (cf->PP.empty() || (!cf->PP.empty() && cf->sP.size() > 1)) {
                for (size_t ii=0; ii<cf->sP.size(); ii++) {
                    circle (*pImg, cf->sP[ii], 5, SF_GREEN, -1);
                    if (ii) { line (*pImg, cf->sP[ii-1], cf->sP[ii], SF_GREEN, 1, cv::LINE_AA); }
                }
            }
            // also draw a thin line to the cursor position
            if (cf->ct == InputMode::Bezier) {
                if ((cf->PP.empty() && cf->P.size() > 0)
                    || (!cf->PP.empty() && cf->P.size() > 1)) {
                    line (*pImg, cf->P[cf->P.size()-1], pt, SF_GREEN, 1, cv::LINE_AA);
                }
            } else if (cf->ct == InputMode::ReverseBezier) {
                if ((cf->PP.empty() && cf->sP.size() > 0)
                    || (!cf->PP.empty() && cf->sP.size() > 1)) {
                    line (*pImg, cf->sP[0], pt, SF_GREEN, 1, cv::LINE_AA);
                }
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

    void draw_boundary (const std::vector<cv::Point>& vp, cv::Mat* _pImg, const cv::Scalar& colour)
    {
        double alpha = 0.8;
        for (size_t ii=0; ii<vp.size(); ii++) {
            cv::Point cur = vp[ii];
            cv::Mat roi = (*_pImg)(cv::Rect(cur.x, cur.y, 1, 1));
            cv::Mat color(roi.size(), CV_8UC3, colour);
            cv::addWeighted (color, alpha, roi, 1.0 - alpha , 0.0, roi, CV_8UC3);
        }
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
            if (!cf->FLB[j].empty()) {
                draw_boundary (cf->FLB[j], pImg, SF_BLUE);
            }
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

    //! Draw the axis marks
    void draw_axismarks (const cv::Point& pt)
    {
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        FrameData* cf = _this->gcf();

        // circle under the cursor
        if (cf->ct == InputMode::Axismark) {
            circle (*pImg, pt, 7, SF_ORANGE, 1);
        }

        // Draw circles for the axismarks, with a number next to each one.
        cv::Point toffset(8,5); // a text offset
        for (size_t ii=0; ii<cf->AM.size(); ii++) {
            circle (*pImg, cf->AM[ii], 5, SF_ORANGE, -1);
            std::stringstream ss;
            ss << 'a' << (1+ii);
            cv::Point tpt(cf->AM[ii]);
            putText (*pImg, ss.str(), tpt+toffset, cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
        }
    }

    //! Draw the global land marks. The order of global landmarks is defined by order in
    //! which the user placed them.
    void draw_global_landmarks (const cv::Point& pt)
    {
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        FrameData* cf = _this->gcf();

        // circle under the cursor
        if (cf->ct == InputMode::GlobalLandmark) { circle (*pImg, pt, 7, SF_BLUEISH, 1); }

        // Iterate through the frames, to find out what starting index is for the global
        // landmarks on this slice
        size_t start_gl = 1;
        for (const std::pair<unsigned int, unsigned int>& glm : this->globalLandmarks) {
            int glm_fidx = glm.first-1;
            if (glm_fidx == cf->idx) { break; }
            start_gl += this->vFrameData[glm_fidx].GLM.size();
        }

        // Draw circles for the axismarks, with a number next to each one.
        cv::Point toffset(8,5); // a text offset
        for (size_t ii=0; ii<cf->GLM.size(); ii++) {
            circle (*pImg, cf->GLM[ii], 5, SF_BLUEISH, -1);
            std::stringstream ss;
            ss << "gl" << (start_gl+ii);
            cv::Point tpt(cf->GLM[ii]);
            putText (*pImg, ss.str(), tpt+toffset, cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
        }
    }

    //! Input mode for drawing the numbered alignment marks
    void draw_landmarks (const cv::Point& pt)
    {
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        FrameData* cf = _this->gcf();

        // circle under the cursor
        if (cf->ct == InputMode::Landmark) { circle (*pImg, pt, 7, SF_BLACK, 1); }

        // Draw circles for the landmarks, with a number next to each one.
        cv::Point toffset(8,5); // a text offset
        int lm_ok = cf->landmarkCheck(); // True if all landmarks are present and correct
        for (size_t ii=0; ii<cf->LM.size(); ii++) {
            circle (*pImg, cf->LM[ii], 5, (static_cast<int>(ii)<lm_ok ? SF_BLACK : SF_RED), -1);

            std::stringstream ss;
            ss << (1+ii);
            cv::Point tpt(cf->LM[ii]);
            putText (*pImg, ss.str(), tpt+toffset, cv::FONT_HERSHEY_SIMPLEX, 0.8, SF_BLACK, 1, cv::LINE_AA);
        }
    }

    //! Input mode for entering landmarks by defining 3 points on a circle.
    void draw_circlemarks (const cv::Point& pt)
    {
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        FrameData* cf = _this->gcf();

        // circle under the cursor
        if (cf->ct == InputMode::Circlemark) {
            circle (*pImg, pt, 7, SF_PURPLE, 1);
            // Cross hares too
            line (*pImg, pt-cv::Point(0,6), pt+cv::Point(0,6), SF_BLACK, 1);
            line (*pImg, pt-cv::Point(6,0), pt+cv::Point(6,0), SF_BLACK, 1);
        }

        // Draw any entries in CM
        auto cmi = cf->CM.begin();
        while (cmi != cf->CM.end()) {
            std::vector<cv::Point> pts = cmi->second;
            for (size_t ii = 0; ii < pts.size(); ++ii) {
                circle (*pImg, pts[ii], 4, SF_PURPLE, -1);
            }
            ++cmi;
        }

        // Draw any entries in CM_points
        for (size_t ii = 0; ii < cf->CM_points.size(); ++ii) {
            circle (*pImg, cf->CM_points[ii], 4, SF_RED, -1);
        }
    }

    void removeLastThing()
    {
        // If we're removing a global landmark, have to delete from globalLandmarks first.
        FrameData* cf = this->gcf();
        if (cf->ct == InputMode::GlobalLandmark) {
            if (!cf->GLM.empty()) {
                std::pair<unsigned int, unsigned int> cmp = std::make_pair (this->getFrameNum(), cf->GLM.size()-1);
                std::vector<std::pair<unsigned int, unsigned int>>::iterator gli = this->globalLandmarks.begin();
                while (gli != this->globalLandmarks.end()) {
                    if (*gli == cmp) {
                        std::cout << "Erasing global landmark " << cmp.first << "," << cmp.second << std::endl;
                        gli = this->globalLandmarks.erase (gli);
                    } else {
                        gli++;
                    }
                }
            }
        }
        cf->removeLastThing();
    }

    //! Actions to take on a mouse user-interface event
    static void onmouse (int event, int x, int y, int flags, void* param)
    {
        // Make copies of pointers to neaten up the code, below
        DM* _this = DM::i();
        cv::Mat* pImg = _this->getImg();
        cv::Mat* sImg = _this->getSImg();
        FrameData* cf = _this->gcf();

        // If we're in no file mode, don't respond to the mouse
        if (_this->noFiles == true) {
            imshow (_this->winName, *pImg);
            return;
        }

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
            } else if (cf->ct == InputMode::ReverseBezier) {
                // If we have some already-registered curves in PP, and sP is empty, we
                // have to add *2* points.
                if (!cf->PP.empty() && cf->sP.empty()) { cf->sP.push_front (cf->PP.front().front()); }
                cf->sP.push_front (pt);
                cf->setShowUsers(true);
            } else if (cf->ct == InputMode::Freehand) {
                cf->addToFL (pt);
            } else if (cf->ct == InputMode::Landmark) {
                cf->LM.push_back (pt);
            } else if (cf->ct == InputMode::GlobalLandmark) {
                cf->addGlobalLandmark (pt);
                _this->globalLandmarks.push_back (std::make_pair (_this->getFrameNum(), cf->GLM.size()-1));
            } else if (cf->ct == InputMode::Circlemark) {
                cf->addCirclepoint (pt);
            } else if (cf->ct == InputMode::Axismark) {
                cf->addAxismark (pt);
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
        _this->draw_global_landmarks (pt);
        _this->draw_axismarks (pt);
        _this->draw_circlemarks (pt);

        std::stringstream ss;
        int xh = 30;
        ss << "Frame: " << _this->getFrameNum() << "/" << _this->getNumFrames()
           << " " << cf->getFitInfo() << ".";
        ss << " Range: " << cf->frame_maxmin.second << "," << cf->frame_maxmin.first;
        ss.precision(3);
        // Float font size based on image size (frame width)?
        float fwidth = (float)cf->frame.cols;
        float fontsz = fwidth / 1727.0f;
        putText (*pImg, ss.str(), cv::Point(xh,30), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);

        std::stringstream ss2;
        ss2.precision(3);
        ss2 << " Signal range: " << cf->frame_signal_maxmin.second << "," << cf->frame_signal_maxmin.first
            << " (using blur offset: " << _this->bgBlurSubtractionOffset << ")";
        putText (*sImg, ss2.str(), cv::Point(xh,30), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_WHITE, 1, cv::LINE_AA);

        if (_this->clearAllPending == true) {
            std::stringstream ss3;
            ss3 << "Clear curves on ALL frames? (press 'C' to confirm, 'Esc' to cancel)";
            putText (*pImg, ss3.str(), cv::Point(xh,80), cv::FONT_HERSHEY_SIMPLEX, 1.2*fontsz, SF_BLACK, 1, cv::LINE_AA);
            putText (*sImg, ss3.str(), cv::Point(xh,80), cv::FONT_HERSHEY_SIMPLEX, 1.2*fontsz, SF_WHITE, 1, cv::LINE_AA);
        }
        else if (_this->exportPending == true) {
            std::stringstream ss3;
            ss3 << "Export data? (press key again to confirm, 'Esc' to cancel)";
            putText (*pImg, ss3.str(), cv::Point(xh,80), cv::FONT_HERSHEY_SIMPLEX, 1.2*fontsz, SF_BLACK, 1, cv::LINE_AA);
            putText (*sImg, ss3.str(), cv::Point(xh,80), cv::FONT_HERSHEY_SIMPLEX, 1.2*fontsz, SF_WHITE, 1, cv::LINE_AA);
        }
        else if (_this->importPending == true) {
            std::stringstream ss3;
            ss3 << "IMPORT data? (press key again to confirm, 'Esc' to cancel)";
            putText (*pImg, ss3.str(), cv::Point(xh,80), cv::FONT_HERSHEY_SIMPLEX, 1.2*fontsz, SF_BLACK, 1, cv::LINE_AA);
            putText (*sImg, ss3.str(), cv::Point(xh,80), cv::FONT_HERSHEY_SIMPLEX, 1.2*fontsz, SF_WHITE, 1, cv::LINE_AA);
        }

        // h for help
        std::string hs("Press 'h' for help");
        putText (*pImg, hs, cv::Point(cf->frame.cols-300,cf->frame.rows-20), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);

        // Draw text with cursor coordinates on bottom left of screen in a small font
        std::stringstream css;
        css << "px(y=" << x << ", z=" << y << ") "
            << "mm[x=" << cf->layer_x << ", y=" << (x/cf->pixels_per_mm) << ", z=" << (y/cf->pixels_per_mm) << "]";
        putText (*pImg, css.str(), cv::Point(xh,cf->frame.rows-20), cv::FONT_HERSHEY_SIMPLEX, fontsz/2.0f, SF_BLACK, 1, cv::LINE_AA);

        int yh = 90;
        int yinc = 40;
        if (_this->flags.test(AppShowHelp)) {
            putText (*pImg, std::string("Use the sliders to control the bin parameters"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("1:   Toggle Bezier controls    2: Toggle user points"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("3:   Toggle the fit line    4: Toggle the bins"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("Spc: Next curve"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("c:   Cancel last point/freehand region"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("C:   Delete ALL curves"),
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
            hh << "w:   Save to file: " << _this->datafile << " (W: Save this and also the 2D map to ./map.h5)";
            putText (*pImg, hh.str(),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("o:   Draw mode (Curve/freehand/landmark)"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("s:   Toggle add points to curve at start/end"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("k:   Export points to files in /tmp"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("p:   Export landmark OR curves OR freehand to file in /tmp (depends on Draw mode)"),
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("l:   Import landmarks from ") + _this->lm_exportfile,
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("i:   Import curve points from ") + _this->cp_exportfile,
                     cv::Point(xh,yh), cv::FONT_HERSHEY_SIMPLEX, fontsz, SF_BLACK, 1, cv::LINE_AA);
            yh += yinc;
            putText (*pImg, std::string("j:   Import freehand loops from ") + _this->fh_exportfile,
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
        std::cout << __FUNCTION__ << " called\n";
        FrameData* cf = DM::i()->gcf();
        unsigned int nbt = DM::i()->nBinsTarg;
        if (nbt < 2) { nbt = 2; }
        cf->setBins (nbt);
        cf->setShowBoxes (true);
        cf->updateFit();
        cf->refreshBoxes (-cf->binA, cf->binB);
        DM::onmouse (cv::EVENT_MOUSEMOVE, -1, -1, 0, NULL);
    }

    static void createTrackbars()
    {
        // Set up trackbars. Have to do this for each frame
        std::string tbBinA = "Box A";
        std::string tbBinB = "Box B";
        std::string tbNBins = "Num bins";
        DM* _this = DM::i();
        // binA slider goes from -200 to +200
        cv::createTrackbar (tbBinA, _this->winName, &_this->binA, 400, ontrackbar_boxes);
        cv::setTrackbarPos (tbBinA, _this->winName, _this->binA);
        // binB slider from 0 to 300
        cv::createTrackbar (tbBinB, _this->winName, &_this->binB, 300, ontrackbar_boxes);
        cv::setTrackbarPos (tbBinB, _this->winName, _this->binB);
        cv::createTrackbar (tbNBins, _this->winName, &_this->nBinsTarg, 200, ontrackbar_nbins);
        cv::setTrackbarPos (tbNBins, _this->winName, _this->nBinsTarg);
    }

    static void updateTrackbars()
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
