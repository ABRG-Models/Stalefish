#pragma once

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <stdexcept>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <iostream>
#include <utility>
#include <bitset>
#include <limits>
#include "tools.h"
#include <morph/BezCurvePath.h>
#include <morph/BezCurve.h>
#include <morph/BezCoord.h>
#include <morph/HdfData.h>
#include <morph/MathAlgo.h>
#include <morph/NM_Simplex.h>
#include <morph/MathConst.h>
#include <morph/Winder.h>

// This is no longer really "curve type" but rather "input mode". So in Bezier mode, you
// add points for the curve fitting; in freehand mode you draw a loop, and in align
// mode, you give landmarks for slice alignment.
enum class InputMode {
    Poly,      // Variable order polynomial (note: the existence of this mode is ignored)
    Bezier,    // Cubic Bezier
    Freehand,  // A freehand drawn loop enclosing a region // This will go - Freehand loops to be visible over curve.
    Landmark   // User provides alignment landmark locations on each slice
};

// What sort of colour model is in use?
enum class ColourModel {
    Greyscale,
    AllenDevMouse
};

enum FrameFlag {
    ShowBoxes, // Show the yellow boxes?
    ShowUsers, // Show the user points?
    ShowCtrls, // Show the ctrl points of the fits?
    ShowFits,  // Show the fits?
    Mirrored,  // Is the image mirrored L-R?
    Flipped    // Is the image flipped U-D?
};

/*!
 * A class to hold a cortical section image, user-supplied cortex edge points and the
 * resulting fit (either polynomial or Bezier curved). Also stores information about
 * freehand drawn loops (whose content mean luminance can be saved out) and
 * user-supplied landmarks, allowing for the alignment of multiple brain slices.
 */
class FrameData
{
    //! Private attributes
private:
    //! The 'previous' frame in the stack of frames (index into)
    int previous = -1;
    std::vector<FrameData>* parentStack;

    //! Number of bins to create for the fit (one less than nFit)
    int nBins;
    //! Number of points to create in the fit
    int nFit;

    //! Public attributes
public:
    //! What curve type?
    InputMode ct = InputMode::Bezier;

    // Polynomial fit specific attributes

    //! Order of the polynomial fit
    int polyOrder;
    //! A rotation is applied to the points before calling polyfit, after which fitted
    //! points are rotated back again.
    double theta = 0.0;
    //! Min and max value used as input to polyfit
    double minX, maxX;
    //! The parameters of a polyfit, as returned by polyfit()
    std::vector<double> pf;

    // Bezier curve attributes

    //! A Bezier curve path to fit the cortex.
    morph::BezCurvePath<double> bcp;

    // Attributes which pertain either to polynomial or Bezier curves

    //! The vector of user-supplied points from which to make a curve fit.
    std::vector<cv::Point> P;
    //! A vector of vectors of points for multi-section Bezier curves
    std::vector<std::vector<cv::Point>> PP;
    //! Index into PP
    int pp_idx = 0;

    // Landmark attributes

    //! The landmark points for this frame.
    std::vector<cv::Point> LM;

    //! The means computed for the boxes. This is now "mean_signal" really, as the pixel values
    //! (monochrome or colour) are now converted with the colour space parameters into a signal.
    std::vector<double> means;
    //! The raw values for each box as a vector of floating points numbers for each box.
    std::vector<std::vector<float>> boxes_raw;
    //! Raw colours of boxes in RGB
    std::vector<std::vector<std::array<float, 3>>> boxes_raw_bgr;

    //! Target number of bins; used by bins slider
    int nBinsTarg;
    //! The bin lengths, set with a slider.
    int binA = 0;
    int binB = 100;

    //! A set of points created from the fit
    std::vector<cv::Point> fitted;
    //! The centroid of fitted.
    cv::Point2d fit_centroid;
    //! This holds offset and scaled fitted points: (fitted - fit_centroid) * pixels_per_mm
    std::vector<cv::Point2d> fitted_offset;
    //! This holds the fitted_offset points after they have been rotated to be in line
    //! with fitted_rotated points in the 'previous' frame. 'in line' means the
    //! smallest sum-of-square distances between the two fitted_rotated sets. Depends
    //! on nBins being the same in each.
    std::vector<cv::Point2d> fitted_rotated;
    //! For point in fitted, the tangent at that location
    std::vector<cv::Point2d> tangents;
    //! For point in fitted, the normal at that location
    std::vector<cv::Point2d> normals;
    //! The axes for the polynomial fit
    std::vector<cv::Point> axis;
    //! Coefficients for the polynomial fit axis
    std::vector<double> axiscoefs;
    //! origins for the lines making the box sides (some distance from the curve)
    std::vector<cv::Point> pointsInner;
    //! endpoints for the lines making the box sides (a greater distance from the curve)
    std::vector<cv::Point> pointsOuter;
    //! The boxes that are drawn and from which to sample the gene expression
    std::vector<std::vector<cv::Point> > boxes;

    //! A vector of user-supplied points for the Freehand drawn loop
    std::vector<cv::Point> FL;
    std::array<cv::Point, 2> extents_FL; // Extnents of the loop FL
    //! vector of vectors containing the points enclosed by the path FL
    std::vector<std::vector<cv::Point>> FLE;

    //! The mean luminance of each freehand loop enclosed region in FLE.
    std::vector<double> FL_means;
    //! The raw values for each region as a vector of floating points numbers for each one.
    std::vector<std::vector<float>> FL_raw;
    //! Raw values for each region in colour
    std::vector<std::vector<std::array<float, 3>>> FL_raw_bgr;

    //! A bit set containing flags
    std::bitset<8> flags;
    //! The image data, required when sampling the image in one of the boxes. CV_8UC3
    cv::Mat frame;
    //! A copy of the image data in float format (in range 0 to 1, rather than 0 to 255). CV_32FC3
    cv::Mat frameF;
    //! A blurred copy of the image data. CV_32FC3.
    cv::Mat blurred;
    //! Sets the width of the blurring Gaussian's sigma
    double bgBlurScreenProportion = 0.1667;
    //! An offset used when subtracting blurred BG from image.
    float bgBlurSubtractionOffset = 255.0f;
    //! The frame, with the blurred background offset. CV_32FC3.
    cv::Mat frame_bgoff;
    //! The frame image filename from which frame was loaded. Stored so it can be
    //! recorded when writing out.
    std::string filename;
    //! The 'x' position of the brain slice in this frame (coordinates in the plane of
    //! the slice are y/z (units: mm).
    float layer_x = 0.0f;
    //! The thickness of this slice (mm)
    float thickness = 0.05;
    //! The scaling factor
    double pixels_per_mm = 100.0;
    //! The index of the frame
    int idx;

    // Colour space parameters.

    //! What kind of colour model is in use?
    ColourModel cmodel = ColourModel::Greyscale;
    //! Colour space rotation to apply to [b g r] colour vectors
    std::array<float, 9> colour_rot;
    //! Colour space pre-translation. [x y z] is [b g r]
    std::array<float, 3> colour_trans;
    //! red-green ellipse for "elliptical tube of expressing colours
    std::array<float, 2> ellip_axes;
    //! The slope of the linear luminosity vs signal fit.
    float luminosity_factor;
    //! at what luminosity does the signal cut off to zero?
    float luminosity_cutoff;

public:
    FrameData()
    {
        throw std::runtime_error ("Default constructor is not allowed");
    }
    //! Constructor initializes default values
    FrameData (const cv::Mat& fr,
               const double _bgBlurScreenProportion,
               const float _bgBlurSubtractionOffset)
    {
        // init previous to null.
        this->previous = -1;
        this->parentStack = (std::vector<FrameData>*)0;
        this->frame = fr.clone();
        // Scale and convert frame to float format
        this->frame.convertTo (this->frameF, CV_32FC3, 1/255.0);
        this->axiscoefs.resize (2, 0.0);
        this->axis.resize (2);
        // NB: Init these before the next three resize() calls
        this->setBins (100);
        this->polyOrder = 3;
        // Init flags
        this->flags.set (ShowFits);
        this->flags.set (ShowUsers);
        this->flags.set (ShowBoxes);

        // Make a blurred copy of the floating point format frame, for estimating lighting background
        this->bgBlurScreenProportion = _bgBlurScreenProportion;
        this->blurred = cv::Mat::zeros (this->frameF.rows, this->frameF.cols, CV_32FC3);
        cv::Size ksz;
        std::cout << "FrameData constructor: bgBlurScreenProportion = "
                  << this->bgBlurScreenProportion << std::endl;
        ksz.width = this->frameF.cols * 2.0 * this->bgBlurScreenProportion;
        ksz.width += (ksz.width%2 == 1) ? 0 : 1; // ensure ksz.width is odd
        ksz.height = this->frameF.rows/3;
        ksz.height += (ksz.height%2 == 1) ? 0 : 1;
        double sigma = (double)this->frameF.cols * this->bgBlurScreenProportion;
        cv::GaussianBlur (this->frameF, this->blurred, ksz, sigma);

        // Now subtract the blur from the original
        cv::Mat tmp1;
        // An offset so that we don't lose very small amounts of signal above the
        // background when we subtract the blurred version of the image. To become a
        // parameter for the user to modify at application level.
        if (_bgBlurSubtractionOffset < 0.0f || _bgBlurSubtractionOffset > 255.0f) {
            throw std::runtime_error ("The bg blur subtraction offset should be in range [0,255]");
        }
        this->bgBlurSubtractionOffset = _bgBlurSubtractionOffset;
        cv::subtract (this->bgBlurSubtractionOffset/255.0f, this->blurred,
                      tmp1, cv::noArray(), CV_32FC3);
        // show max mins (For debugging)
        this->showMaxMin (this->blurred, "this->blurred");
        this->showMaxMin (tmp1, "tmp1 (const-blurred)");
        // Add tmp1 to this->frameF to get frame_bgoff.
        cv::add (this->frameF, tmp1, this->frame_bgoff, cv::noArray(), CV_32FC3);
        this->showMaxMin (this->frame_bgoff, "frame_bgoff");
    }

    //! Set the number of bins and update the size of the various containers
    void setBins (unsigned int num)
    {
        if (num > 5000) {
            throw std::runtime_error ("Too many bins...");
        }
        this->nBins = num;
        this->nBinsTarg = num;
        this->nFit = num + 1;
        this->fitted.resize (this->nFit);
        this->fitted_offset.resize (this->nFit);
        this->fitted_rotated.resize (this->nFit);
        this->pointsInner.resize (this->nFit);
        this->pointsOuter.resize (this->nFit);
        this->tangents.resize (this->nFit);
        this->normals.resize (this->nFit);
    }

    //! Show the max and the min of a
    void showMaxMin (const cv::Mat& m, const std::string& matlabel = "(unknown)")
    {
        float minm = 100.0f;
        float maxm = -100.0f;
        for (int r = 0; r < m.rows; ++r) {
            for (int c = 0; c < m.cols; ++c) {
                //std::cout << "Frame("<<r<<","<<c<<") = " << m.at< cv::Vec<float, 3> >(r,c) << "\n";
                float val = (float)m.at< cv::Vec<float, 3> >(r,c)[0];
                minm = val < minm ? val : minm;
                maxm = val > maxm ? val : maxm;
            }
        }
        std::cout << "The matrix " << matlabel << "  has min/max: " << minm << "/" << maxm << std::endl;
    }

    //! Getter for the blurred image
    cv::Mat* getBlur() { return &this->blurred; }
    cv::Mat* getFrameOffs() { return &this->frame_bgoff; }

    //! Getter for nBins
    int getNBins() { return this->nBins; }
    int getNBinsTarg() { return this->nBinsTarg; }

    //! Setter for previous
    void setPrevious (int prev) { this->previous = prev; }
    void setParentStack (std::vector<FrameData>* parentSt) { this->parentStack = parentSt; }

    //! Get information about the fit
    std::string getFitInfo() const
    {
        std::stringstream ss;
        if (this->ct == InputMode::Poly) {
            ss << "Poly order: " << this->polyOrder << ", Bins: " << this->nBins;

        } else /*if (this->ct == InputMode::Bezier)*/ {
            std::stringstream bb;
            bool first = true;
            for (auto cv : this->bcp.curves) {
                if (first) {
                    bb << cv.getOrder();
                    first = false;
                } else {
                    bb << "/" << cv.getOrder();
                }
            }
            ss << "Bezier order: " << bb.str() << ", Bins: " << this->nBins;
        }

        if (this->ct == InputMode::Poly || this->ct == InputMode::Bezier) {
            ss << ". Curve mode";
        } else if (this->ct == InputMode::Freehand) {
            // Get any fit info for a freehand loop (e.g. is it contiguous; how many pixels)
            ss << ". Freehand mode";
        } else if (this->ct == InputMode::Landmark) {
            ss << ". Landmark mode";
        } else {
            ss << ". unknown mode";
        }
        return ss.str();
    }

    //! This is a candidate for MathAlgo; filling squares in between two randomly
    //! chosen squares on a grid. Each square is 1x1 on a grid, with its centre
    //! specified by its cv::Point.  Draw a line between firstSquare and
    //! endSquare. Fill in all pixels which are crossed by the line. Do this with a
    //! recursive algorithm, as it's easy and we're very unlikely to exceed the
    //! recursion limit.
    void fillFL (cv::Point& firstSquare, const cv::Point& endSquare)
    {
        // Finished when the last element of FL is pt.
        if (!this->FL.empty() && this->FL.back() == endSquare) {
            return;
        }

        // Start at firstSquare
        if (firstSquare == endSquare) {
            // There's nothing to do
            return;
        }

        // firstSquare is not the same as endSquare, so fill in between them.  Which of the 8
        // adjoining squares contains the line specified by firstSquare.xy and m?  What's
        // the intersection of the line and the square perimeter of firstSquare? Note, we
        // also need to track currentSquare.

        // Don't need to test if .x==.x AND .y==.y

        // First, push back the firstSquare itself
        this->FL.push_back (firstSquare);

        int xdiff = std::abs(firstSquare.x - endSquare.x);
        int ydiff = std::abs(firstSquare.y - endSquare.y);

        if (firstSquare.y == endSquare.y) {
            // Dirn is east or west
            firstSquare.x += (endSquare.x > firstSquare.x ? 1 : -1);
            return this->fillFL (firstSquare, endSquare);

        } else if (firstSquare.x == endSquare.x) {
            // Dirn is n or s
            firstSquare.y += (endSquare.y > firstSquare.y ? 1 : -1);
            return this->fillFL (firstSquare, endSquare);

        } else if (firstSquare.y < endSquare.y) {
            // SouthWest or SouthEast
            if (firstSquare.x < endSquare.x) {
                // SouthEast
                if (xdiff > ydiff) {
                    // Mark E and move SE
                    firstSquare.x += 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.y += 1;

                } else if (xdiff < ydiff) {
                    // Mark S and move SE
                    firstSquare.y += 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.x += 1;

                } else { // xdiff == ydiff
                    // Move SE only
                    firstSquare.x += 1;
                    firstSquare.y += 1;
                }

            } else {
                // SouthWest
                if (xdiff > ydiff) {
                    // Mark W and move SW
                    firstSquare.x -= 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.y += 1;

                } else if (xdiff < ydiff) {
                    // Mark S and move SW
                    firstSquare.y += 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.x -= 1;

                } else { // xdiff == ydiff
                    // Move SE only
                    firstSquare.x -= 1;
                    firstSquare.y += 1;
                }
            }
        } else {
            // std::cout << "Dirn is NW/NE, computing...\n";
            // North West/East
            if (firstSquare.x < endSquare.x) {
                // NorthEast
                if (xdiff > ydiff) {
                    // Mark E and move NE
                    firstSquare.x += 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.y -= 1;

                } else if (xdiff < ydiff) {
                    // Mark N and move NE
                    firstSquare.y -= 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.x += 1;

                } else { // xdiff == ydiff
                    // Move NE only
                    firstSquare.x += 1;
                    firstSquare.y -= 1;
                }

            } else {
                // NorthWest
                if (xdiff > ydiff) {
                    // Mark W and move NW
                    firstSquare.x -= 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.y -= 1;

                } else if (xdiff < ydiff) {
                    // Mark N and move NW
                    firstSquare.y -= 1;
                    this->FL.push_back (firstSquare);
                    firstSquare.x -= 1;

                } else { // xdiff == ydiff
                    // Move NW only
                    firstSquare.x -= 1;
                    firstSquare.y -= 1;
                }
            }
        }

        // Recurse
        return this->fillFL (firstSquare, endSquare);
    }

    //! Find rectangular region enclosing loop, return as two coordinates of top left
    //! (minimum x and y, with y going down) and bottom right (max x and y) corners
    std::array<cv::Point, 2> getExtents (const std::vector<cv::Point>& loop)
    {
        std::array<cv::Point, 2> extents;
        extents[0] =  {10000000, 10000000};   // MIN values for x and y
        extents[1] = {-10000000,-10000000}; // MAX values for x and y
        for (auto p : loop) {
            if (p.x < extents[0].x) {
                extents[0].x = p.x;
            }
            if (p.y < extents[0].y) {
                extents[0].y = p.y;
            }
            if (p.x > extents[1].x) {
                extents[1].x = p.x;
            }
            if (p.y > extents[1].y) {
                extents[1].y = p.y;
            }
        }
        return extents;
    }

    //! Set true if a loop has been drawn, then completed. Reset to false once mouse button is released
    bool loopFinished = false;

    //! Find all pixels enclosed by the pixels in this->FL which define a loop
    std::vector<cv::Point> getEnclosedByFL()
    {
        // FIXME: Prefer not to have to uniquify here:
        auto last = std::unique(this->FL.begin(), this->FL.end());
        this->FL.erase(last, this->FL.end());

        std::vector<cv::Point> rtn;

        // First, find extents of the loop
        this->extents_FL = this->getExtents (this->FL);

        // Create a winder object to compute winding numbers
        morph::Winder w (this->FL);

        // It's perhaps inefficient to compute the winding number of EVERY pixel here,
        // but I'll leave it for now (computers are fast).
        for (int x = this->extents_FL[0].x; x <= this->extents_FL[1].x; ++x) {
            for (int y = this->extents_FL[0].y; y <= this->extents_FL[1].y; ++y) {
                cv::Point px (x, y);
                auto inloop = std::find (this->FL.begin(), this->FL.end(), px);
                if (inloop == this->FL.end()) {
                    // Compute winding number
                    int winding_number = w.wind (px);
                    //std::cout << "Winding number of pixel " << px << " = " << winding_number << std::endl;
                    if (winding_number != 0) {
                        rtn.push_back (px);
                    }
                } // else current pixel is a member of the loop itself.
            }
        }
        return rtn;
    }

    //! Add a pixel that was under the mouse pointer to the freehand points FL. Also
    //! add the pixels between the last pixel in FL and endSquare.
    void addToFL (cv::Point& pt)
    {
        // If FL empty then task is simple:
        if (this->FL.empty()) {
             this->FL.push_back (pt);
             return;
        }

        auto existing = std::find (this->FL.begin(), this->FL.end(), pt);
        if (existing != this->FL.end() && existing != this->FL.begin()) {
            // If pt is in FL already, then we closed the loop
            // std::cout << "An existing point has been found, so close the loop.\n";
            // We want to cut off any extraneous pixels in FL.
            // Delete from this->FL from start to one before existing.
            this->FL.erase (this->FL.begin(), existing);
            // Fill the gap between last point on curve and the existing curve point.
            cv::Point fp = this->FL.back();
            this->fillFL (fp, pt);
            // Get enclosed and add to FLE:
            std::vector<cv::Point> inside = this->getEnclosedByFL();
            this->FLE.push_back (inside);
            this->FL.clear();
            this->loopFinished = true;

        } else { // The new point pt is NOT in FL already

            std::pair<float, float> p1, p2;
            float snap_threshold = 3.0;

            // Otherwise, check if we can close the loop...
            p1.first = (float)this->FL.begin()->x;
            p1.second = (float)this->FL.begin()->y;
            p2.first = (float)pt.x;
            p2.second = (float)pt.y;

            if (this->FL.size() > 20 // Avoid joining a just-started loop
                && morph::MathAlgo::distance<float>(p1, p2) < snap_threshold) {
                // The point pt close to the start of the freehand-drawn loop, so "autosnap" to it.
                pt = *(this->FL.begin());
                //std::cout << "Joining the loop!\n";
                cv::Point fp = this->FL.back();
                this->fillFL (fp, pt);
                // Get the enclosed points and add:
                std::vector<cv::Point> inside = this->getEnclosedByFL();
                this->FLE.push_back (inside);
                this->FL.clear();
                this->loopFinished = true;

            } else {
                // Can't close the loop; just add to it
                cv::Point fp = this->FL.back();
                this->fillFL (fp, pt);
            }
        }
    }

    //! Remove the last freehand drawn region or the last point, depending on mode
    void removeLastThing()
    {
        if (this->ct == InputMode::Freehand) {
            this->removeLastRegion();
        } else if (this->ct == InputMode::Landmark) {
            this->removeLastLandmark();
        } else {
            this->removeLastPoint();
        }
    }

    //! Remove the last freehand drawn region
    void removeLastRegion()
    {
        // If there's a half-finished boundary, get rid of that first:
        if (!this->FL.empty()) {
            this->FL.clear();
        } else {
            // Otherwise, remove the last completed freehand drawn region
            if (!this->FLE.empty()) {
                this->FLE.pop_back();
            }
        }
    }

    //! Remove the last user point
    void removeLastPoint()
    {
        if (this->PP.empty() && this->P.size() == 1) {
            // Normal behaviour, just remove point from P
            this->P.pop_back();

        } else if (!this->PP.empty() && this->P.size() == 1) {
            // Remove point from P and...
            this->P.pop_back();
            // Because it is the same locn, the last point from PP.back(), too
            this->removeLastPoint();

        } else if (!this->P.empty()) {
            this->P.pop_back();

        } else {
            // P is empty, go to previous curve and remove a point from that
            if (this->ct == InputMode::Bezier && this->pp_idx>0) {
                this->P = this->PP[--this->pp_idx];
                this->PP.pop_back();
                this->P.pop_back();
            }
        }
    }

    //! Remove the last landmark coordinate
    void removeLastLandmark()
    {
        if (!this->LM.empty()) {
            this->LM.pop_back();
        }
    }

    //! In Bezier mode, store the current set of user points (P) into PP and clear P.
    void nextCurve()
    {
        if (this->ct == InputMode::Poly) {
            // no op.
            return;
        }
        // Don't add unless at least 3 points to fit:
        if (this->P.size() < 3) {
            return;
        }
        this->PP.push_back (this->P);
        this->P.clear();
        this->P.push_back (this->PP.back().back());
        this->pp_idx++;
    }

    //! Read important data from file
    void read (morph::HdfData& df, bool oldformat=false)
    {
        // Note this file assumes idx has been set for the frame.
        std::string frameName = this->getFrameName();

        if (oldformat == true) {
            throw std::runtime_error ("Note: there is currently no old format conversion code.");
            // NB: This code is left as a place holder in case we need to read in an old
            // format and save in a new format.
        }

        std::string dname = frameName + "/class/polyOrder";

        df.read_val (dname.c_str(), this->polyOrder);

        dname = frameName + "/class/P";
        df.read_contained_vals (dname.c_str(), this->P);

        dname = frameName + "/class/PP_n";
        unsigned int pp_size = 0;
        df.read_val (dname.c_str(), pp_size);

        this->PP.resize(pp_size);
        for (size_t i = 0; i<pp_size; ++i) {
            std::stringstream ss;
            ss << frameName + "/class/PP";
            ss.width(3);
            ss.fill('0');
            ss << i;
            df.read_contained_vals (ss.str().c_str(), this->PP[i]);
        }

        dname = frameName + "/class/pp_idx";
        df.read_val (dname.c_str(), this->pp_idx);

        // Freehand-drawn regions
        dname = frameName + "/class/FL";
        df.read_contained_vals (dname.c_str(), this->FL);
        dname = frameName + "/class/FLE_n";
        unsigned int fle_size = 0;
        df.read_val (dname.c_str(), fle_size);
        this->FLE.resize(fle_size);
        for (size_t i = 0; i<fle_size; ++i) {
            std::stringstream ss;
            ss << frameName + "/class/FLE";
            ss.width(3);
            ss.fill('0');
            ss << i;
            //std::cout << "Reading into FLE["<<i<<"] with name " << ss.str() << "\n";
            df.read_contained_vals (ss.str().c_str(), this->FLE[i]);
            //std::cout << "FLE[i] now has size " << this->FLE[i].size() << "\n";
        }

        dname = frameName + "/class/nBinsTarg";
        df.read_val (dname.c_str(), this->nBinsTarg);
        this->setBins (this->nBinsTarg);
        dname = frameName + "/class/binA";
        df.read_val (dname.c_str(), this->binA);
        dname = frameName + "/class/binB";
        df.read_val (dname.c_str(), this->binB);
        dname = frameName + "/class/flags";
        df.read_val (dname.c_str(), this->flags);
        dname = frameName + "/class/filename";
        df.read_string (dname.c_str(), this->filename);

        // Don't read back from H5 file what is specified in json; too confusing.
        //dname = frameName + "/class/layer_x";
        //df.read_val (dname.c_str(), this->layer_x);
        //dname = frameName + "/class/thickness";
        //df.read_val (dname.c_str(), this->thickness);
        //dname = frameName + "/class/pixels_per_mm";
        //df.read_val (dname.c_str(), this->pixels_per_mm);
    }

    //! Write the data out to an HdfData file \a df.
    void write (morph::HdfData& df)
    {
        // Update box means. not const
        this->computeBoxMeans();

        // And any freehand regions. not const
        this->computeFreehandMeans();

        std::string frameName = this->getFrameName();

        // Write out essential information to re-load state of the application and the
        // user's work saving points etc.
        std::string dname = frameName + "/class/polyOrder";
        df.add_val (dname.c_str(), this->polyOrder);
        dname = frameName + "/class/P";
        df.add_contained_vals (dname.c_str(), this->P);

        dname = frameName + "/class/PP_n";
        unsigned int pp_size = this->PP.size();
        df.add_val (dname.c_str(), pp_size);
        for (size_t i = 0; i<pp_size; ++i) {
            std::stringstream ss;
            ss << frameName + "/class/PP";
            ss.width(3);
            ss.fill('0');
            ss << i;
            df.add_contained_vals (ss.str().c_str(), this->PP[i]);
        }

        // Freehand drawn regions
        dname = frameName + "/class/FL";
        df.add_contained_vals (dname.c_str(), this->FL);
        dname = frameName + "/class/FLE_n";
        unsigned int fle_size = this->FLE.size();
        df.add_val (dname.c_str(), fle_size);
        for (size_t i = 0; i<fle_size; ++i) {
            std::stringstream ss;
            ss << frameName + "/class/FLE";
            ss.width(3);
            ss.fill('0');
            ss << i;
            df.add_contained_vals (ss.str().c_str(), this->FLE[i]);
        }

        dname = frameName + "/class/pp_idx";
        df.add_val (dname.c_str(), this->pp_idx);
        dname = frameName + "/class/nBinsTarg";
        df.add_val (dname.c_str(), this->nBinsTarg);
        dname = frameName + "/class/binA";
        df.add_val (dname.c_str(), this->binA);
        dname = frameName + "/class/binB";
        df.add_val (dname.c_str(), this->binB);
        dname = frameName + "/class/flags";
        df.add_val (dname.c_str(), this->flags);
        dname = frameName + "/class/filename";
        df.add_string (dname.c_str(), this->filename);
        // The cv::Mat image data could be saved in the h5 file, though at the cost of
        // quite a bit of storage:
        // dname = frameName + "/class/frame";
        // df.add_contained_vals (dname.c_str(), this->frame);
        dname = frameName + "/class/layer_x";
        df.add_val (dname.c_str(), this->layer_x);
        dname = frameName + "/class/thickness";
        df.add_val (dname.c_str(), this->thickness);
        dname = frameName + "/class/pixels_per_mm";
        df.add_val (dname.c_str(), this->pixels_per_mm);
        dname = frameName + "/class/idx";
        df.add_val (dname.c_str(), this->idx);
        // Write the BG blurring parameters
        dname = frameName + "/class/bg_blur_screen_proportion";
        df.add_val (dname.c_str(), this->bgBlurScreenProportion);
        dname = frameName + "/class/bg_blur_subtraction_offset";
        df.add_val (dname.c_str(), this->bgBlurSubtractionOffset);

        /*
         * The rest of the methods write out data that WON'T be read by the
         * FrameData::read method (these would all be re-computed before being
         * re-written in a later run of the program).
         */
        for (size_t bi = 0; bi < this->boxes_raw.size(); ++bi) {
            dname = frameName + "/box" + std::to_string(bi);
            df.add_contained_vals (dname.c_str(), this->boxes_raw[bi]);
        }
        dname = frameName + "/nboxes";
        df.add_val (dname.c_str(), static_cast<unsigned int>(this->boxes_raw.size()));

        dname = frameName + "/means";
        df.add_contained_vals (dname.c_str(), this->means);

        // Autoscale means and save a copy
        dname = frameName + "/means_autoscaled";
        // this->means is vector<double>
        std::vector<double> means_autoscaled = morph::MathAlgo::autoscale (this->means, 0.0, 1.0);
        df.add_contained_vals (dname.c_str(), means_autoscaled);

        // Freehand drawn regions - results
        for (size_t ri = 0; ri < this->FL_raw.size(); ++ri) {
            dname = frameName + "/freehand" + std::to_string(ri);
            //std::cout << "Writing out: " << dname << std::endl;
            df.add_contained_vals (dname.c_str(), this->FL_raw[ri]);
        }
        dname = frameName + "/nfreehand";
        df.add_val (dname.c_str(), static_cast<unsigned int>(this->FL_raw.size()));

        dname = frameName + "/freehand_means";
        df.add_contained_vals (dname.c_str(), this->FL_means);

        // Add the centroid of the freehand regions (in the y-z or 'in-slice' plane)
        for (size_t i = 0; i<fle_size; ++i) {
            cv::Point cntroid = morph::MathAlgo::centroid (this->FLE[i]);
            std::stringstream cntss;
            cntss << frameName + "/freehand" << std::to_string(i) << "_centroid";

            // Offset and scale cntroid suitably (from screen pixels to mm in the slice
            // plane), before saving
            std::cout << "centroid in screen pix: " << cntroid << std::endl;
            cv::Point2d coff = this->offsetPoint (cntroid);
            std::cout << "centroid in scaled pix: " << coff << std::endl;
            cv::Point2d coffrot = this->rotate (this->theta, coff);
            std::cout << "centroid in scaled pix, rotated: " << coffrot << std::endl;

            df.add_contained_vals (cntss.str().c_str(), coffrot);
        }

        // Need to get from fitted to y and z. Note that fitted is in (integer) pixels...
        // vector<cv::Point> fitted;
        //
        // Make up the boxes. A box (in 3d space) can be a vector of 12 floats. Thus
        // we should be able to write a vector of boxes as a vector<vector<float>>
        // These are "surface_boxes" because they're the box thats in the plane of the
        // cortical sheet (roughly xy) rather than the box in the slice plane (yz).
        std::vector<std::array<float,12>> surface_boxes;
#if 0
        std::vector<std::array<float,12>> smooth_boxes; // smoothed surface
#endif
        std::vector<std::array<float,3>> surface_box_centroids;
        std::array<float, 12> sbox;
        //std::cout << "Surface boxes extend from " << layer_x << " to " << (layer_x + thickness) << std::endl;
        for (int i = 1; i < this->nFit; ++i) {
            // c1 x,y,z
            sbox[0] = this->layer_x;                 // x
            sbox[1] = this->fitted_rotated[i-1].x;  // y
            sbox[2] = this->fitted_rotated[i-1].y; // z
            // c2 x,y,z
            sbox[3] = this->layer_x;               // x
            sbox[4] = this->fitted_rotated[i].x;  // y
            sbox[5] = this->fitted_rotated[i].y; // z
            // c3 x,y,z
            sbox[6] = this->layer_x+this->thickness; // x
            sbox[7] = this->fitted_rotated[i].x;     // y
            sbox[8] = this->fitted_rotated[i].y;    // z
            // c4 x,y,z
            sbox[9] = this->layer_x+this->thickness; // x
            sbox[10] = this->fitted_rotated[i-1].x;  // y
            sbox[11] = this->fitted_rotated[i-1].y; // z

            std::array<float, 3> sbox_centroid = morph::MathAlgo::centroid3D (sbox);//<float>
            surface_boxes.push_back (sbox);
            surface_box_centroids.push_back (sbox_centroid);
        }
#if 0
        for (int i = 1; i < this->nFit; ++i) {
            // c1 x,y,z
            sbox[0] = this->layer_x;                 // x
            sbox[1] = this->fitted_rotated[i-1].x;  // y
            sbox[2] = this->fitted_rotated[i-1].y; // z
            // c2 x,y,z
            sbox[3] = this->layer_x;               // x
            sbox[4] = this->fitted_rotated[i].x;  // y
            sbox[5] = this->fitted_rotated[i].y; // z
            // c3 x,y,z
            sbox[6] = this->layer_x+this->thickness; // x
            sbox[7] = this->fitted_rotated[i].x;     // y
            sbox[8] = this->fitted_rotated[i].y;    // z
            // c4 x,y,z
            sbox[9] = this->layer_x+this->thickness; // x
            sbox[10] = this->fitted_rotated[i-1].x;  // y
            sbox[11] = this->fitted_rotated[i-1].y; // z

            smooth_boxes.push_back (sbox);
        }
#endif
        dname = frameName + "/fitted";
        df.add_contained_vals (dname.c_str(), this->fitted);

        dname = frameName + "/fitted_offset";
        df.add_contained_vals (dname.c_str(), this->fitted_offset);

        dname = frameName + "/fitted_rotated";
        df.add_contained_vals (dname.c_str(), this->fitted_rotated);

        // sboxes are 'surface boxes' - they lay in the plan of the cortical surface
        // and are not to be confused with the yellow boxes drawn in the UI in the y-z
        // plane.
        dname = frameName + "/sboxes";
        df.add_contained_vals (dname.c_str(), surface_boxes);

        dname = frameName + "/sbox_centers";
        df.add_contained_vals (dname.c_str(), surface_box_centroids);

        // From surface_box_centroids, can compute linear distance along curve. Could
        // be useful for making naive maps that unroll the cortex in one dimension.
        dname = frameName + "/sbox_linear_distance";
        float total_linear_distance = 0.0f;
        std::vector<float> linear_distances (this->nBins, 0.0f);
        for (int i=1; i<this->nBins; ++i) {
            // Compute distance from Previous to current
            float d = morph::MathAlgo::distance<float> (surface_box_centroids[i-1],
                                                        surface_box_centroids[i]);
            total_linear_distance += d;
            linear_distances[i] = total_linear_distance;
        }
        // Now offset the linear distances so that the middle is 0.
        float halftotal = total_linear_distance / 2.0f;
        for (int i=0; i<this->nBins; ++i) {
            linear_distances[i] -= halftotal;
        }
        df.add_contained_vals (dname.c_str(), linear_distances);

        std::cout << "write() completed for one frame." << std::endl;
    }

    //! Mirror the image and mark in the flags that it was mirrored
    void mirror()
    {
        this->mirror_image_only();
        this->flags.flip (Mirrored);
    }
    //! Carry out the actual mirroring operation on its own, leaving flags unchanged
    void mirror_image_only()
    {
        cv::Mat mirrored (this->frame.rows, this->frame.cols, this->frame.type());
        cv::flip (this->frame, mirrored, 1);
        this->frame = mirrored;
    }

    //! Flip the image & mark as such in flags
    void flip()
    {
        this->flip_image_only();
        this->flags.flip (Flipped);
    }
    //! Flip the image without marking as flipped in flags.
    void flip_image_only()
    {
        cv::Mat flipped (this->frame.rows, this->frame.cols, this->frame.type());
        cv::flip (this->frame, flipped, 1);
        this->frame = flipped;
    }

    //! Recompute the fit
    void updateFit()
    {
        if (this->ct == InputMode::Poly) {
            this->updateFitPoly();
        } else if (this->ct == InputMode::Bezier) {
            this->updateFitBezier();
        } else if (this->ct == InputMode::Freehand) {
            // What to do? Find all the pixels inside?
        } else {
            return;
        }

        if (this->ct == InputMode::Poly || this->ct == InputMode::Bezier) {
            // Scale
            this->offsetScaleFit();
            // Rotate
            this->rotateFitOptimally();
            std::cout << "At end of updateFit(void). binA/binB: " << binA << "," << binB << std::endl;
        }
    }

    //! Re-compute the boxes from the curve (taking ints)
    void refreshBoxes (const int lenA, const int lenB) { this->refreshBoxes ((double)lenA, (double)lenB); }

    //! Re-compute the boxes from the curve (double version)
    void refreshBoxes (const double lenA, const double lenB)
    {

        // Don't refresh boxes for Freehand mode
        if (this->ct == InputMode::Freehand) {
            // Or perhaps hide stuff? Delete points in P? or leave the points in P,
            // and have a separate store of the points in a freehand loop.
            return;
        }

        if (this->ct == InputMode::Poly) {
            this->pointsInner = PolyFit::rotate (PolyFit::tracePolyOrth (this->pf, this->minX, this->maxX,
                                                                         this->nFit, lenA),
                                                 this->theta);
            this->pointsOuter = PolyFit::rotate (PolyFit::tracePolyOrth (this->pf, this->minX, this->maxX,
                                                                         this->nFit, lenB),
                                                 this->theta);
        } else {
            for (int i=0; i<this->nFit; i++) {
                cv::Point2d normLenA = this->normals[i]*lenA;
                cv::Point2d normLenB = this->normals[i]*lenB;
                this->pointsInner[i] = this->fitted[i] + cv::Point2i((int)normLenA.x, (int)normLenA.y);
                this->pointsOuter[i] = this->fitted[i] + cv::Point2i((int)normLenB.x, (int)normLenB.y);
            }
        }
        // Make the boxes from pointsInner and pointsOuter
        this->boxes.resize (this->nBins);
        for (int i=0; i<this->nBins; i++) {
            std::vector<cv::Point> pts(4);
            pts[0] = this->pointsInner[i];
            pts[1] = this->pointsInner[i+1];
            pts[2] = this->pointsOuter[i+1];
            pts[3] = this->pointsOuter[i];
            this->boxes[i] = pts;
        }
    }

    //! Toggle between curve fitting, freehand loop drawing or alignment mark (landmark) input.
    void toggleInputMode()
    {
        // Note: InputMode::Poly is ignored now
        if (this->ct == InputMode::Landmark) {
            this->ct = InputMode::Bezier;
        } else if (this->ct == InputMode::Bezier) {
            this->ct = InputMode::Freehand;
        } else if (this->ct == InputMode::Freehand) {
            this->ct = InputMode::Landmark;
        } else {
            // Shouldn't get here...
            this->ct = InputMode::Bezier;
        }
    }

    // Toggle controls

    void toggleShowBoxes() { this->flags[ShowBoxes] = this->flags.test(ShowBoxes) ? false : true; }
    void setShowBoxes (bool t) { this->flags[ShowBoxes] = t; }

    void toggleShowFits() { this->flags[ShowFits] = this->flags.test(ShowFits) ? false : true; }
    void setShowFits (bool t) { this->flags[ShowFits] = t; }

    void toggleShowUsers() { this->flags[ShowUsers] = this->flags.test(ShowUsers) ? false : true; }
    void setShowUsers (bool t) { this->flags[ShowUsers] = t; }

    void toggleShowCtrls() { this->flags[ShowCtrls] = this->flags.test(ShowCtrls) ? false : true; }
    void setShowCtrls (bool t) { this->flags[ShowCtrls] = t; }

private:

    //! Blur the image and use the blurred result to estimate the overall background
    //! luminance, which may vary due to the arrangement of lighting when the brain
    //! slices were imaged.
    void estimateBackgroundLuminance()
    {
        cv::Mat blurred = cv::Mat::zeros (this->frame.rows, this->frame.cols, CV_8UC3);
        cv::Size ksz;
        ksz.width = 21;
        ksz.height = 21;
        cv::GaussianBlur (this->frame, blurred, ksz, 2.0);
    }

    //! Find a value for each pixel of image \a frame within the box defined by \a pp
    //! and return this in a vector of floats, without conversion
    std::vector<float> getBoxedPixelVals (const std::vector<cv::Point> pp)
    {
        cv::Point pts[4] = {pp[0],pp[1],pp[2],pp[3]};

        // Convert frame_bgoff to uchar, multiplying by 255 on the way:
        cv::Mat frame_bgoffU;
        this->frame_bgoff.convertTo (frame_bgoffU, CV_8UC3, 255.0);

        cv::Mat mask = cv::Mat::zeros(frame_bgoffU.rows, frame_bgoffU.cols, CV_8UC3);
        cv::fillConvexPoly (mask, pts, 4, cv::Scalar(255,255,255));
        cv::Mat maskGray;
        cv::cvtColor (mask, maskGray, cv::COLOR_BGR2GRAY);
        std::vector<cv::Point2i> maskpositives;
        cv::findNonZero (maskGray, maskpositives);
        cv::Mat result, resultGray;
        frame_bgoffU.copyTo (result, mask);
        cv::cvtColor (result, resultGray, cv::COLOR_BGR2GRAY);
        std::vector<cv::Point2i> positives;
        cv::findNonZero (resultGray, positives);
        // As long as we make boxedPixelVals the size of the non-zeros in mask, the
        // numbers will work out!
        std::vector<float> boxedPixelVals (maskpositives.size());
        for (size_t j=0; j<maskpositives.size(); j++) {
            cv::Scalar pixel = resultGray.at<uchar>(maskpositives[j]);
            boxedPixelVals[j] = (float)pixel.val[0];
        }
        return boxedPixelVals;
    }

    //! In a box, obtain colour values as BGR float triplets
    std::vector<std::array<float, 3>> getBoxedPixelColour (const std::vector<cv::Point> pp)
    {
        cv::Point pts[4] = {pp[0],pp[1],pp[2],pp[3]};
        cv::Mat mask = cv::Mat::zeros(frame.rows, frame.cols, CV_8UC3);
        cv::fillConvexPoly (mask, pts, 4, cv::Scalar(255,255,255));
        //
        cv::Mat maskGray;
        cv::cvtColor (mask, maskGray, cv::COLOR_BGR2GRAY);
        std::vector<cv::Point2i> maskpositives;
        cv::findNonZero (maskGray, maskpositives);

        cv::Mat result, resultFloat, resultGray;
        // Fixme: Can probably copy frameF to result, which is already in floating point format.
        frame.copyTo (result, mask); // Note: 'result' will be in BGR format
        cv::cvtColor (result, resultGray, cv::COLOR_BGR2GRAY);
        result.convertTo (resultFloat, CV_32FC3);
        std::vector<cv::Point2i> positives;
        cv::findNonZero (resultGray/*Float*/, positives); // only for CV_8UC1 :( I want 'findNonBlack' on result
        std::vector<std::array<float, 3>> boxedPixelVals (maskpositives.size());
        for (size_t j=0; j<maskpositives.size(); j++) {
            cv::Vec3f pixel = resultFloat.at<cv::Vec3f>(maskpositives[j]);
            // NB: This assumes image is in BGR format and we return in BGR format.
            boxedPixelVals[j][0] = pixel.val[0];
            boxedPixelVals[j][1] = pixel.val[1];
            boxedPixelVals[j][2] = pixel.val[2];
        }

        return boxedPixelVals;
    }

    //! Get the raw pixel values from the region
    std::vector<float> getRegionPixelVals (const std::vector<cv::Point>& region)
    {
        throw std::runtime_error ("FIXME - get same data as getBoxedPixelVals. Also display value");
        cv::Mat mask = cv::Mat::zeros(this->frame.rows, this->frame.cols, CV_8UC3);
        std::cout << "mask rows: " << mask.rows << ", mask cols: " << mask.cols << std::endl;
        for (auto px : region) {
            int _col = px.x;
            int _row = px.y;
            mask.at<cv::Scalar>(_col, _row) = cv::Scalar(255,255,255);
        }

        // From here, same as getBoxedPixelVals
        cv::Mat result, resultGray;
        this->frame.copyTo (result, mask);
        cv::cvtColor (result, resultGray, cv::COLOR_BGR2GRAY);
        std::vector<cv::Point2i> positives;
        cv::findNonZero (resultGray, positives);
        std::vector<float> regionPixelVals (positives.size());
        for (size_t j=0; j<positives.size(); j++) {
            cv::Scalar pixel = resultGray.at<uchar>(positives[j]);
            regionPixelVals[j] = (float)pixel.val[0];
        }
        return regionPixelVals;
    }

    //! For each freehand drawn loop, compute the mean luminance within the loop,
    //! storing in this->FL_means
    void computeFreehandMeans()
    {
        // Loop through FLE. For each set of points, output the points as a list and
        // also compute the mean.
        this->FL_means.resize (this->FLE.size());
        this->FL_raw.resize (this->FLE.size());
        this->FL_raw_bgr.resize (this->FLE.size());
        for (size_t i=0; i<this->FLE.size(); i++) {
            // region is FLE[i]
            this->FL_means[i] = 0.0;

            if (this->cmodel == ColourModel::AllenDevMouse) {
                throw std::runtime_error ("AllenDevMouse ColourModel is not implemented for freehand regions right now");

            } else { // Default is ColourModel::Greyscale
                this->FL_raw[i] = this->getRegionPixelVals (this->FLE[i]);
                for (size_t j=0; j<this->FL_raw[i].size(); j++) {
                    // Signal conversion: x - x_0:
                    float signal = this->FL_raw[i][j] - this->luminosity_cutoff;
                    // m * (x - x_0):
                    signal *= this->luminosity_factor;
                    // Any signal <0 is 0.
                    this->FL_means[i] += (double)(signal > 0.0f ? signal : 0.0f);
                }
                // Divide the signal by the number of pixels in the region
                this->FL_means[i] /= (double)this->FL_raw[i].size();
            }
        }
    }

    //! Compute the mean values for the bins. Not const. But means don't need to be a
    //! member as they're only computed to be written out to file.
    void computeBoxMeans()
    {
        this->boxes_raw.resize (this->boxes.size());
        this->boxes_raw_bgr.resize (this->boxes.size());
        this->means.resize (this->boxes.size());
        for (size_t i=0; i<this->boxes.size(); i++) {

            // Zero the means value
            this->means[i] = 0.0;

            // if luminance value only/greyscale:
            if (this->cmodel == ColourModel::AllenDevMouse) {
                // But we'll have to pass parameters for transforming the colours and
                // determining if they're on the "expressing" axis. This will include
                // a translate matrix, a rotation matrix and ellipse parameters,
                // obtained from the octave script plotcolour.m
                this->boxes_raw_bgr[i] = this->getBoxedPixelColour (this->boxes[i]);

                float ellip_maj_sq = ellip_axes[0] * ellip_axes[0];
                float ellip_min_sq = ellip_axes[1] * ellip_axes[1];
                // std::cout << "box " << i << " has " << this->boxes_raw_bgr[i].size() << " pixels" << std::endl;
                for (size_t j=0; j<this->boxes_raw_bgr[i].size(); j++) {
                    // Perform colour transform here, so that we get a transformed blue value
                    float b = boxes_raw_bgr[i][j][0];
                    float g = boxes_raw_bgr[i][j][1];
                    float r = boxes_raw_bgr[i][j][2];
                    //std::cout << "bgr: " << b << "," << g << "," << r << std::endl;

                    // 1. Translate rgb colour. NB: It's this->colour_trans
                    float b_t = b - colour_trans[0];
                    float g_t = g - colour_trans[1];
                    float r_t = r - colour_trans[2];

                    // 2. Rotate colour. NB: It's this->colour_rot
                    float b_r = colour_rot[0]*b_t + colour_rot[1]*g_t + colour_rot[2]*r_t;
                    float g_r = colour_rot[3]*b_t + colour_rot[4]*g_t + colour_rot[5]*r_t;
                    float r_r = colour_rot[6]*b_t + colour_rot[7]*g_t + colour_rot[8]*r_t;
                    //std::cout << "bgr transformed: " << b_r << "," << g_r << "," << r_r << std::endl;

                    // if (r_r,g_r) lies inside ellipse (given by ellip_axes) then blue_transf equals b_r
                    float blue_transf = 0.0f;
                    float erad = ((g_r*g_r)/ellip_maj_sq) + ((r_r*r_r)/ellip_min_sq);
                    if (erad <= 1.0f) {
                        // Inside ellipse:
                        blue_transf = b_r;
                    }

                    // Now apply signal conversion
                    float signal = blue_transf - this->luminosity_cutoff;
                    // m * (x - x_0):
                    signal *= this->luminosity_factor;
                    // Any signal <0 is 0.
                    // std::cout << "signal value is " << signal << std::endl;
                    this->means[i] += (double)(signal > 0.0f ? signal : 0.0f);
                }
                // Divide the signal by the number of pixels in the box
                this->means[i] /= (double)this->boxes_raw_bgr[i].size();

            } else { // Default is ColourModel::Greyscale

                this->boxes_raw[i] = this->getBoxedPixelVals (this->boxes[i]);

                this->estimateBackgroundLuminance();

                for (size_t j=0; j<this->boxes_raw[i].size(); j++) {
                    // Signal conversion is simple for monochrome pixel values:
                    // x - x_0:
                    //std::cout << "boxes_raw[i="<<i<<"][j="<<j<<"]="<<this->boxes_raw[i][j];
                    float signal = this->boxes_raw[i][j] - this->luminosity_cutoff;
                    //std::cout << ", signal=" << signal << std::endl;
                    // m * (x - x_0):
                    signal *= this->luminosity_factor;
                    //std::cout << "signal*luminosity_factor=" << signal << std::endl;
                    // Any signal <0 is 0.
                    this->means[i] += (double)(signal > 0.0f ? signal : 0.0f);
                }
                // Divide the signal by the number of pixels in the box
                this->means[i] /= (double)this->boxes_raw[i].size();
            }
        }
    }

    //! Update the fit, but don't rotate. Used by rotateFitOptimally()
    void updateFit_norotate()
    {
        if (this->ct == InputMode::Poly) {
            this->updateFitPoly();
        } else {
            this->updateFitBezier();
        }
    }

    //! Update the fit, scale and rotate by \a _theta. Used by rotateFitOptimally()
    void updateFit (double _theta)
    {
        if (this->ct == InputMode::Poly) {
            this->updateFitPoly();
        } else {
            this->updateFitBezier();
        }
        // Scale
        this->offsetScaleFit();
        // Rotate
        this->rotate (_theta);
    }

    //! Recompute the Bezier fit
    void updateFitBezier()
    {
        if (this->PP.empty() && this->P.size() < 2) {
            std::cout << "Too few points to fit" << std::endl;
            return;
        }

        this->bcp.reset();

        // Loop over PP first
        for (auto _P : this->PP) {
            std::vector<std::pair<double,double>> user_points;
            user_points.clear();
            for (auto pt : _P) {
                user_points.push_back (std::make_pair(pt.x, pt.y));
            }

            morph::BezCurve<double> bc;
            if (this->bcp.isNull()) {
                // No previous curves; fit just on user_points
                bc.fit (user_points);
                this->bcp.addCurve (bc);
            } else {
                // Have previous curve, use last control of previous curve to make
                // smooth transition.
                morph::BezCurve<double> last = this->bcp.curves.back();
                bc.fit (user_points, last);
                this->bcp.removeCurve();
                this->bcp.addCurve (last);
                this->bcp.addCurve (bc);

            }
        }

        if (this->P.size()>2) {
            std::vector<std::pair<double,double>> user_points;
            user_points.clear();
            for (auto pt : this->P) {
                user_points.push_back (std::make_pair(pt.x, pt.y));
            }
            morph::BezCurve<double> bc;
            if (this->bcp.isNull()) {
                // No previous curves; fit just on user_points
                bc.fit (user_points);
                this->bcp.addCurve (bc);
            } else {
                morph::BezCurve<double> last = this->bcp.curves.back();
                bc.fit (user_points, last);
                this->bcp.removeCurve();
                this->bcp.addCurve (last);
                this->bcp.addCurve (bc);
            }
        }

        // Update this->fitted
        this->bcp.computePoints (static_cast<unsigned int>(this->nFit));
        std::vector<morph::BezCoord<double>> coords = this->bcp.getPoints();
        std::vector<morph::BezCoord<double>> tans = this->bcp.getTangents();
        std::vector<morph::BezCoord<double>> norms = this->bcp.getNormals();
        //Point2d fitsum;
        for (int i = 0; i < this->nFit; ++i) {
            this->fitted[i] = cv::Point(coords[i].x(),coords[i].y());
            //fitsum += Point2d(coords[i].x(),coords[i].y());
            this->tangents[i] = cv::Point2d(tans[i].x(),tans[i].y());
            this->normals[i] = cv::Point2d(norms[i].x(),norms[i].y());
        }
    }

    void offsetScaleFit()
    {
        cv::Point2d fitsum;
        for (int i = 0; i < this->nFit; ++i) {
            fitsum += cv::Point2d(this->fitted[i].x, this->fitted[i].y);
        }
        this->fit_centroid = fitsum/this->nFit;
        std::cout << "Fit centroid: " << this->fit_centroid << std::endl;
        // Apply offset and scale
        for (int i = 0; i < this->nFit; ++i) {
            cv::Point2d fd(this->fitted[i]);
            this->fitted_offset[i] = (fd - this->fit_centroid) / this->pixels_per_mm;
        }
    }

    //! Move the point \a pt in 'screen pixels', offsetting by the centroid of the
    //! fitted Bezier curve (if applicable) and scaled by this->pixels_per_mm.
    cv::Point2d offsetPoint (const cv::Point& pt)
    {
        // Apply offset and scale
        cv::Point2d fd;
        fd.x = (double)pt.x;
        fd.y = (double)pt.y;
        if (this->pixels_per_mm == 0) {
            std::cerr << "WARNING: pixels_per_mm is 0, will get divide by zero..." << std::endl;
        }
        cv::Point2d rtn = (fd - this->fit_centroid) / this->pixels_per_mm;
        return rtn;
    }

    //! Compute a sum of squared distances between the points in this fit and the
    //! points in the previous fit
    double computeSosWithPrev (double _theta)
    {
        if (this->previous < 0) {
            return std::numeric_limits<double>::max();
        }

        if ((*this->parentStack)[this->previous].getNBins() != this->nBins) {
            // Number of bins has to be same
            return std::numeric_limits<double>::max();
        }

        this->rotate (_theta);

        double sos = 0.0;
        for (int i = 0; i < this->nFit; ++i) {
            // Distance^2 from this->fitted_rotated[i] to this->previous->fitted_rotated[i]
            double xdiff = this->fitted_rotated[i].x - (*this->parentStack)[this->previous].fitted_rotated[i].x;
            double ydiff = this->fitted_rotated[i].y - (*this->parentStack)[this->previous].fitted_rotated[i].y;
            double d_ = (xdiff*xdiff + ydiff*ydiff);
            //std::cout << "sos += " << d_ << ", ";
            sos += d_;
        }
        //std::cout << "\nFor rotation angle " << _theta << " returning sos=" << sos << std::endl;
        return sos;
    }

    //! Rotate the points in this->fitted_offset by theta about the origin and store
    //! the result in this->fitted_rotated
    void rotate (double _theta)
    {
        if (_theta == 0.0) {
            for (int i = 0; i < this->nFit; ++i) {
                this->fitted_rotated[i] = this->fitted_offset[i];
            }
            return;
        }

        for (int i = 0; i < this->nFit; ++i) {
            double xi = this->fitted_offset[i].x;
            double yi = this->fitted_offset[i].y;
            double sin_theta = sin (_theta);
            double cos_theta = cos (_theta);
            this->fitted_rotated[i].x = xi * cos_theta - yi * sin_theta;
            this->fitted_rotated[i].y = xi * sin_theta + yi * cos_theta;
        }
    }

    //! Rotate the point \a pt by theta about the origin, returning rotated point
    cv::Point2d rotate (double _theta, cv::Point2d pt)
    {
        cv::Point2d rtn = pt;

        if (_theta == 0.0) {
            return rtn;
        }

        double xi = pt.x;
        double yi = pt.y;
        double sin_theta = sin (_theta);
        double cos_theta = cos (_theta);
        rtn.x = xi * cos_theta - yi * sin_theta;
        rtn.y = xi * sin_theta + yi * cos_theta;

        return rtn;
    }

    //! Rotate the fit until we get the best one.
    void rotateFitOptimally()
    {
        // If there's no previous frame, then fitted_rotated should be same as fitted_offset
        if (this->previous < 0) {
            std::cout << "No previous frame, so just copying fitted_offset to fitted_rotated..." << std::endl;
            for (int i = 0; i < this->nFit; ++i) {
                this->fitted_rotated[i] = this->fitted_offset[i];
            }
            return;
        }

        std::cout << "rotateFitOptimally: DO have previous frame" << std::endl;

        // Now check if the previous frame has different number of bins
        int nBinsSave = this->nBins;
        int nBinsTmp = (*this->parentStack)[this->previous].getNBins();
        std::cout << "  nBinsSave(this->nBins) = "<< nBinsSave << std::endl;
        std::cout << "  nBinsTmp(this->previous->nBins) = "<< nBinsTmp << std::endl;
        if (nBinsTmp != nBinsSave) {
            // This temporarily re-computes THIS frame's fit with nBinsTmp
            std::cout << "set bins to nBinsTmp = " << nBinsTmp << std::endl;
            this->setBins (nBinsTmp);
            this->updateFit_norotate();
        }

        // Now we have a fit which has same number of bins as the previous frame,
        // this means we can compute SOS objective function
        double thet1 = 0.0;
        double thet2 = 0.5;
        morph::NM_Simplex<double> simp (thet1, thet2);
        // Set a termination threshold for the SD of the vertices of the simplex
        simp.termination_threshold = 2.0 * std::numeric_limits<double>::epsilon();
        // Set a 10000 operation limit, in case the above threshold can't be reached
        simp.too_many_operations = 1000;

        while (simp.state != morph::NM_Simplex_State::ReadyToStop) {

            if (simp.state == morph::NM_Simplex_State::NeedToComputeThenOrder) {
                // 1. apply objective to each vertex
                for (unsigned int i = 0; i <= simp.n; ++i) {
                    simp.values[i] = this->computeSosWithPrev (simp.vertices[i][0]);
                }
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToOrder) {
                simp.order();

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeReflection) {
                double val = this->computeSosWithPrev (simp.xr[0]);
                simp.apply_reflection (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeExpansion) {
                double val = this->computeSosWithPrev (simp.xe[0]);
                simp.apply_expansion (val);

            } else if (simp.state == morph::NM_Simplex_State::NeedToComputeContraction) {
                double val = this->computeSosWithPrev (simp.xc[0]);
                simp.apply_contraction (val);
            }
        }
        std::vector<double> vP = simp.best_vertex();
        double min_sos = simp.best_value();

        std::cout << "Best sos value: " << min_sos << " and best theta: " << vP[0] << std::endl;

        if (nBinsTmp != nBinsSave) {
            // Need to reset bins and update the fit again, but this time rotating by vP[0]
            std::cout << "reset bins to nBinsSave = " << nBinsSave << std::endl;
            this->setBins (nBinsSave);
            std::cout << "Update fit with rotation " << vP[0] << std::endl;
            this->updateFit (vP[0]);
        } else {
            // There was no need to change bin; rotate by the best theta, vP[0]
            std::cout << "Update unchanged fit with rotation " << vP[0] << std::endl;
            this->rotate (vP[0]);
        }
    }

    //! Recompute the polynomial fit
    void updateFitPoly()
    {
        this->axiscoefs = PolyFit::polyfit (this->P, 1);
        this->axis = PolyFit::tracePoly (this->axiscoefs, 0, this->frame.cols, 2);
        this->theta = atan (this->axiscoefs[1]);
        std::vector<cv::Point> rotated = PolyFit::rotate (this->P, -this->theta);

        this->maxX = -1e9;
        this->minX = +1e9;
        for (size_t i=0; i<rotated.size(); i++) {
            if (rotated[i].x > this->maxX) { this->maxX = rotated[i].x; }
            if (rotated[i].x < this->minX) { this->minX = rotated[i].x; }
        }
        this->pf = PolyFit::polyfit (rotated, this->polyOrder);

        this->fitted = PolyFit::rotate (PolyFit::tracePoly (this->pf, this->minX, this->maxX, this->nFit),
                                        this->theta);
    }

    //! Common code to generate the frame name
    std::string getFrameName() const
    {
        std::stringstream ss;
        ss << "/Frame";
        ss.width(3);
        ss.fill('0');
        ss << (1+this->idx); // Count from 1 in the data file
        return ss.str();
    }

}; // FrameData
