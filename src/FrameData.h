#pragma once

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using cv::Mat;
using cv::Point;
using cv::Point2d;
using cv::Point2i;
using cv::Scalar;
#include <stdexcept>
using std::runtime_error;
using std::exception;
#include <vector>
using std::vector;
#include <array>
using std::array;
#include <string>
using std::string;
using std::to_string;
#include <sstream>
using std::stringstream;
#include <iostream>
using std::cout;
using std::endl;
using std::flush;
#include <utility>
using std::make_pair;
using std::pair;
#include <bitset>
using std::bitset;
#include <limits>
using std::numeric_limits;
#include "tools.h"
#include <morph/BezCurvePath.h>
using morph::BezCurvePath;
#include <morph/BezCurve.h>
using morph::BezCurve;
#include <morph/BezCoord.h>
using morph::BezCoord;
#include <morph/HdfData.h>
using morph::HdfData;
#include <morph/MathAlgo.h>
using morph::MathAlgo;
#include <morph/NM_Simplex.h>
using morph::NM_Simplex;
using morph::NM_Simplex_State;
#include <morph/MathConst.h>

enum class CurveType {
    Poly,  // Variable order polynomial
    Bezier // Cubic Bezier
};

// What sort of colour model is in use?
enum class ColourModel {
    Greyscale,
    AllenDevMouse
};

enum Flag {
    ShowBoxes, // Show the yellow boxes?
    ShowUsers, // Show the user points?
    ShowCtrls, // Show the ctrl points of the fits?
    ShowFits,  // Show the fits?
    Mirrored,  // Is the image mirrored L-R?
    Flipped    // Is the image flipped U-D?
};

/*!
 * A class to hold a cortical section image, user-supplied cortex edge points and the
 * resulting fit (either polynomial or Bezier curved).
 */
class FrameData
{
    //! Private attributes
private:
    //! The 'previous' frame in the stack of frames (index into)
    int previous = -1;
    vector<FrameData>* parentStack;

    //! Number of bins to create for the fit (one less than nFit)
    int nBins;
    //! Number of points to create in the fit
    int nFit;

    //! Public attributes
public:
    //! What curve type?
    CurveType ct = CurveType::Bezier;

    //! What kind of colour model is in use?
    ColourModel cmodel = ColourModel::Greyscale;

    //! Polynomial fit specific attributes
    //@{
    //! Order of the polynomial fit
    int polyOrder;
    //! A rotation is applied to the points before calling polyfit, after which fitted
    //! points are rotated back again.
    double theta;
    //! Min and max value used as input to polyfit
    double minX, maxX;
    //! The parameters of a polyfit, as returned by polyfit()
    vector<double> pf;
    //@}

    //! Bezier curve attributes
    //@{
    //! A Bezier curve path to fit the cortex.
    BezCurvePath<double> bcp;
    //@}

    //! Attributes which pertain either to polynomial or Bezier curves
    //@{
    //! The vector of user-supplied points from which to make a curve fit
    vector<Point> P;
    //! A vector of vectors of points for multi-section Bezier curves
    vector<vector<Point>> PP;
    //! Index into PP
    int pp_idx = 0;
    //! The means computed for the boxes.
    vector<double> means;
    //! The raw values for each box as a vector of doubles for each box.
    vector<vector<double>> boxes_raw;
    //! Target number of bins; used by bins slider
    int nBinsTarg;
    //! The bin lengths, set with a slider.
    int binA = 0;
    int binB = 100;
    //! A set of points created from the fit
    vector<Point> fitted;
    //! The centroid of fitted.
    Point2d fit_centroid;
    //! This holds offset and scaled fitted points: (fitted - fit_centroid) * pixels_per_mm
    vector<Point2d> fitted_offset;
    //! This holds the fitted_offset points after they have been rotated to be in line
    //! with fitted_rotated points in the 'previous' frame. 'in line' means the
    //! smallest sum-of-square distances between the two fitted_rotated sets. Depends
    //! on nBins being the same in each.
    vector<Point2d> fitted_rotated;
    //! For point in fitted, the tangent at that location
    vector<Point2d> tangents;
    //! For point in fitted, the normal at that location
    vector<Point2d> normals;
    //! The axes for the polynomial fit
    vector<Point> axis;
    //! Coefficients for the polynomial fit axis
    vector<double> axiscoefs;
    //! origins for the lines making the box sides (some distance from the curve)
    vector<Point> pointsInner;
    //! endpoints for the lines making the box sides (a greater distance from the curve)
    vector<Point> pointsOuter;
    //! The boxes that are drawn and from which to sample the gene expression
    vector<vector<Point> > boxes;
    //! A bit set containing flags
    bitset<8> flags;
    //! The image data, required when sampling the image in one of the boxes.
    Mat frame;
    //! The frame image filename from which frame was loaded. Stored so it can be
    //! recorded when writing out.
    string filename;
    //! The 'x' position of the brain slice in this frame (coordinates in the plane of
    //! the slice are y/z (units: mm).
    float layer_x = 0.0f;
    //! The thickness of this slice (mm)
    float thickness = 0.05;
    //! The scaling factor
    double pixels_per_mm = 100.0;
    //! The index of the frame
    int idx;
    //@}

public:
    FrameData () {
        throw runtime_error ("Default constructor is not allowed");
    }
    //! Constructor initializes default values
    FrameData (const Mat& fr) {
        // init previous to null.
        this->previous = -1;
        this->parentStack = (vector<FrameData>*)0;
        this->frame = fr.clone();
        this->axiscoefs.resize (2, 0.0);
        this->axis.resize (2);
        // NB: Init these before the next three resize() calls
        this->setBins (100);
        this->polyOrder = 3;
        // Init flags
        this->flags.set (ShowFits);
        this->flags.set (ShowUsers);
        this->flags.set (ShowBoxes);
    };

    //! Set the number of bins and update the size of the various containers
    void setBins (unsigned int num) {
        if (num > 5000) {
            throw runtime_error ("Too many bins...");
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

    //! Getter for nBins
    int getNBins (void) {
        return this->nBins;
    }
    int getNBinsTarg (void) {
        return this->nBinsTarg;
    }

    //! Setter for previous
    void setPrevious (int prev) {
        this->previous = prev;
    }
    void setParentStack (vector<FrameData>* parentSt) {
        this->parentStack = parentSt;
    }

    //! Get information about the fit
    string getFitInfo (void) const {
        stringstream ss;
        if (this->ct == CurveType::Poly) {
            ss << "Poly order: " << this->polyOrder << ", Bins: " << this->nBins;
        } else if (this->ct == CurveType::Bezier) {
            stringstream bb;
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
        } else {
            ss << "unknown";
        }
        return ss.str();
    }

    //! Remove the last user point
    void removeLastPoint (void) {
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
            if (this->ct == CurveType::Bezier && this->pp_idx>0) {
                this->P = this->PP[--this->pp_idx];
                this->PP.pop_back();
                this->P.pop_back();
            }
        }
    }

    //! In Bezier mode, store the current set of user points (P) into PP and clear P.
    void nextCurve (void) {
        if (this->ct == CurveType::Poly) {
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

    //! Compute the mean values for the bins
    void getBoxMeans (void) {
        this->boxes_raw.resize (this->boxes.size());
        this->means.resize (this->boxes.size());
        for (size_t i=0; i<this->boxes.size(); i++) {
            // if luminance value only/greyscale:
            this->cmodel = ColourModel::AllenDevMouse; // HACK
            if (this->cmodel == ColourModel::AllenDevMouse) {
                // But we'll have to pass parameters for transforming the colours and
                // determining if they're on the "expressing" axis. This will include
                // a translate matrix, a rotation matrix and ellipse parameters,
                // obtained from the octave script plotcolour.m
                this->boxes_raw[i] = StaleUtil::getAllenPixelVals (this->frame, this->boxes[i]);
            } else { // Default is ColourModel::Greyscale
                this->boxes_raw[i] = StaleUtil::getBoxedPixelVals (this->frame, this->boxes[i]);
            }

            this->means[i] = 0.0;
            for (size_t j=0; j<this->boxes_raw[i].size(); j++) {
                this->means[i] += this->boxes_raw[i][j];
            }
            this->means[i] /= (double)this->boxes_raw[i].size();
        }
    }

    //! Read important data from file
    void read (HdfData& df) {
        // Note this file assumes idx has been set for the frame.
        string frameName = this->getFrameName();

        string dname = frameName + "/class/polyOrder";

        df.read_val (dname.c_str(), this->polyOrder);

        dname = frameName + "/class/P";
        df.read_contained_vals (dname.c_str(), this->P);

        dname = frameName + "/class/PP_n";
        unsigned int pp_size = 0;
        df.read_val (dname.c_str(), pp_size);

        this->PP.resize(pp_size);
        for (size_t i = 0; i<pp_size; ++i) {
            stringstream ss;
            ss << frameName + "/class/PP";
            ss.width(3);
            ss.fill('0');
            ss << i;
            df.read_contained_vals (ss.str().c_str(), this->PP[i]);
        }

        dname = frameName + "/class/pp_idx";
        df.read_val (dname.c_str(), this->pp_idx);

        dname = frameName + "/class/nBinsTarg";
        df.read_val (dname.c_str(), this->nBinsTarg);
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

    //! Write the data out to an HdfData file @df.
    void write (HdfData& df) const {
        string frameName = this->getFrameName();

        // Write out essential information to re-load state of the application and the
        // user's work saving points etc.
        string dname = frameName + "/class/polyOrder";
        df.add_val (dname.c_str(), this->polyOrder);
        dname = frameName + "/class/P";
        df.add_contained_vals (dname.c_str(), this->P);

        dname = frameName + "/class/PP_n";
        unsigned int pp_size = this->PP.size();
        df.add_val (dname.c_str(), pp_size);
        for (size_t i = 0; i<pp_size; ++i) {
            stringstream ss;
            ss << frameName + "/class/PP";
            ss.width(3);
            ss.fill('0');
            ss << i;
            df.add_contained_vals (ss.str().c_str(), this->PP[i]);
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

        /*
         * The rest of the methods write out data that WON'T be read by the
         * FrameData::read method (these would all be re-computed before being
         * re-written in a later run of the program).
         */
        for (size_t bi = 0; bi < this->boxes_raw.size(); ++bi) {
            dname = frameName + "/box" + to_string(bi);
            df.add_contained_vals (dname.c_str(), this->boxes_raw[bi]);
        }
        dname = frameName + "/nboxes";
        df.add_val (dname.c_str(), static_cast<unsigned int>(this->boxes_raw.size()));

        dname = frameName + "/means";
        df.add_contained_vals (dname.c_str(), this->means);

        // Autoscale means and save a copy
        dname = frameName + "/means_autoscaled";
        vector<double> means_autoscaled = MathAlgo<double>::autoscale (this->means);
        df.add_contained_vals (dname.c_str(), means_autoscaled);

        // Need to get from fitted to y and z. Note that fitted is in (integer) pixels...
        // vector<Point> fitted;
        //
        // Make up the boxes. A box (in 3d space) can be a vector of 12 floats. Thus
        // we should be able to write a vector of boxes as a vector<vector<float>>
        // These are "surface_boxes" because they're the box thats in the plane of the
        // cortical sheet (roughly xy) rather than the box in the slice plane (yz).
        vector<array<float,12>> surface_boxes;
#if 0
        vector<array<float,12>> smooth_boxes; // smoothed surface
#endif
        vector<array<float,3>> surface_box_centroids;
        array<float, 12> sbox;
        cout << "Surface boxes extend from " << layer_x << " to " << (layer_x + thickness) << endl;
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

            array<float, 3> sbox_centroid = MathAlgo<float>::centroid3D (sbox);
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
        vector<float> linear_distances (this->nBins, 0.0f);
        for (int i=1; i<this->nBins; ++i) {
            // Compute distance from Previous to current
            float d = MathAlgo<float>::distance (surface_box_centroids[i-1],
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

        cout << "write() completed." << endl;
    }

    //! Mirror the image and mark in the flags that it was mirrored
    void mirror (void) {
        this->mirror_image_only();
        this->flags.flip (Mirrored);
    }
    //! Carry out the actual mirroring operation on its own, leaving flags unchanged
    void mirror_image_only (void) {
        Mat mirrored (this->frame.rows, this->frame.cols, this->frame.type());
        cv::flip (this->frame, mirrored, 1);
        this->frame = mirrored;
    }

    //! Flip the image & mark as such in flags
    void flip (void) {
        this->flip_image_only();
        this->flags.flip (Flipped);
    }
    //! Flip the image without marking as flipped in flags.
    void flip_image_only (void) {
        Mat flipped (this->frame.rows, this->frame.cols, this->frame.type());
        cv::flip (this->frame, flipped, 1);
        this->frame = flipped;
    }

    //! Recompute the fit
    void updateFit (void) {
        if (this->ct == CurveType::Poly) {
            this->updateFitPoly();
        } else {
            this->updateFitBezier();
        }
        // Scale
        this->offsetScaleFit();
        // Rotate
        this->rotateFitOptimally();

        cout << "At end of updateFit(void). binA/binB: " << binA << "," << binB << endl;
    }

    //! Re-compute the boxes from the curve (taking ints)
    void refreshBoxes (const int lenA, const int lenB) {
        this->refreshBoxes ((double)lenA, (double)lenB);
    }

    //! Re-compute the boxes from the curve (double version)
    void refreshBoxes (const double lenA, const double lenB) {
        if (this->ct == CurveType::Poly) {
            this->pointsInner = PolyFit::rotate (PolyFit::tracePolyOrth (this->pf, this->minX, this->maxX,
                                                                         this->nFit, lenA),
                                                 this->theta);
            this->pointsOuter = PolyFit::rotate (PolyFit::tracePolyOrth (this->pf, this->minX, this->maxX,
                                                                         this->nFit, lenB),
                                                 this->theta);
        } else {
            for (int i=0; i<this->nFit; i++) {
                Point2d normLenA = this->normals[i]*lenA;
                Point2d normLenB = this->normals[i]*lenB;
                this->pointsInner[i] = this->fitted[i] + Point2i((int)normLenA.x, (int)normLenA.y);
                this->pointsOuter[i] = this->fitted[i] + Point2i((int)normLenB.x, (int)normLenB.y);
            }
        }
        // Make the boxes from pointsInner and pointsOuter
        this->boxes.resize (this->nBins);
        for (int i=0; i<this->nBins; i++) {
            vector<Point> pts(4);
            pts[0] = this->pointsInner[i];
            pts[1] = this->pointsInner[i+1];
            pts[2] = this->pointsOuter[i+1];
            pts[3] = this->pointsOuter[i];
            this->boxes[i] = pts;
        }
    }

    //! Toggle between polynomial and Bezier curve fitting
    void toggleCurveType (void) {
        if (this->ct == CurveType::Poly) {
            this->ct = CurveType::Bezier;
        } else {
            this->ct = CurveType::Poly;
        }
    }

    //! Toggle controls
    //@{
    void toggleShowBoxes (void) {
        this->flags[ShowBoxes] = this->flags.test(ShowBoxes) ? false : true;
    }
    void setShowBoxes (bool t) {
        this->flags[ShowBoxes] = t;
    }

    void toggleShowFits (void) {
        this->flags[ShowFits] = this->flags.test(ShowFits) ? false : true;
    }
    void setShowFits (bool t) {
        this->flags[ShowFits] = t;
    }

    void toggleShowUsers (void) {
        this->flags[ShowUsers] = this->flags.test(ShowUsers) ? false : true;
    }
    void setShowUsers (bool t) {
        this->flags[ShowUsers] = t;
    }

    void toggleShowCtrls (void) {
        this->flags[ShowCtrls] = this->flags.test(ShowCtrls) ? false : true;
    }
    void setShowCtrls (bool t) {
        this->flags[ShowCtrls] = t;
    }
    //@}

private:
    /*!
     * Private methods
     */
    //! Update the fit, but don't rotate. Used by rotateFitOptimally()
    void updateFit_norotate (void) {
        if (this->ct == CurveType::Poly) {
            this->updateFitPoly();
        } else {
            this->updateFitBezier();
        }
    }

    //! Update the fit, scale and rotate by @_theta. Used by rotateFitOptimally()
    void updateFit (double _theta) {
        if (this->ct == CurveType::Poly) {
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
    void updateFitBezier (void) {

        if (this->PP.empty() && this->P.size() < 2) {
            cout << "Too few points to fit" << endl;
            return;
        }

        this->bcp.reset();

        // Loop over PP first
        for (auto _P : this->PP) {
            vector<pair<double,double>> user_points;
            user_points.clear();
            for (auto pt : _P) {
                user_points.push_back (make_pair(pt.x, pt.y));
            }

            BezCurve<double> bc;
            if (this->bcp.isNull()) {
                // No previous curves; fit just on user_points
                bc.fit (user_points);
                this->bcp.addCurve (bc);
            } else {
                // Have previous curve, use last control of previous curve to make
                // smooth transition.
                BezCurve<double> last = this->bcp.curves.back();
                bc.fit (user_points, last);
                this->bcp.removeCurve();
                this->bcp.addCurve (last);
                this->bcp.addCurve (bc);

            }
        }

        if (this->P.size()>2) {
            vector<pair<double,double>> user_points;
            user_points.clear();
            for (auto pt : this->P) {
                user_points.push_back (make_pair(pt.x, pt.y));
            }
            BezCurve<double> bc;
            if (this->bcp.isNull()) {
                // No previous curves; fit just on user_points
                bc.fit (user_points);
                this->bcp.addCurve (bc);
            } else {
                BezCurve<double> last = this->bcp.curves.back();
                bc.fit (user_points, last);
                this->bcp.removeCurve();
                this->bcp.addCurve (last);
                this->bcp.addCurve (bc);
            }
        }

        // Update this->fitted
        this->bcp.computePoints (static_cast<unsigned int>(this->nFit));
        vector<BezCoord<double>> coords = this->bcp.getPoints();
        vector<BezCoord<double>> tans = this->bcp.getTangents();
        vector<BezCoord<double>> norms = this->bcp.getNormals();
        //Point2d fitsum;
        for (int i = 0; i < this->nFit; ++i) {
            this->fitted[i] = Point(coords[i].x(),coords[i].y());
            //fitsum += Point2d(coords[i].x(),coords[i].y());
            this->tangents[i] = Point2d(tans[i].x(),tans[i].y());
            this->normals[i] = Point2d(norms[i].x(),norms[i].y());
        }
    }

    void offsetScaleFit (void) {
        Point2d fitsum;
        for (int i = 0; i < this->nFit; ++i) {
            fitsum += Point2d(this->fitted[i].x, this->fitted[i].y);
        }
        this->fit_centroid = fitsum/this->nFit;
        cout << "Fit centroid: " << this->fit_centroid << endl;
        // Apply offset and scale
        for (int i = 0; i < this->nFit; ++i) {
            Point2d fd(this->fitted[i]);
            this->fitted_offset[i] = (fd - this->fit_centroid) / this->pixels_per_mm;
        }
    }

    //! Compute a sum of squared distances between the points in this fit and the
    //! points in the previous fit
    double computeSosWithPrev (double theta_) {

        if (this->previous < 0) {
            return numeric_limits<double>::max();
        }

        if ((*this->parentStack)[this->previous].getNBins() != this->nBins) {
            // Number of bins has to be same
            return numeric_limits<double>::max();
        }

        this->rotate (theta_);

        double sos = 0.0;
        for (int i = 0; i < this->nFit; ++i) {
            // Distance^2 from this->fitted_rotated[i] to this->previous->fitted_rotated[i]
            double xdiff = this->fitted_rotated[i].x - (*this->parentStack)[this->previous].fitted_rotated[i].x;
            double ydiff = this->fitted_rotated[i].y - (*this->parentStack)[this->previous].fitted_rotated[i].y;
            double d_ = (xdiff*xdiff + ydiff*ydiff);
            //cout << "sos += " << d_ << ", ";
            sos += d_;
        }
        //cout << "\nFor rotation angle " << theta_ << " returning sos=" << sos << endl;
        return sos;
    }

    //! Rotate the points in this->fitted_offset by theta about the origin and store
    //! the result in this->fitted_rotated
    void rotate (double theta) {
        if (theta == 0.0) {
            for (int i = 0; i < this->nFit; ++i) {
                this->fitted_rotated[i] = this->fitted_offset[i];
            }
            return;
        }

        for (int i = 0; i < this->nFit; ++i) {
            double xi = this->fitted_offset[i].x;
            double yi = this->fitted_offset[i].y;
            double sin_theta = sin (theta);
            double cos_theta = cos (theta);
            this->fitted_rotated[i].x = xi * cos_theta - yi * sin_theta;
            this->fitted_rotated[i].y = xi * sin_theta + yi * cos_theta;
        }
    }

    //! Rotate the fit until we get the best one.
    void rotateFitOptimally (void) {

        // If there's no previous frame, then fitted_rotated should be same as fitted_offset
        if (this->previous < 0) {
            cout << "No previous frame, so just copying fitted_offset to fitted_rotated..." << endl;
            for (int i = 0; i < this->nFit; ++i) {
                this->fitted_rotated[i] = this->fitted_offset[i];
            }
            return;
        }

        cout << "rotateFitOptimally: DO have previous frame" << endl;

        // Now check if the previous frame has different number of bins
        int nBinsSave = this->nBins;
        int nBinsTmp = (*this->parentStack)[this->previous].getNBins();
        cout << "  nBinsSave(this->nBins) = "<< nBinsSave << endl;
        cout << "  nBinsTmp(this->previous->nBins) = "<< nBinsTmp << endl;
        if (nBinsTmp != nBinsSave) {
            // This temporarily re-computes THIS frame's fit with nBinsTmp
            cout << "set bins to nBinsTmp = " << nBinsTmp << endl;
            this->setBins (nBinsTmp);
            this->updateFit_norotate();
        }

        // Now we have a fit which has same number of bins as the previous frame,
        // this means we can compute SOS objective function
        double thet1 = 0.0;
        double thet2 = 0.5;
        NM_Simplex<double> simp (thet1, thet2);
        // Set a termination threshold for the SD of the vertices of the simplex
        simp.termination_threshold = 2.0 * numeric_limits<double>::epsilon();
        // Set a 10000 operation limit, in case the above threshold can't be reached
        simp.too_many_operations = 1000;

        while (simp.state != NM_Simplex_State::ReadyToStop) {

            if (simp.state == NM_Simplex_State::NeedToComputeThenOrder) {
                // 1. apply objective to each vertex
                for (unsigned int i = 0; i <= simp.n; ++i) {
                    simp.values[i] = this->computeSosWithPrev (simp.vertices[i][0]);
                }
                simp.order();

            } else if (simp.state == NM_Simplex_State::NeedToOrder) {
                simp.order();

            } else if (simp.state == NM_Simplex_State::NeedToComputeReflection) {
                double val = this->computeSosWithPrev (simp.xr[0]);
                simp.apply_reflection (val);

            } else if (simp.state == NM_Simplex_State::NeedToComputeExpansion) {
                double val = this->computeSosWithPrev (simp.xe[0]);
                simp.apply_expansion (val);

            } else if (simp.state == NM_Simplex_State::NeedToComputeContraction) {
                double val = this->computeSosWithPrev (simp.xc[0]);
                simp.apply_contraction (val);
            }
        }
        vector<double> vP = simp.best_vertex();
        double min_sos = simp.best_value();

        cout << "Best sos value: " << min_sos << " and best theta: " << vP[0] << endl;

        if (nBinsTmp != nBinsSave) {
            // Need to reset bins and update the fit again, but this time rotating by vP[0]
            cout << "reset bins to nBinsSave = " << nBinsSave << endl;
            this->setBins (nBinsSave);
            cout << "Update fit with rotation " << vP[0] << endl;
            this->updateFit (vP[0]);
        } else {
            // There was no need to change bin; rotate by the best theta, vP[0]
            cout << "Update unchanged fit with rotation " << vP[0] << endl;
            this->rotate (vP[0]);
        }
    }

    //! Recompute the polynomial fit
    void updateFitPoly (void) {
        this->axiscoefs = PolyFit::polyfit (this->P, 1);
        this->axis = PolyFit::tracePoly (this->axiscoefs, 0, this->frame.cols, 2);
        this->theta = atan (this->axiscoefs[1]);
        vector<Point> rotated = PolyFit::rotate (this->P, -this->theta);

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
    string getFrameName (void) const {
        stringstream ss;
        ss << "/Frame";
        ss.width(3);
        ss.fill('0');
        ss << this->idx;
        return ss.str();
    }

}; // FrameData
