#pragma once

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using cv::Mat;
using cv::Point;
using cv::Point2d;
using cv::Point2i;
using cv::Scalar;
#include <vector>
using std::vector;
#include <array>
using std::array;
#include <string>
using std::string;
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

enum class CurveType {
    Poly,  // Variable order polynomial
    Bezier // Cubic Bezier
};

enum Flag {
    ShowBoxes, // Show the yellow boxes?
    ShowUsers, // Show the user points?
    ShowCtrls, // Show the ctrl points of the fits?
    ShowFits // Show the fits?
};

/*!
 * A class to hold a cortical section image, user-supplied cortex edge points and the
 * resulting fit (either polynomial or Bezier curved).
 */
class FrameData
{
public:
    //! What curve type?
    CurveType ct = CurveType::Bezier;

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
private:
    //! Number of bins to create for the fit (one less than nFit)
    int nBins;
    //! Number of points to create in the fit
    int nFit;
public:
    //! Target number of bins; used by bins slider
    int nBinsTarg;
    //! The bin lengths, set with a slider.
    int binA;
    int binB;
    //! A set of points created from the fit
    vector<Point> fitted;
    //! The centroid of fitted.
    Point2d fit_centroid;
    //! fitted - fit_centroid
    vector<Point2d> fitted_offset;
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

    //! Constructor initializes default values
    FrameData (const Mat& fr) {
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

    void setBins (unsigned int num=1) {
        this->nBins = num;
        this->nBinsTarg = num;
        this->nFit = num + 1;
        this->fitted.resize (this->nFit);
        this->fitted_offset.resize (this->nFit);
        this->pointsInner.resize (this->nFit);
        this->pointsOuter.resize (this->nFit);
        this->tangents.resize (this->nFit);
        this->normals.resize (this->nFit);
    }

    void incBins (unsigned int num=1) {
        this->setBins (this->nBins + num);
    }

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
        cout << "nextCurve; pp_idx is now " << pp_idx << endl;
    }

    void getBoxMeans (void) {
        this->means.resize (this->boxes.size());
        for (size_t i=0; i<this->boxes.size(); i++) {
            vector<double> boxVals = StaleUtil::getBoxedPixelVals (this->frame, this->boxes[i]);
            this->means[i] = 0.0;
            for (size_t j=0; j<boxVals.size(); j++) {
                this->means[i] += boxVals[j];
            }
            this->means[i] /= (double)boxVals.size();
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
        df.add_val (dname.c_str(), this->pp_idx);

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

        dname = frameName + "/class/layer_x";
        df.read_val (dname.c_str(), this->layer_x);
        dname = frameName + "/class/thickness";
        df.read_val (dname.c_str(), this->thickness);
        dname = frameName + "/class/pixels_per_mm";
        df.read_val (dname.c_str(), this->pixels_per_mm);
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

        dname = frameName + "/means";
        df.add_contained_vals (dname.c_str(), this->means);

        // Need to get from fitted to y and z. Note that fitted is in (integer) pixels...
        // vector<Point> fitted;
        //
        // Make up the boxes. A box (in 3d space) can be a vector of 12 floats. Thus
        // we should be able to write a vector of boxes as a vector<vector<float>>
        // These are "surface_boxes" because they're the box thats in the plane of the
        // cortical sheet (roughly xy) rather than the box in the slice plane (yz).
        vector<array<float,12>> surface_boxes;
        vector<array<float,3>> surface_box_centroids;
        array<float, 12> sbox;
        for (int i = 1; i < this->nFit; ++i) {
            // c1 x,y,z
            sbox[0] = this->layer_x;                // x
            sbox[1] = this->fitted_offset[i-1].x;  // y
            sbox[2] = this->fitted_offset[i-1].y; // z
            // c2 x,y,z
            sbox[3] = this->layer_x;              // x
            sbox[4] = this->fitted_offset[i].x;  // y
            sbox[5] = this->fitted_offset[i].y; // z
            // c3 x,y,z
            sbox[6] = this->layer_x+this->thickness; // x
            sbox[7] = this->fitted_offset[i].x;     // y
            sbox[8] = this->fitted_offset[i].y;    // z
            // c4 x,y,z
            sbox[9] = this->layer_x+this->thickness; // x
            sbox[10] = this->fitted_offset[i-1].x;  // y
            sbox[11] = this->fitted_offset[i-1].y; // z

            array<float, 3> sbox_centroid = MathAlgo<float>::centroid3D (sbox);
            surface_boxes.push_back (sbox);
            surface_box_centroids.push_back (sbox_centroid);
        }

        dname = frameName + "/fitted";
        df.add_contained_vals (dname.c_str(), fitted);

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
    }

    void printMeans (void) const {
        cout << "[";
        for (size_t j=0; j<this->means.size(); j++) {
            cout << this->means[j] << ",";
        }
        cout << "]" << endl << flush;
    }

    //! Recompute the polynomial fit
    void updateFit (void) {
        if (this->ct == CurveType::Poly) {
            this->updateFitPoly();
        } else {
            this->updateFitBezier();
        }
    }

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
                //cout << "Adding a PP point " << pt.x <<"," << pt.y << endl;
                user_points.push_back (make_pair(pt.x, pt.y));
            }

            BezCurve<double> bc;
            if (this->bcp.isNull()) {
                // No previous curves; fit just on user_points
                bc.fit (user_points);
                //cout << "fit with no previous curve..." << endl;
                this->bcp.addCurve (bc);
            } else {
                // Have previous curve, use last control of previous curve to make
                // smooth transition.
                //cout << "fit with previous curve..." << endl;
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
                //cout << "Adding a P point " << pt.x <<"," << pt.y << endl;
                user_points.push_back (make_pair(pt.x, pt.y));
            }
            BezCurve<double> bc;
            if (this->bcp.isNull()) {
                // No previous curves; fit just on user_points
                bc.fit (user_points);
                //cout << "fit P with no previous curve..." << endl;
                this->bcp.addCurve (bc);
            } else {
                BezCurve<double> last = this->bcp.curves.back();
                bc.fit (user_points, last);
                //cout << "fit P with previous curve..." << endl;
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
        Point2d fitsum;
        for (int i = 0; i < this->nFit; ++i) {
            this->fitted[i] = Point(coords[i].x(),coords[i].y());
            fitsum += Point2d(coords[i].x(),coords[i].y());
            this->tangents[i] = Point2d(tans[i].x(),tans[i].y());
            this->normals[i] = Point2d(norms[i].x(),norms[i].y());
        }
        this->fit_centroid = fitsum/this->nFit;
        cout << "Fit centroid: " << this->fit_centroid << endl;
        // Apply offset and scale
        for (int i = 0; i < this->nFit; ++i) {
            Point2d fd(this->fitted[i]);
            //cout << "Fit point: " << fd << endl;
            //cout << "Fit point - centroid: " << (fd - this->fit_centroid)
            this->fitted_offset[i] = (fd - this->fit_centroid) / this->pixels_per_mm;
        }
    }

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

    //! Re-compute the boxes from the curve
    void refreshBoxes (const int lenA, const int lenB) {
        this->refreshBoxes ((double)lenA, (double)lenB);
    }
    void refreshBoxes (const double lenA, const double lenB) {

        // cout << "Called, lenA=" << lenA << ", lenB=" << lenB << endl;
        if (this->ct == CurveType::Poly) {
            this->pointsInner = PolyFit::rotate (PolyFit::tracePolyOrth (this->pf, this->minX, this->maxX,
                                                                     this->nFit, lenA),
                                             this->theta);
            this->pointsOuter = PolyFit::rotate (PolyFit::tracePolyOrth (this->pf, this->minX, this->maxX,
                                                                      this->nFit, lenB),
                                              this->theta);
        } else {
            for (int i=0; i<this->nFit; i++) {
                //cout << "fitted[i]: " << fitted[i] << endl;
                //cout << "normals[i]: " << this->normals[i] << endl;
                //cout << "normals[i]*lenA: " << (this->normals[i]*lenA) << endl;
                Point2d normLenA = this->normals[i]*lenA;
                Point2d normLenB = this->normals[i]*lenB;
                //cout << "normals[i]*lenB: " << normLenB << endl;
                this->pointsInner[i] = this->fitted[i] + Point2i((int)normLenA.x, (int)normLenA.y);
                //cout << "pointsInner[i]: " << this->pointsInner[i] << endl;
                this->pointsOuter[i] = this->fitted[i] + Point2i((int)normLenB.x, (int)normLenB.y);
                //cout << "pointsOuter[i]: " << this->pointsOuter[i] << endl;
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
            //cout << "Define a box of width " << pointsInner[i] - pointsOuter[i]
            //     << ", length: " << pointsInner[i] - pointsInner[i+1] << endl;
            this->boxes[i] = pts;
        }
    }

    void toggleCurveType (void) {
        if (this->ct == CurveType::Poly) {
            this->ct = CurveType::Bezier;
        } else {
            this->ct = CurveType::Poly;
        }
    }

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

private:
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
