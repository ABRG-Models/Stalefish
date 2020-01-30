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
#include "tools.h"
#include <morph/BezCurvePath.h>
using morph::BezCurvePath;
#include <morph/BezCurve.h>
using morph::BezCurve;
#include <morph/BezCoord.h>
using morph::BezCoord;

enum class CurveType {
    Poly,  // Variable order polynomial
    Bezier // Cubic Bezier
};

/*!
 * A class to hold a cortical section image, user-supplied cortex edge points and the
 * resulting fit (either polynomial, or, to come, Bezier curved).
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
    //! An auto-fitting Bezier curve
    BezCurve<double> bc;
    //! But the auto fit curve isn't good enough to do the whole cortical curve, so
    //! we'll need several in a BezCurvePath
    BezCurvePath<double> bcp;
    //! And then an index into bcp...
    int curvePathIndex;
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
    //! Number of bins to create for the fit (one less than nFit)
    int nBins;
    //! Number of points to create in the fit
    int nFit;
    //! A set of points created from the fit
    vector<Point> fitted;
    //! For point in fitted, the tangent at that location
    vector<Point2d> tangents;
    //! For point in fitted, the normal at that location
    vector<Point2d> normals;
    //! The axes for the fit
    vector<Point> axis;
    //! Coefficients for the axis
    vector<double> axiscoefs;
    //! origins for the lines making the box sides (some distance from the curve)
    vector<Point> pointsInner;
    //! endpoints for the lines making the box sides (a greater distance from the curve)
    vector<Point> pointsOuter;
    //! The boxes that are drawn and from which to sample the gene expression
    vector<vector<Point> > boxes;
    //! Do we show the fit lines or not?
    bool showFit = true;
    //! The image data, required when sampling the image in one of the boxes.
    Mat frame;
    //@}

    //! Constructor initializes default values
    FrameData (const Mat& fr) {
        this->frame = fr.clone();
        this->axiscoefs.resize (2, 0.0);
        this->axis.resize (2);
        // NB: Init these before the next three resize() calls
        this->nFit = 101;
        this->nBins = nFit-1;
        this->fitted.resize (this->nFit);
        this->pointsInner.resize (this->nFit);
        this->pointsOuter.resize (this->nFit);
        this->tangents.resize (this->nFit);
        this->normals.resize (this->nFit);
        this->polyOrder = 3;
    };

    string getFitInfo (void) const {
        stringstream ss;
        if (this->ct == CurveType::Poly) {
            ss << "Poly order: " << this->polyOrder << ", Bins: " << this->nBins;
        } else if (this->ct == CurveType::Bezier) {
            ss << "Bezier order: " << this->bc.getOrder() << ", Bins: " << this->nBins;
        } else {
            ss << "unknown";
        }
        return ss.str();
    }

    void removeLastPoint (void) {
        if (!this->P.empty()) {
            cout << "P.pop_back()" << endl;
            this->P.pop_back();
            cout << "P.popped_back()" << endl;
        } else {
            // P is empty. Go to previous curve...
            cout << "check prev. curve, pp_idx=" << pp_idx << endl;
            if (this->ct == CurveType::Bezier && this->pp_idx>0) {
                cout << "Set P to PP[" << pp_idx-1 << "]" << endl;
                this->P = this->PP[--this->pp_idx];
                cout << "PP.pop_back()" << endl;
                this->PP.pop_back();
                cout << "PP.popped_back()" << endl;
                //cout << "P.pop_back()" << endl;
                //this->P.pop_back(); // Why was this here?
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
                cout << "Adding a PP point " << pt.x <<"," << pt.y << endl;
                user_points.push_back (make_pair(pt.x, pt.y));
            }

            if (this->bcp.isNull()) {
                // No previous curves; fit just on user_points
                this->bc.fit (user_points);
            } else {
                // Have previous curve, use last control of previous curve to make
                // smooth transition.
                this->bc.fit (user_points, this->bcp.curves.back());
            }

            // Update this->fitted
            this->bcp.addCurve (this->bc);
        }

        if (this->P.size()>2) {
            vector<pair<double,double>> user_points;
            user_points.clear();
            for (auto pt : this->P) {
                cout << "Adding a P point " << pt.x <<"," << pt.y << endl;
                user_points.push_back (make_pair(pt.x, pt.y));
            }
            if (this->bcp.isNull()) {
                // No previous curves; fit just on user_points
                this->bc.fit (user_points);
            } else {
                this->bc.fit (user_points, this->bcp.curves.back());
            }

            this->bcp.addCurve (this->bc);
        }

        // Update this->fitted
        this->bcp.computePoints (static_cast<unsigned int>(this->nFit));
        vector<BezCoord<double>> coords = this->bcp.getPoints();
        vector<BezCoord<double>> tans = this->bcp.getTangents();
        vector<BezCoord<double>> norms = this->bcp.getNormals();
        for (int i = 0; i < this->nFit; ++i) {
            this->fitted[i] = Point(coords[i].x(),coords[i].y());
            this->tangents[i] = Point2d(tans[i].x(),tans[i].y());
            this->normals[i] = Point2d(norms[i].x(),norms[i].y());
            //cout << "norms[i]=" << norms[i] << endl;
            //cout << "normals[i]=" << this->normals[i] << endl; // ok
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
    void refreshBoxes (const double lenA, const double lenB) {

        cout << "Called, lenA=" << lenA << ", lenB=" << lenB << endl;
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

    void toggleShowFit (void) {
        this->showFit = this->showFit ? false : true;
    }


}; // FrameData
