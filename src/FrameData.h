#pragma once

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using cv::Mat;
using cv::Point;
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
    CurveType ct = CurveType::Poly;

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
    BezCurve bc;
    //! But the auto fit curve isn't good enough to do the whole cortical curve, so
    //! we'll need several in a BezCurvePath
    BezCurvePath bcp;
    //! And then an index into bcp...
    int curvePathIndex;
    //@}

    //! Attributes which pertain either to polynomial or Bezier curves
    //@{
    //! The vector of user-supplied points from which to make a curve fit
    vector<Point> P;
    //! The means computed for the boxes.
    vector<double> means;
    //! Number of bins to create for the fit
    int nBins;
    //! Number of points to create from the fit
    int nFit;
    //! A set of points created from the fit
    vector<Point> fitted;
    //! The axes for the fit
    vector<Point> axis;
    //! Coefficients for the axis
    vector<double> axiscoefs;
    //! origins for the lines making the box sides (some distance from the curve)
    vector<Point> origins;
    //! endpoints for the lines making the box sides (a greater distance from the curve)
    vector<Point> tangents;
    //! The boxes that are drawn and from which to sample the gene expression
    vector<vector<Point> > boxes;
    //! The image data, required when sampling the image in one of the boxes.
    Mat frame;
    //@}

    //! Constructor initializes default values
    FrameData (const Mat& fr) {
        this->frame = fr.clone();
        this->axiscoefs.resize (2, 0.0);
        this->axis.resize (2);
        // NB: Init these before the next three resize() calls
        this->nFit = 50;
        this->nBins = 50;
        this->fitted.resize (this->nFit);
        this->origins.resize (this->nBins+1);
        this->tangents.resize (this->nBins+1);
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
            this->P.pop_back();
        }
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
        if (this->P.size() < 2) {
            cout << "Too few points to fit" << endl;
            return;
        }
        vector<pair<float,float>> user_points;
        user_points.clear();
        for (auto pt : this->P) {
            user_points.push_back (make_pair(pt.x, pt.y));
        }
        this->bc.fit (user_points);

        // Update this->fitted
        this->bcp.addCurve (this->bc);
        this->bcp.computePoints (static_cast<unsigned int>(this->nFit));
        vector<BezCoord> coords = this->bcp.getPoints();
        int i = 0;
        for (BezCoord& bcoord : coords) {
            if (i>=this->nFit) { break; }
            this->fitted[i++] = Point(bcoord.x(),bcoord.y());
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

        if (this->ct == CurveType::Poly) {
            this->origins = PolyFit::rotate (PolyFit::tracePolyOrth (this->pf, this->minX, this->maxX,
                                                                     this->nBins+1, lenA),
                                             this->theta);
            this->tangents = PolyFit::rotate (PolyFit::tracePolyOrth (this->pf, this->minX, this->maxX,
                                                                      this->nBins+1, lenB),
                                              this->theta);
        } else {
            // We have this->fitted as a series of points, need to create origins and tangents...
            vector<BezCoord> tans = this->bcp.getTangents();
            int i = 0;
            for (BezCoord& bcoord : tans) {
                if (i>this->nBins) { break; }
                this->tangents[i] = this->fitted[i] + Point(bcoord.x(),bcoord.y());
                i++;
            }
        }

        // Make the boxes from origins and tangents
        this->boxes.resize (this->nBins);
        for (int i=0; i<this->nBins; i++) {
            vector<Point> pts(4);
            pts[0] = this->origins[i];
            pts[1] = this->origins[i+1];
            pts[2] = this->tangents[i+1];
            pts[3] = this->tangents[i];
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

}; // FrameData
