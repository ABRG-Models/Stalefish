#pragma once

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using cv::Mat;
using cv::Point;
using cv::Point2i;
using cv::Scalar;
#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::endl;
using std::flush;
#include "tools.h"

/*!
 * A class to hold a cortical section image, user-supplied cortex edge points and the
 * resulting fit (either polynomial, or, to come, Bezier curved).
 */
class FrameData
{
public:

    int polyOrder;
    int nBins;
    int nFit;

    double theta, minX, maxX;
    vector<double> means;
    vector<Point> P;
    vector<Point> fitted;
    vector<Point> axis;
    vector<Point> origins, tangents;
    vector<vector<Point> > boxes;
    vector<double> C;
    Mat frame;
    vector<double> axiscoefs;

    FrameData (const Mat fr) {
        this->frame = fr;
        this->axiscoefs.resize (2, 0.0);
        this->axis.resize (2);
        this->fitted.resize (this->nFit);
        this->origins.resize (this->nBins+1);
        this->tangents.resize (this->nBins+1);
        this->polyOrder = 3;
        this->nFit = 50;
        this->nBins = 50;
    };

    void removeLastPoint (void) {
        if (!this->P.empty()) {
            this->P.pop_back();
        }
    }

    void getBoxMeans (void) {
        this->means.resize (this->boxes.size());
        for (size_t i=0; i<this->boxes.size(); i++) {
            vector<double> boxVals = getPolyPixelVals (this->frame, this->boxes[i]);
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
        this->axiscoefs = polyfit (this->P, 1);
        this->axis = tracePoly (this->axiscoefs, 0, this->frame.cols, 2);
        this->theta = atan(this->axiscoefs[1]);
        vector<Point> rotated = rotate (this->P, -this->theta);

        this->maxX = -1e9;
        this->minX = +1e9;
        for (size_t i=0; i<rotated.size(); i++) {
            if (rotated[i].x > this->maxX) { this->maxX = rotated[i].x; }
            if (rotated[i].x < this->minX) { this->minX = rotated[i].x; }
        }
        this->C = polyfit (rotated, this->polyOrder);

        this->fitted = rotate (tracePoly (this->C, this->minX, this->maxX, this->nFit), this->theta);
    }

    void refreshBoxes (const double lenA, const double lenB) {

        this->origins = rotate (tracePolyOrth (this->C, this->minX, this->maxX,
                                               this->nBins+1, lenA),
                                this->theta);
        this->tangents = rotate (tracePolyOrth (this->C, this->minX, this->maxX,
                                                this->nBins+1, lenB),
                                 this->theta);

        this->boxes.resize(this->nBins);

        for (int i=0; i<this->nBins; i++) {
            vector<Point> pts(4);
            pts[0] = this->origins[i];
            pts[1] = this->origins[i+1];
            pts[2] = this->tangents[i+1];
            pts[3] = this->tangents[i];
            this->boxes[i] = pts;
        }
    }

}; // FrameData
