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

    FrameData (const Mat frame) {
        this->frame = frame;
        axiscoefs.resize(2,0.);
        axis.resize(2);
        fitted.resize(nFit);
        origins.resize(nBins+1);
        tangents.resize(nBins+1);
        polyOrder = 3;
        nFit = 50;
        nBins = 50;
    };

    void removeLastPoint (void) {
        if (!P.empty()) {
            P.pop_back();
        }
    }

    void getBoxMeans (void) {
        means.resize(boxes.size());
        for(size_t i=0;i<boxes.size();i++){
            vector<double> boxVals = getPolyPixelVals(frame,boxes[i]);
            means[i] = 0.;
            for(size_t j=0;j<boxVals.size();j++){
                means[i] += boxVals[j];
            }
            means[i] /= (double)boxVals.size();
        }
    }

    void printMeans (void) const {
        cout<<"[";
        for(size_t j=0;j<means.size();j++){
            cout<<means[j]<<",";
        }
        cout<<"]"<<endl<<flush;
    }

    //! Recompute the polynomial fit
    void updateFit (void) {
        axiscoefs = polyfit(P,1);
        axis = tracePoly(axiscoefs,0,frame.cols,2);
        theta = atan(axiscoefs[1]);
        vector<Point> rotated = rotate(P,-theta);

        maxX = -1e9;
        minX = +1e9;
        for(size_t i=0;i<rotated.size();i++){
            if(rotated[i].x>maxX){ maxX = rotated[i].x; }
            if(rotated[i].x<minX){ minX = rotated[i].x; }
        }
        C = polyfit(rotated,polyOrder);

        fitted = rotate(tracePoly(C,minX,maxX,nFit),theta);
    }

    void refreshBoxes (const double lenA, const double lenB) {

        origins = rotate(tracePolyOrth(C,minX,maxX,nBins+1,lenA),theta);
        tangents = rotate(tracePolyOrth(C,minX,maxX,nBins+1,lenB),theta);

        boxes.resize(nBins);
        for(int i=0;i<nBins;i++){
            vector<Point> pts(4);
            pts[0]=origins[i];
            pts[1]=origins[i+1];
            pts[2]=tangents[i+1];
            pts[3]=tangents[i];
            boxes[i]=pts;
        }
    }

}; // FrameData
