/*
 * Polynomial fit algorithms. We should move these into morphologica as they're
 * generally useful.
 */

#pragma once

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using cv::Mat;
using cv::Vec3f;
using cv::Point;
using cv::Scalar;
#include <vector>
using std::vector;
#include <math.h>
#include <iostream>
using std::cout;
using std::endl;

//! A set of functions for polynomial fitting
class PolyFit
{
public:
    //! Compute a polynomial fit of order n on the user-supplied points P. Return the
    //! parameters of the fit.
    static vector<double>
    polyfit (const vector<Point> P, int n) {
        int N = P.size();

        double X[2*n+1];
        for (int i=0;i<2*n+1;i++){
            X[i]=0;
            for (int j=0;j<N;j++){
                X[i]=X[i]+pow(P[j].x,i);
            }
        }
        double B[n+1][n+2];
        vector<double> a(n+1,0.);
        for (int i=0;i<=n;i++){
            for (int j=0;j<=n;j++){
                B[i][j]=X[i+j];
            }
        }
        double Y[n+1];
        for (int i=0;i<n+1;i++){
            Y[i]=0;
            for (int j=0;j<N;j++){
                Y[i]=Y[i]+pow(P[j].x,i)*P[j].y;
            }
        }
        for (int i=0;i<=n;i++){
            B[i][n+1]=Y[i];
        }
        n++;

        for (int i=0;i<n;i++){
            for (int k=i+1;k<n;k++){
                if (B[i][i]<B[k][i]){
                    for (int j=0;j<=n;j++){
                        double temp=B[i][j];
                        B[i][j]=B[k][j];
                        B[k][j]=temp;
                    }
                }
            }
        }
        for (int i=0;i<n-1;i++){
            for (int k=i+1;k<n;k++) {
                double t=B[k][i]/B[i][i];
                for (int j=0;j<=n;j++){
                    B[k][j]=B[k][j]-t*B[i][j];
                }
            }
        }
        for (int i=n-1;i>=0;i--){
            a[i]=B[i][n];
            for (int j=0;j<n;j++){
                if (j!=i){ a[i]=a[i]-B[i][j]*a[j]; }
            }
            a[i]=a[i]/B[i][i];
        }

        return a;
    }

    //! Compute @n points of the polynomial curve f(x) whose parameters are passed in
    //! with @a. x is varied from @minX to @maxX.
    static vector<Point>
    tracePoly (const vector<double> a,
               const double minX, const double maxX, const int n) {
        vector<Point> P(n);

        double increment = (maxX-minX)/(double)(n-1);
        double x = minX;
        double y = 0.;
        for (int i=0; i<n; i++) {
            y = 0.;
            for (size_t j=0; j<a.size(); j++) {
                y += a[j]*pow(x,j);
            }
            P[i] = Point(x,y);
            x += increment;
        }
        return P;
    }

    static vector<Point>
    tracePolyOrth (const vector<double> a,
                   const double minX, const double maxX,
                   const int n, const double len) {
        vector<Point> P(n);

        double increment = (maxX-minX)/(double)(n-1);
        double x = minX;
        for (int i=0;i<n;i++){
            double t=0.;
            double y=0.;
            double p=0.;
            for (size_t j=0; j<a.size(); j++) {
                t += (double)j*a[j]*p;
                p = pow(x,j);
                y += a[j]*p;
            }
            P[i] = Point (x+len*cos(atan(t)+M_PI*0.5),
                          y+len*sin(atan(t)+M_PI*0.5));
            x += increment;
        }
        return P;
    }

    //! Carry out a rotational transformation on the points P, returning the result.
    static vector<Point>
    rotate (const vector<Point> P, const double theta) {
        vector<Point> P2(P.size());
        double cosTheta = cos(theta);
        double sinTheta = sin(theta);
        for (size_t i=0; i<P.size(); i++) {
            P2[i] = Point(P[i].x*cosTheta-P[i].y*sinTheta,
                          P[i].x*sinTheta+P[i].y*cosTheta);
        }
        return P2;
    }
};

//! Some additional functions, useful for StaleFish.
class StaleUtil
{
public:
    //! Find a value for each pixel of image @frame within the box defined by @pp and
    //! return this in a vector of floats, without conversion
    static vector<float>
    getBoxedPixelVals (const Mat frame, const vector<Point> pp) {
        Point pts[4] = {pp[0],pp[1],pp[2],pp[3]};
        Mat mask = Mat::zeros(frame.rows, frame.cols, CV_8UC3);
        fillConvexPoly (mask, pts, 4, Scalar(255,255,255));
        Mat result, resultGray;
        frame.copyTo (result, mask);
        cvtColor (result, resultGray, cv::COLOR_BGR2GRAY);
        vector<cv::Point2i> positives;
        findNonZero (resultGray, positives);
        vector<float> boxedPixelVals (positives.size());
        for (size_t j=0; j<positives.size(); j++) {
            Scalar pixel = resultGray.at<uchar>(positives[j]);
            boxedPixelVals[j] = (float)pixel.val[0]/* /255.0f */; // fixme: Better to do the scaling elsewhere
        }
        return boxedPixelVals;
    }

    //! In a box, obtain colour values as BGR float triplets
    static vector<std::array<float, 3>>
    getBoxedPixelColour (const Mat frame, const vector<Point> pp) {
        Point pts[4] = {pp[0],pp[1],pp[2],pp[3]};
        Mat mask = Mat::zeros(frame.rows, frame.cols, CV_8UC3);
        fillConvexPoly (mask, pts, 4, Scalar(255,255,255));
        Mat result, resultFloat, resultGray;
        frame.copyTo (result, mask); // Note: 'result' will be in BGR format
        cvtColor (result, resultGray, cv::COLOR_BGR2GRAY);
        result.convertTo (resultFloat, CV_32FC3);
        vector<cv::Point2i> positives;
        findNonZero (resultGray/*Float*/, positives); // only for CV_8UC1 :( I want 'findNonBlack' on result
        vector<std::array<float, 3>> boxedPixelVals (positives.size());
        for (size_t j=0; j<positives.size(); j++) {
            Vec3f pixel = resultFloat.at<Vec3f>(positives[j]);
            // NB: This assumes image is in BGR format and we return in BGR format.
            boxedPixelVals[j][0] = pixel.val[0];
            boxedPixelVals[j][1] = pixel.val[1];
            boxedPixelVals[j][2] = pixel.val[2];
        }

        return boxedPixelVals;
    }
};
