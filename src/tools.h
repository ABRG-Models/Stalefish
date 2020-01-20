#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace cv;
using namespace std;

vector<double> polyfit(vector<Point> P, int n){

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

vector<Point> rotate(vector<Point> P, double theta){
    vector<Point> P2(P.size());
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);
    for(int i=0;i<P.size();i++){
        P2[i] = Point(P[i].x*cosTheta-P[i].y*sinTheta,P[i].x*sinTheta+P[i].y*cosTheta);
    }
    return P2;
}

vector<Point> tracePoly(vector<double> a,double minX,double maxX,int n){
    vector<Point> P(n);

    double increment = (maxX-minX)/(double)(n-1);
    double x = minX;
    double y = 0.;
    for (int i=0;i<n;i++){
        y = 0.;
        for(int j=0;j<a.size();j++){
            y += a[j]*pow(x,j);
        }
        P[i] = Point(x,y);
        x += increment;
    }
    return P;
}

vector<Point> tracePolyOrth(vector<double> a,double minX,double maxX,int n, double len){
    vector<Point> P(n);

    double increment = (maxX-minX)/(double)(n-1);
    double x = minX;
    double y = 0.;
    for (int i=0;i<n;i++){
        double t=0.;
        double y=0.;
        double p=0.;
        for(int j=0;j<a.size();j++){
            t += (double)j*a[j]*p;
            p = pow(x,j);
            y += a[j]*p;
        }
        P[i] = Point(x+len*cos(atan(t)+M_PI*0.5),y+len*sin(atan(t)+M_PI*0.5));
        x += increment;
    }
    return P;
}

vector<double> getPolyPixelVals(Mat frame, vector<Point> pp){

    Point pts[4] = {pp[0],pp[1],pp[2],pp[3]};
    Mat mask = Mat::zeros(frame.rows, frame.cols, CV_8UC3);
    fillConvexPoly(mask, pts, 4, cv::Scalar(255,255,255) );
    Mat result, resultGray;
    frame.copyTo(result,mask);
    cvtColor(result,resultGray,CV_BGR2GRAY);
    vector<Point2i> positives;
    findNonZero(resultGray, positives);
    vector<double> polyPixelVals(positives.size());
    for(int j=0;j<positives.size();j++){
        Scalar pixel = resultGray.at<uchar>(positives[j]);
        polyPixelVals[j] = (double)pixel.val[0]/255.;
    }
    return polyPixelVals;
}



class frameData{

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

    frameData(Mat frame){
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

    void removeLastPoint(void){
        if(P.size()){
            P.pop_back();
        }
    }

    void getBoxMeans(void){
        means.resize(boxes.size());
        for(int i=0;i<boxes.size();i++){
            vector<double> boxVals = getPolyPixelVals(frame,boxes[i]);
            means[i] = 0.;
            for(int j=0;j<boxVals.size();j++){
                means[i] += boxVals[j];
            }
            means[i] /= (double)boxVals.size();
        }
    }

    void printMeans(void){
        cout<<"[";
        for(int j=0;j<means.size();j++){
            cout<<means[j]<<",";
        }
        cout<<"]"<<endl<<flush;
    }

    void updateFit(void){
        double length = 40.;
        axiscoefs = polyfit(P,1);
        axis = tracePoly(axiscoefs,0,frame.cols,2);
        theta = atan(axiscoefs[1]);
        vector<Point> rotated = rotate(P,-theta);

        maxX = -1e9;
        minX = +1e9;
        for(int i=0;i<rotated.size();i++){
            if(rotated[i].x>maxX){ maxX = rotated[i].x; }
            if(rotated[i].x<minX){ minX = rotated[i].x; }
        }
        C = polyfit(rotated,polyOrder);

        fitted = rotate(tracePoly(C,minX,maxX,nFit),theta);
    }

    void refreshBoxes(double lenA, double lenB){

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


};

