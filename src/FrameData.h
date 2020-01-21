#pragma once

class FrameData{

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

    FrameData(Mat frame){
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
        for(size_t i=0;i<boxes.size();i++){
            vector<double> boxVals = getPolyPixelVals(frame,boxes[i]);
            means[i] = 0.;
            for(size_t j=0;j<boxVals.size();j++){
                means[i] += boxVals[j];
            }
            means[i] /= (double)boxVals.size();
        }
    }

    void printMeans(void){
        cout<<"[";
        for(size_t j=0;j<means.size();j++){
            cout<<means[j]<<",";
        }
        cout<<"]"<<endl<<flush;
    }

    void updateFit(void){
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
