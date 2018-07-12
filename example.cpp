#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "SCC/SCC.h"
#include <sstream>
#include <string>
#include <memory>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include<time.h> 
using namespace Eigen;
using namespace std;
using namespace cv;
void showimage(Mat &image,vector<Point2f> &points,ArrayXd &label)
{
    vector<cv::Scalar> color;
    color.push_back(cv::Scalar(255,0,0));//B
    color.push_back(cv::Scalar(0,255,0));//G
    color.push_back(cv::Scalar(0,0,255));//R
    color.push_back(cv::Scalar(255,255,0));
    color.push_back(cv::Scalar(0,255,255));
    color.push_back(cv::Scalar(84,46,8));
    color.push_back(cv::Scalar(240,32,160));
    //check size 
    if(points.size()!=label.rows())
    {
        cout<<"the number of keyPoints is not equal to the size of label"<<endl;
    }
    else //show
    {
        for(size_t i=0;i<points.size();i++)
        {
            Point2f pt=points[i];
            double id=label(i);
            cv::rectangle(image,cv::Point2f(pt.x-2,pt.y-2),cv::Point2f(pt.x+2,pt.y+2),color.at(size_t(id)),-1);  
        }
       imshow("SCC_for_C++",image);
       waitKey(0);
    }
}
int main()
{
    ifstream f;
     ifstream point;
    f.open("../data/data.txt");
    point.open("../data/point.txt");
    if(!f.is_open()){
        cout<<"open error"<<endl;
        return -1;
    }
     MatrixXf Y(40,450);
     string temp;
     stringstream ss;
    int row=0;
     while(std::getline(f,temp,'\n'))
     {
         ss<<temp;
         float number=0;

         int col=0;
         while(ss>>number)
         {
             Y.row(row)(col)=number;
             col=col+1;
         }
         ss.clear();
         ss.str("");
         row=row+1;
         
     }
     //point
     vector<Point2f> pp;
     pp.reserve(Y.cols());
     ss.clear();
     ss.str("");
    float x=0;
    float y=0;
     while(std::getline(point,temp,'\n'))
     {
         ss<<temp;
         ss>>x;
         ss>>y;
         pp.push_back(Point2f(x,y));
         ss.clear();
         ss.str("");
     }
//SCC
     double start_time=clock();
    shared_ptr<ygz::SCC> scc=shared_ptr<ygz::SCC>(new ygz::SCC(Y));
    scc->AdmmLasso();

    MatrixXf CKsym(Y.cols(),Y.cols());
    scc->BuildAdjacency(CKsym);

    ArrayXXd centerpoint= Eigen::ArrayXXd::Zero(5,5);
    ArrayXd label=Eigen::ArrayXd::Zero(CKsym.rows());;
    scc->SpectralClustering(CKsym,5,centerpoint,label);
    double finsh_time=clock();
    cout<<"run time: "<<(finsh_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;
    //show image and point 
    Mat image=imread("../data/scc_test.png");
    showimage(image,pp,label);
}
