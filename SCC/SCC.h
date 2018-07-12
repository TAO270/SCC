#ifndef SCC_H
#define SCC_H
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include<Your path/libigl/include/igl/cotmatrix.h>
#include<iostream>
#include<vector>
#include<time.h> 
namespace ygz{
using namespace Eigen;
using namespace std;
class SCC
{
private:
    const MatrixXf _Y; //Matrix<double, Dynamic, Dynamic>
private:
    //initial parameters
    int _alpha1;
    int _alpha2;
    int _alpha3;
    float mu1;
    float mu2;
   unsigned char _Iter;
    float _thr;
    const int mD;//rows
    const int mN;//cols
    MatrixXf mC;
    MatrixXf mA;
public:   
    SCC(MatrixXf& Y, int alpha1 = 800, int alpha2 = 800, int alpha3 = 800, float thr = 5*1e-4, unsigned char MaxIter =200);
    SCC();
    void AdmmLasso();
    void BuildAdjacency(MatrixXf &CKsym,float ro=0.7);
    void SpectralClustering(MatrixXf &CKsym,int group,ArrayXXd &centerpoint,ArrayXd &lable,int method=1);
    void initParameters();
};
}
#endif
