#include "SCC.h"
#include<iostream>
//kmeans++
#include "kmeans/KMeansRexCore.cpp"
namespace ygz{

  SCC::SCC ( MatrixXf &Y, int alpha1, int alpha2 , int alpha3,float thr , unsigned char MaxIter): _Y(Y),mD(Y.rows()),mN(Y.cols()),mC(mN,mN),mA(mN,mN)
{
    _alpha1=alpha1;
    _alpha2=alpha2;
    _alpha3=alpha3;
    _thr=thr;
    _Iter=MaxIter;
}

void SCC::initParameters()
{
    
    const MatrixXf &Y=_Y;
    float gamma=_alpha3/Y.cwiseAbs().colwise().sum().maxCoeff();//1-norm
    MatrixXf P(mD,mD+mN);
    P.block(0,0,mD,mN)=Y;
    P.block(0,mN,mD,mD)=Eigen::MatrixXf::Identity(mD,mD)*(1/gamma);//P=[Y eye(D)/gamma]
//computeLambda_mat
    MatrixXf T(mD+mN,mN);
    T=P.transpose()*Y;//T = P' * Y;
    T.block(0,0,mN,mN).diagonal().setZero();//T(1:N,:) = T(1:N,:) - diag(diag(T(1:N,:)));
    float lambda=T.cwiseAbs().colwise().maxCoeff().minCoeff();//lambda = min(max(T,[],1));
    mu1=_alpha1*(1/lambda);
    mu2=_alpha2*1;
    mA=(mu1*(_Y.transpose()*_Y)+mu2*Eigen::MatrixXf::Identity(mN,mN)).inverse();//N*N
}
void SCC::AdmmLasso()
{
    initParameters();
    //initialization 
    MatrixXf C1(MatrixXf::Zero(mN,mN));// zeros(N+D,N);
    //MatrixXf Lambda1(MatrixXf::Zero(D,N));//zeros(D,N);
    MatrixXf Lambda2(MatrixXf::Zero(mN,mN));//zeros(N+D,N);
    
    vector<float> error1; error1.reserve(_Iter); 
    error1.push_back(10*_thr);// err1 = 10*thr1
    
  //  vector<float> error2; error2.reserve(_Iter); 
  //  error2.push_back(10*_thr);// err2 = 10*thr2
    
    int i=0;
    MatrixXf Z(mN,mN);
    MatrixXf C2(mN,mN);
    const MatrixXf Ones_N_N(MatrixXf::Ones(mN,mN));
    MatrixXf Zeros_N_N(MatrixXf::Zero(mN,mN));
    //MatrixXf temp(mN,mN);
    //ADMM iterations
    while(error1[i]>_thr&&i<_Iter)
    {
        //updating Z
        Z=mA*( mu1*(_Y.transpose()*_Y) +mu2*(C1-Lambda2*(1/mu2))  );//N*N
        Z.diagonal().setZero();//Z = Z - diag(diag(Z));
       //updating C // C2= max( 0,(abs(Z+Lambda2/mu2) - 1/mu2*ones(N)) ) .* sign(Z+Lambda2/mu2);
        C2=(Zeros_N_N.cwiseMax( ( (Z+Lambda2*(1/mu2)).cwiseAbs().array()-(1/mu2) ).matrix() ) ).array()*( (Z+Lambda2*(1/mu2)).array().sign() );
        C2.diagonal().setZero();//C2 = C2 - diag(diag(C2));
        // updating Lagrange multipliers
        Lambda2 = Lambda2 + mu2 * (Z - C2);
        // computing errors
        float err=(Z-C2).array().cwiseAbs().maxCoeff();
         error1.push_back(err); 
        //update iterations
          C1=C2;
          i++;
    }
    //get Matrix mC
        mC=C2;
}
void SCC::BuildAdjacency(MatrixXf &CKsym,float ro)
{
    int N=mC.cols();
    MatrixXf Cp(MatrixXf::Zero(N,N));
    MatrixXf S(N,N);
    MatrixXi Ind(N,N);
    //discard 1-mRo data
    if(ro<1)
    {
    igl::sort(mC.cwiseAbs(),1,false,S,Ind);//sort cols-wise --descend
    float cL1=0;
    for(int i=0;i<N;i++)
        {
            cL1=S.col(i).sum();
            bool stop=false;
            float cSum=0;
            int t=0;
            while(!stop)
            {
                cSum=cSum+S(t,i);
                if(cSum>=ro*cL1){
                    stop=true; 
                   for(int idx=0;idx<t+1;idx++)
                   {
                       Cp.col(i)(Ind(idx,i))=mC.col(i)(Ind(idx,i));
                   }
                }
                t++;
            }
        }
    }
    else{
        Cp=mC;
    }
    Cp=Cp.cwiseAbs();
    for(int i=0;i<N;i++)
    {
        Cp.col(i)=Cp.col(i)*(1/(S.col(i)(1) + 1e-6) );//normalized
    }
    CKsym=Cp+Cp.transpose();
}
void SCC::SpectralClustering(MatrixXf &CKsym,int group,ArrayXXd &centerpoint,ArrayXd &lable,int method)
{
    const int N=CKsym.cols();
    MatrixXf Diag_D(N,N);
    Diag_D=MatrixXf((1/CKsym.colwise().sum().cwiseSqrt().array()).matrix().asDiagonal());
    MatrixXf LapN=MatrixXf::Identity(N,N)-Diag_D*CKsym*Diag_D;
    //SVD 
    MatrixXf kerN;
    if(N<=16)
    {
        //JacobiSVD
        JacobiSVD<MatrixXf> svd(LapN,ComputeThinV);
        kerN=svd.matrixV().block(0,N-group,N,group);
    }
    else
    {
        //BDCSVD
        BDCSVD<MatrixXf> svd(LapN,ComputeThinV);
        kerN=svd.matrixV().block(0,N-group,N,group);
    }
    //normalized kerN based rows  kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
    kerN.rowwise().normalize();
    Eigen::ArrayXXd Ker=kerN.cast<double>().array();
    if(method==1)
    {
    char intit[]="plusplus";
    RunKMeans(Ker.data(),Ker.rows(),Ker.cols(),group,1000,20,intit,centerpoint.data(),lable.data());
    }
    else
    {
     char intit[]="random";
    RunKMeans(Ker.data(),Ker.rows(),Ker.cols(),group,1000,20,intit,centerpoint.data(),lable.data());
    }
    
}

}
