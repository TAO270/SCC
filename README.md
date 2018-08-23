# SCC(Sparse subspace clustering)
A re-implementation of MATLAB code with c++.<br>
A sparse subspace clustering method, and more details can be found here [paper](https://ieeexplore.ieee.org/abstract/document/6482137/) [MATLAB code](http://www.vision.jhu.edu/code/).<br>
## Dependency
You need to install eigen3 library for using because the data is calculated and stored in a matrix form.<br>
* Eigen3(Version 3.3.4): [Installation tutorial](http://eigen.tuxfamily.org/index.php?title=Main_Page)<br>
* libigl: a header-only library [Homepage](https://github.com/libigl/libigl).<br>
### In order to use the example, you need to add an OpenCV dependency.<br>
* OpenCV(Version 2.4)<br>
* CMake(Version 3.5)<br>
## Compile
    cd SCC
    cd build
    cmake ..
    make
You can find the Kmeans algorithm used in the code [here](https://github.com/michaelchughes/KMeansRex).

---
### There is a simple example classifying different features of images based on their trajectories over time  in images.<br>
The result shows that features of pedestrian's area of this image show different motion trajectories, which results in these feature points being classified into the same cluster.<br>
<br>
![](https://github.com/Markbess/SCC/blob/master/data/result.png)
