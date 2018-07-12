# SCC(Sparse subspace clustering)
A sparse subspace clustering method, and more details can be found here [paper](https://ieeexplore.ieee.org/abstract/document/6482137/) [MATLAB code](http://www.vision.jhu.edu/code/).<br>
A simple case in the paper is implemented with C++.<br>
## Dependency
You need to install eigen3 library for using because the data is calculated and stored in a matrix form.<br>
* Eigen3(version 3.3.4):[here](http://eigen.tuxfamily.org/index.php?title=Main_Page)<br>
* libigl:a header-only library.[here](https://github.com/libigl/libigl)<br>
In order to use the example, you need to add an OpenCV dependency.<br>
* OpenCV(version 2.4)<br>
* cmake(version 3.5)<br>
## Compile
    cd SCC
    cd build
    cmake ..
    make
You can find the Kmeans algorithm used in the code [here](https://github.com/michaelchughes/KMeansRex).
