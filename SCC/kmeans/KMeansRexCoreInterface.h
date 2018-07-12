#pragma once

extern "C" {
	void RunKMeans(double *X_IN,  int N,  int D, int K, int Niter, \
	           int seed, char* initname, double *center, double *label);
	void SampleRowsPlusPlus(double *X_IN,  int N,  int D, int K, int seed, double *Mu_OUT);
}
