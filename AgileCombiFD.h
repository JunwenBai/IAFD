#ifndef AgileFD_H
#define AgileFD_H


#include <cstdio>
#include <time.h>
#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <string>
#include <armadillo>
#include <fstream>
#include "assert.h"
#include "spline.h"
#include <set>
#include <map>
#include "GaussianCompare.h"

using namespace std;
using namespace arma;

// record composition coordinates
class Coordinate {
public:
    string s;
    vector<double> v;
};

// record instance file
class Instance {
public:
    int N, L;
    mat XrdMat;
    vector<double> qvalues;
    vector<Coordinate> comp_xyz;
    string uuid;
};

// a class for phases
class Phase {
public:
    vector<double> xrd;
};

// a class recording proportion, shift, etc.
class PhaseMixture {
public:
    bool isInvolved;
    double proportion;
    double shift;
};

class Slice {
public:
    set<int> samples;
    double lowcost;
    double highcost;
    int phaseIndex;
};

// arguments and rets for multi-thread parallelization
class arg_struct {
public:
    int K;
    int M;
    int sample;
    double sumcost;
    double mipgap;
    vector<colvec> h;
};

// parameters of AgileCombiFD
class AgileCombiFDParameters {
public:
    int M, K;
    mat initW, initH;
    double MatchSigmaTol, MatchShiftTol; 
	vector<double> cvgRates;
	double convergenceRate;
	vector<int> shiftsNo;

    AgileCombiFDParameters(int _M, int _K) {
        M = _M;
        K = _K;
		cvgRates.clear();
    }
    
	AgileCombiFDParameters(int _M, int _K, double shift, double sigma) {
        M = _M;
        K = _K;
		MatchSigmaTol = sigma;
		MatchShiftTol = shift;
		cvgRates.clear();
    }

	AgileCombiFDParameters(int _M, int _K, double shift, double sigma, vector<int> _shiftsNo) {
		M = _M;
        K = _K;
        MatchSigmaTol = sigma;
        MatchShiftTol = shift;
		shiftsNo = _shiftsNo;
        cvgRates.clear();
	}

    AgileCombiFDParameters() {
        M = 10;
        K = 8;
    }
};

// a class recording initialization information
class SeedFreeze {
public:
    bool whetherValueInit;
    vector<int> valueFreeze;
    vector<vector<double> > valueSeeds;
    bool whetherSampleInit;
    vector<int> sampleFreeze;
    vector<int> sampleSeeds;
    vector<set<int> > samples;
    vector<double> Q;

	SeedFreeze(bool w1, vector<double> q, vector<int> f1, vector<vector<double> > s1, bool w2, vector<int> f2, vector<int> s2, vector<set<int> > s) {
		whetherValueInit = w1;
		valueFreeze = f1;
		valueSeeds = s1;
		whetherSampleInit = w2;
		sampleFreeze = f2;
		sampleSeeds = s2;
		samples = s;
        Q = q;
	}
};

// a class for AgileFD
class AgileCombiFD {
private:
    constexpr static double eps = 1e-9; // threshold
    AgileCombiFDParameters param;

public:
	
	AgileCombiFDParameters getParam() {
		return param;
	}

    // constructor
    AgileCombiFD(AgileCombiFDParameters _param);

    // print shifts information into an output file
    void printShiftsInfo(vector<vector<PhaseMixture> > phase_mixtures, vector<Coordinate> xyz, int M, vector<double> Qlogsteps, vector<mat> H, char*);

    // when we shift a matrix, we use zeros to fill the blank rows or column
    mat InsertRows(mat, int, double, bool);

    // reconstruct signals
    mat Rec(mat, vector< mat >, int);

    // calculate the KL divergence
	double KLdiv(mat, mat, bool);

	// calculate the L1
	double L1(mat, mat, bool);

	// calculate the L2
	double L2(mat, mat, bool);

    // apply updating rules to W and H
    mat KLupdate(mat, mat &, vector<mat > &, vector<int>,mat &, mat, vector<int>, bool, int);

    // initialize W
    void InitW(mat &);

    // initialize H
    void InitH(vector<mat> &);

    // nomalize W
    mat NormalizeW(mat W);

    // calculate the L1 loss of two respective columns of two matrices
    double diff(mat, int, mat, int);

    // at the end of this algorithm, normalize H and adjust W
    void NormalizeWH(mat &W, vector<mat > &H);
    
    // solver (main function)
    void Solve(Instance, double, int, char*, mat, bool, bool, SeedFreeze, double, mat &, vector<mat> &, mat &, bool, time_t, int, bool, vector<double>&, vector<double>&, bool, Mat<int>&, char* paramBuffer, vector<int> visitMarks, char* sticksfile);

    // print reconstructed signals (mainly for debugging)
    void printRec(mat D, vector<double> Q, mat Dlog, vector<double> Qlogsteps, mat DRec, char *);

    // print all the necessary information to the output file (phases, reconstructed signals, concentrations, shifts)
    void generateSolTxt(vector<Phase> phases, vector<vector<PhaseMixture> > phase_mixtures, int N, int M, int K, int L, int len, vector<double> Q, vector<double> Qlogsteps, char* sol_file, vector<mat> H, mat, mat, mat, mat, int, bool, Mat<int>&, char* paramBuffer, vector<int> visitMarks, char* sticksfile, string uuid);

};

#endif
