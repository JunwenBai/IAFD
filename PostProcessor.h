#include <ilcplex/ilocplex.h>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <iterator>
#include <algorithm>
#include "AgileCombiFD.h"

using namespace std;
using namespace arma;

class PostProcessor {

private:
    constexpr static double eps = 1e-9;

public:
    int N, K, M;
    double mipgap, tol;
    bool whethernb; //indicates if file is read correctly
    vector<vector<int> > neighbors;
    arg_struct args[400];
    mat avgshifts;
    bool first;
    vector<double> Q;
    vector<double> Qlogsteps;
    vector<int> visitMarks;

    // 1st constructor
    PostProcessor();

    // 2nd constructor
    PostProcessor(int K_, int M_,  double mipgap_, double tol);

    // 3rd constructor
    PostProcessor(char* filename,vector<mat> H ,int N_, int K_, int M_, double mipgap_, double tol, double shifteps, vector<double>, vector<double>);

    // read *_edges.txt file
    void readNeighbors(char* filename, vector<vector<int> > &nb, bool &whethernb);

    // enforce Gibbs phase rule
    void addGibbs(int sample, vector<Coordinate> &xyz, vector<mat> &H, mat & W, mat & D,  vector<int>&, bool onev, bool alloy, double epsilon, double shifteps, bool whetherCC, Mat<int> &CC);

    // calculate # of violates at some sample
    int nodeViol(int sample);

    // delete invalid connected components
    void delCC(int k);

    // check connecitivity
    bool checkCC(vector<mat> &H, int N, int K, int M);

    // enforce phase field connectivity
    void connect(vector<mat> &H, int N, int K, int M, Mat<int> &CC);

    // dfs function used to enforce phase field connectivity
    void govisit(int sample, double &sumv, vector<int> &visit, int cnt, vector<mat> &H, int k, map<set<int>, int> &combmap, map<int, set<int> > &combdict, vector<int> &mark, set<int>& losers);

    // enforce single phase region connecitivity
    void govisitphase(int sample, double &sumv, vector<int> &visit, int cnt, vector<mat> &H, int k);

    // dfs function used to enforce single phase region connectivity
    void connectphase(vector<mat> &H, int N, int K, int M, Mat<int> &CC);

    // correct point
    void correctPoint(int sample, vector<Coordinate> &xyz, vector<mat> &H, mat& W, mat& D,  vector<int> &shiftsNo, bool onev, bool alloy, double epsilon, double alloyeps, set<int> phaseSet);

};
