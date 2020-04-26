#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <tclap/CmdLine.h>
#include <random>
#include <string>
#include "GaussianCompare.h"
#include "AgileCombiFD.h"
#include "PostProcessor.h"
#include <regex>
#include <map>

ILOSTLBEGIN
#define INF 1
#define Max_N 999


using namespace std;
using namespace TCLAP;


int N, n_threads;

mat D, W;
vector<mat> H;


// read instance files
Instance read_inst(char* filename, bool addNoise, double noiseStd) { // this function is mainly a parser
    ifstream f;
    f.open(filename);
    Instance inst; // inst comprises all the information that an instance file has
    char c;
    string buf,labstr,valstr,indstr;
    int sampleNo = 0; // the index of a sample
    string::size_type n;
    //For splitting input lines into substrings for variable labels, indices, and 
    regex line_re("\\s*([A-Za-z]+)_?([\\w\\[\\]*]*)\\s*=\\s*(.*[^\\s/])\\s*(?://.*)?");
    const string coordSetLabel("Composition");
    map<string,int> coordLabels;
    smatch match;
    while (getline(f,buf)) {
        if (regex_search(buf,match,line_re) && match.size()>1){
            labstr = match.str(1);
            indstr = match.str(2);
            valstr = match.str(3);
            if (labstr == "N") { // the number of sample points
                inst.N = stoi(valstr);
                if (inst.L){
                    inst.XrdMat = zeros<mat>(inst.L, inst.N);
                }
            } else if (labstr == "Q") {
                istringstream comma_sep(valstr);
                while(getline(comma_sep,buf,',')){
                    inst.qvalues.push_back(stod(buf));
                }
                inst.L = inst.qvalues.size();
                if (inst.N){
                    inst.XrdMat = zeros<mat>(inst.L, inst.N);
                }
            } else if (labstr == "I") { // XRD patterns
                sampleNo = stoi(indstr)-1;
                istringstream comma_sep(valstr);
                int i=0;
                while(getline(comma_sep,buf,',')){
                    inst.XrdMat(i,sampleNo) = max(0.0,stod(buf));
                    i++;
                }
                if (addNoise) {
                    double max_v = 0.;
                    for (int i = 0; i < inst.L; i++) max_v = std::max(max_v, inst.XrdMat(i, sampleNo));
                    std::default_random_engine generator;
                    std::normal_distribution<double> distribution(0., noiseStd * max_v);
                    for (int i = 0; i < inst.L; i++) {
                        double noise = distribution(generator);
                        double t = inst.XrdMat(i, sampleNo);
                        inst.XrdMat(i, sampleNo) = std::max(t+noise, 0.);
                    }
                }
            } else if (labstr == "UUID"){
                inst.uuid = valstr;
            } else if (labstr == coordSetLabel){
                istringstream comma_sep(valstr);
                int i=0;
                while(getline(comma_sep,buf,',')){
                    Coordinate v;
                    v.s = buf;
                    v.v.resize(inst.N);
                    coordLabels[v.s]=i;
                    inst.comp_xyz.push_back(v);
                    i++;
                }
            } else if (NULL != coordLabels[labstr]){
                istringstream comma_sep(valstr);
                int i=0;
                while(getline(comma_sep,buf,',')){
                    inst.comp_xyz[coordLabels[labstr]].v[i]=stod(buf);
                    i++;
                }
            }
        }
    }

    // print some parameters to make sure the parser processes the data correctly
//    inst.XrdMat /= 10.0;
    cout << "L: " << inst.L << endl; // the length of one XRD pattern
    cout << "N: " << inst.N << endl; // the number of sample points
    cout << "size of qvalues: " << inst.qvalues.size() << endl; // the size of qvalues array (should == L)
    cout << "XrdMat: (" << inst.XrdMat.n_rows << " " << inst.XrdMat.n_cols << ")" << endl; // D matrix: the collection of all the XRD patterns at each sample point
    cout << "UUID: "<<inst.uuid << endl;
    f.close();
    return inst;
}

mat getHumanIntelligence(char* filename, int N, int K, double beta_pow, vector<int> seeds, Slice &slice, bool whetherSlice) { // read human input(# of phases)

    mat beta = ones<mat>(K, N);
    beta = pow(beta, beta_pow);

    if (strlen(filename) == 0) return beta*beta_pow; // if there is no human input, then beta_pow is just a multiplier

    int specific_k = 0;
    for (int i = 0; i < seeds.size(); i++) {
        if (seeds[i] >= 186) specific_k = i;
    }

    FILE* file = freopen(filename, "r", stdin);
    if (file == NULL) return beta*beta_pow;

    double v;
    for (int i = 0; i < N; i++) {
        scanf("%lf", &v);
        for (int k = 0; k < K; k++) {
            if (v < 1e-6) {
                beta(k, i) = beta_pow*6; // if v < eps, the value of beta(k, i) is set to be the upper bound
                continue;
            }
            if (whetherSlice)
                if ( slice.samples.find(i) != slice.samples.end() ) v = 1.0;
            beta(k, i) = 1.0/v;
        }
    }

    // to enlarge the difference
    beta = pow(beta, beta_pow); // beta_pow now is the index

    /*
    std::set<int>::iterator it;
    if (whetherSlice)
    for (it=slice.samples.begin(); it!=slice.samples.end(); ++it) {
        int i = *it;
        for (int k = 0; k < K; k++) {
            if (k == specific_k) beta(k, i) *= slice.lowcost;
            else beta(k, i) *= slice.highcost;
        }
    }
    */

    fclose(stdin);
    return beta;

}

void initial_from_value(char * filename, vector<int> &freeze, vector<vector<double> > &seeds, bool &whetherInit, vector<set<int> > &samples, int N, vector<double> &Q) { // initialization based on given initial values(vectors)
    
    if (strlen(filename) == 0) { // users do not want to use this function
        whetherInit = false;
        return;
    }


    FILE* file = freopen(filename, "r", stdin);
    if (file == NULL) { // no such file
        whetherInit = false;
        return;
    }

    samples.clear();
    freeze.clear();
    whetherInit = true;
    char c;
    int K;
    string buf;

    // get values
    while ((c=getchar()) != EOF) {
        if (c == 10) continue;
        if (c == 13) {getchar();continue;}
        if (c == 'K') {
            char cc = getchar();
            if (cc == '=') cin >> K;
            else {
                getline(cin, buf);
                continue;
            }
        } else if (c == 'Q') { // seeds: initial values for seeding
            char cc = getchar();
            bool isXrd = true;
            bool canContinue = false;
            while (cc != '=') {
                if (cc == 10) {canContinue = true; break;}
                if (cc == 13) {canContinue = true; getchar(); break;}
                if (isdigit(cc)) cc = getchar();
                else {
                    isXrd = false;
                    break;
                }
            }
            if (canContinue) continue;
            if (!isXrd) {
                getline(cin, buf);
                continue;
            }

            Q.clear();
			while (true) {
                double v;
                scanf("%lf", &v);
                Q.push_back(v);
                cc = getchar();
                if (cc == 10) break;
                if (cc == 13) {getchar();break;}
            }
        } else if (c == 'B') { // seeds: initial values for seeding
            char cc = getchar();
            bool isXrd = true;
            bool canContinue = false;
            while (cc != '=') {
                if (cc == 10) {canContinue = true; break;}
                if (cc == 13) {canContinue = true; getchar(); break;}
                if (isdigit(cc)) cc = getchar();
                else {
                    isXrd = false;
                    break;
                }
            }
            if (canContinue) continue;
            if (!isXrd) {
                getline(cin, buf);
                continue;
            }
            
            vector<double> phase;
            while (true) {
                double v;
                scanf("%lf", &v);
                phase.push_back(v);
                cc = getchar();
                if (cc == 10) break;
                if (cc == 13) {getchar();break;}
            }
            seeds.push_back(phase);
        } else if (c == 'F') { // initial values for freezing
            
            char cc = getchar();
            bool isXrd = true;
            bool canContinue = false; // indicate whether this line is valid or not
            while (cc != '=') {
                if (cc == 10) {canContinue = true; break;}
                if (cc == 13) {canContinue = true; getchar(); break;}
                if (isdigit(cc)) cc = getchar();
                else {
                    isXrd = false;
                    break;
                }
            }
            if (canContinue) continue;
            if (!isXrd) {
                getline(cin, buf);
                continue;
            }
            
            int v;
            scanf("%d", &v);
            freeze.push_back(v);
            
        } else if (c == 'S') { // seeds: initial values for seeding

            char cc = getchar();
            bool isXrd = true;
            bool canContinue = false;
            while (cc != '=') {
                if (cc == 10) {canContinue = true; break;}
                if (cc == 13) {canContinue = true; getchar(); break;}
                if (isdigit(cc)) cc = getchar();
                else {
                    isXrd = false;
                    break;
                }
            }
            if (canContinue) continue;
            if (!isXrd) {
                getline(cin, buf);
                continue;
            }
            
            set<int> s;
            s.clear();
            while (true) {
                int v;
                scanf("%d", &v);
                s.insert(v-1);
                cc = getchar();
                if (cc == 10) break;
                if (cc == 13) {getchar();break;}
            }
            samples.push_back(s);

        } else {
            getline(cin, buf);
            continue;
        }
    }
}

double readDouble(string s) {
    char* str = (char*)s.c_str();
    double v;
    sscanf(str, "%lf", &v);
    return v;
}

void readDoubleArray(vector<double> &vec, string buf) {
    vec.clear();
    size_t found = buf.find(',');
    while (found != string::npos) {
        vec.push_back(readDouble(buf.substr(0, found)));
        buf.erase(0, found+1);
        found = buf.find(',');
    }
    vec.push_back(readDouble(buf));
}

void init_H(char* filename, vector<mat> &H, int N, int K, int M, bool &randinitH, char* pointFile) {
    if (strlen(filename) == 0) { // users do not want to use this function
        return;
    }
    FILE* file = freopen(filename, "r", stdin);
    if (file == NULL) { // no such file
        return;
    }
    bool canRand = false;
    if (strlen(pointFile) == 0) canRand = true;

    //----------------------------------------
    
    randinitH = false;
    H.clear();
    for (int m = 0; m < M; m++) {
        mat Hm = zeros<mat>(K,N);
        H.push_back(Hm);
    }
    char c;
    string buf;
    int k = 0, n = 0;

    while ((c = getchar()) != EOF) {
        
        if (c == 10) continue;
        if (c == 13) {getchar(); continue;}
        getline(cin, buf);
        buf = c+buf;

        if (c == 'H') {
            size_t found = buf.find("=");
            if (found == string::npos) continue;
            buf.erase(0, found+1);
            vector<double> temp;

            readDoubleArray(temp, buf);
            for (int m = 0; m < M; m++) {
                H[m](k,n) = temp[m];
                if (H[m](k,n) < 1e-6 and canRand) H[m](k,n) =(1e-3)*(rand() % (Max_N+1)+1.0)/(float)(Max_N+1);
            }
            k++;
            if (k == K) {
                k = 0;
                n++;
            }
        }
    }

}

void printVector(vector<int> v) { // for debugging
    int len = v.size();
    for (int i = 0; i < len; i++) printf("%d ", v[i]); printf("\n");
}

void initial_from_sample(char* filename, vector<int> &freeze, vector<int> &seeds, bool &whetherSampleInit, vector<set<int> > &samples, int N) { // initialize phases based on the XRD patterns at some sample point


    if (strlen(filename) == 0) { // users do not use this functionality
        whetherSampleInit = false;
        return;
    }
    FILE * file = freopen(filename, "r", stdin);
    if (file == NULL) { // no such file
        whetherSampleInit = false;
        return;
    }

    samples.clear();
    freeze.clear();
    seeds.clear();
    whetherSampleInit = true;
    char c;
    int K;
    string buf;

    // get values
    while ((c=getchar()) != EOF) {
        if (c == 10) continue;
        if (c == 13) {getchar();continue;}
        if (c == 'K') {
            char cc = getchar();
            if (cc == '=') cin >> K;
            else {
                getline(cin, buf);
                continue;
            }
        } else if (c == 'B') { // seeds: initial values for seeding
            char cc = getchar();
            bool isXrd = true;
            bool canContinue = false;
            while (cc != '=') {
                if (cc == 10) {canContinue = true; break;}
                if (cc == 13) {canContinue = true; getchar(); break;}
                if (isdigit(cc)) cc = getchar();
                else {
                    isXrd = false;
                    break;
                }
            }
            if (canContinue) continue;
            if (!isXrd) {
                getline(cin, buf);
                continue;
            }
            
            int v;
            scanf("%d", &v);
            seeds.push_back(v-1);

        } else if (c == 'F') { // initial values for freezing
            
            char cc = getchar();
            bool isXrd = true;
            bool canContinue = false; // indicate whether this line is valid or not
            while (cc != '=') {
                if (cc == 10) {canContinue = true; break;}
                if (cc == 13) {canContinue = true; getchar(); break;}
                if (isdigit(cc)) cc = getchar();
                else {
                    isXrd = false;
                    break;
                }
            }
            if (canContinue) continue;
            if (!isXrd) {
                getline(cin, buf);
                continue;
            }
            
            int v;
            scanf("%d", &v);
            freeze.push_back(v);
            
        } else if (c == 'S') { // seeds: initial values for seeding

            char cc = getchar();
            bool isXrd = true;
            bool canContinue = false;
            while (cc != '=') {
                if (cc == 10) {canContinue = true; break;}
                if (cc == 13) {canContinue = true; getchar(); break;}
                if (isdigit(cc)) cc = getchar();
                else {
                    isXrd = false;
                    break;
                }
            }
            if (canContinue) continue;
            if (!isXrd) {
                getline(cin, buf);
                continue;
            }
            
            set<int> s;
            s.clear();
            while (true) {
                int v;
                scanf("%d", &v);
                s.insert(v-1);
                cc = getchar();
                if (cc == 10) break;
                if (cc == 13) {getchar();break;}
            }
            samples.push_back(s);
        } else {
            getline(cin, buf);
            continue;
        }
    }

}

void addGibbs3(mat &W, vector<mat> &H, int N, int K, int M) { // the simplest way to add Gibbs phase rule: choose three phases with highest coefficients in H tensor
    for (int n = 0; n < N; n++) {
        vector<double> vec; // vec[k] represents the sum of H[m](k,n) which means the coefficient of the m-th shifted version of phase k at sample n for all m. Notice that k and n are given
        for (int k = 0; k < K; k++) {
            double sumc = 0.0;
            for (int m = 0; m < M; m++) {
                sumc += H[m](k,n); // calculate the sum
            }
            vec.push_back(sumc);
        }
        int k1, k2, k3; // top 3 phases
        double max = 0.0;
        for (int k = 0; k < K; k++) {
            if (vec[k] > max) {
                max = vec[k];
                k1 = k;
            }
        }
        max = 0.0;
        for (int k = 0; k < K; k++) {
            if (vec[k] > max && k != k1) {
                max = vec[k];
                k2 = k;
            }
        }
        max = 0.0;
        for (int k = 0; k < K; k++) {
            if (vec[k] > max && k != k1 && k != k2) {
                max = vec[k];
                k3 = k;
            }
        }
        for (int m = 0; m < M; m++) {
            for (int k = 0; k < K; k++) {
                if (k != k1 && k != k2 && k != k3) {
                    H[m](k,n) = 0.0;
                }
            }
        }
    }
}


void readSlice(char* filename, Slice &slice, bool &whetherSlice) {
    if (strlen(filename) == 0) { // users do not use this functionality
        whetherSlice = false;
        return;
    }
    FILE * file = freopen(filename, "r", stdin);
    if (file == NULL) { // no such file
        whetherSlice = false;
        return;
    }

    whetherSlice = true;
    slice.samples.clear();
    char c = getchar();
    while (c != '=') c = getchar();
    while (true) {
        int sample;
        scanf("%d", &sample);
        slice.samples.insert(sample-1);

        char c = getchar();
        if (c == 10 || c == EOF) break;
        if (c == 13) {getchar();break;}
    }

    c = getchar();
    while (c != '=') c = getchar();
    scanf("%lf", &slice.lowcost);

    c = getchar();
    while (c != '=') c = getchar();
    scanf("%lf", &slice.highcost);

    c = getchar();
    while (c != '=') c = getchar();
    scanf("%d", &slice.phaseIndex);
    slice.phaseIndex--;

}


void addSlice(char* filename, char* filename2, int N, int K, double beta_pow, Slice &slice, bool whetherSlice, mat& beta) { // read human input(# of phases)

    if (strlen(filename) == 0) return;
    if (strlen(filename2) == 0) {
        beta.fill(slice.lowcost);
        beta = pow(beta, beta_pow);
    }
   
    int specific_k = slice.phaseIndex; 
    std::set<int>::iterator it;
    if (whetherSlice)
    for (it=slice.samples.begin(); it!=slice.samples.end(); ++it) {
        int i = *it;
        for (int k = 0; k < K; k++) {
            if (k == specific_k) beta(k, i) *= slice.lowcost;
            else beta(k, i) *= slice.highcost;
        }
    }

}

bool Correct_point(char* pointFilename, vector<mat> &H, SeedFreeze &sf, int N, int K, int M, char* solfile, char* edgefile, double FLAGS_mipgap, double FLAGS_AlloyTol, Instance data, bool FLAGS_oneVersion, bool FLAGS_Alloy, vector<double> &Q, vector<double> &Qlogsteps, mat &Dlog, vector<int> &shiftsNo, mat &W) {
    if (strlen(pointFilename) == 0) return false;
    if (strlen(solfile) == 0) return false;

    FILE* file = freopen(pointFilename, "r", stdin);
    if (file == NULL) { // no such file
        return false;
    }

    // D matrix: the collection of XRD patterns at each sample point
    mat D = data.XrdMat;

    // convert Q into log space
    int L = data.L;
    double Qlogmin = log(data.qvalues[0]);
    double Qlogmax = log(data.qvalues[L-1]);
    double Qstepsize = (Qlogmax-Qlogmin)/(L-1);
    int rows = sf.valueSeeds.size();
    
    vector<double> Qlog;
    Q.clear();Qlogsteps.clear();Qlog.clear();

    for (int i = 0; i < L; i++) {
        Q.push_back(data.qvalues[i]);
    }
    for (int i = 0; i < L; i++) {
        Qlog.push_back(Qlogmin+i*Qstepsize);
        Qlogsteps.push_back(exp(Qlogmin+i*Qstepsize)); // Qlogsteps is a geometric sequence such that shifting Qlogsteps 1 position to the right means multiplying each qvalue by a constant
    }

    Dlog = zeros<mat>(L, N);

    // since we have new qvalues(Qlogsteps), we need to generate a new D matrix(Dlog) according to these new qvalues
    for (int i = 0; i < N; i++) {
        vector<double> xrd;
        xrd.clear();
        for (int j = 0; j < L; j++) {
            xrd.push_back(D(j, i));
        }
        vector<double> X(Q.begin(), Q.end());
        vector<double> Y(xrd.begin(), xrd.end());
//        Spline s(Q, xrd); // 1-d interpolation
        tk::spline s;
        s.set_points(X, Y);

        for (int j = 0; j < L; j++) {
            Dlog(j, i) = s(Qlogsteps[j]);
            if (Dlog(j, i) < 0.0) Dlog(j, i) = 0.0; // ensure non-negativity
        }
    }
    
    W = zeros<mat>(L, rows);
    for (int r = 0; r < rows; r++) {
        vector<double> xrd;
        xrd.clear();
        for (int l = 0; l < L; l++) {
            xrd.push_back(sf.valueSeeds[r][l]);
        }
        vector<double> X(sf.Q.begin(), sf.Q.end());
        vector<double> Y(xrd.begin(), xrd.end());
//        Spline s(X, Y);
        tk::spline s;
        s.set_points(X,Y);

        for (int l = 0; l < L; l++) {
            W(l,r) = s(Qlogsteps[l]);
            if (W(l,r) < 0.0) W(l,r) = 0.0;
        }
    }


    int num_of_pts;
    for (int n = 0; n < N; n++) {
        int cnt = 0;
        for (int k = 0; k < K; k++) {
            double tmp = 0.0;
            for (int m = 0; m < M; m++) tmp += H[m](k,n);
            if (tmp > 1e-9) cnt++;
        }
        shiftsNo[n] = 3-cnt;
    }
    
    scanf("%d", &num_of_pts);getchar();
    vector<int> samples;
    vector<set<int> > phaseSets;

    for (int rd = 0; rd < num_of_pts; rd++) {
        int sample;
        scanf("%d", &sample);getchar();sample--;
        samples.push_back(sample);
        set<int> phaseSet;
        char c;
        while (1) {
            int tmp;
            scanf("%d", &tmp);c=getchar();tmp--;
            phaseSet.insert(tmp);
            if (c == 10) break;
            if (c == 13) {getchar();break;}
            if (c == EOF) break;
        }
        phaseSets.push_back(phaseSet);
        
    }
    
    PostProcessor pp(edgefile, H, N, K, M, FLAGS_mipgap, FLAGS_AlloyTol, 0.05, Q, Qlogsteps);

    vector<set<int> >::iterator itr;
    int cnt = 0;
    for (itr = phaseSets.begin(); itr != phaseSets.end(); itr++) {
        set<int> phaseSet = *itr;
        int sample = samples[cnt];cnt++;
        pp.correctPoint(sample, data.comp_xyz, H, W, Dlog, shiftsNo, /*oneversion*/FLAGS_oneVersion, /*alloy*/FLAGS_Alloy, FLAGS_AlloyTol, 1e-6, phaseSet);
    }

    return true;
}

int counter = 0;

int main(int argc, char **argv)
{

    CmdLine cmd("AgileFD", ' ', "1.1");
    
	time_t init_time = time(NULL); // record the starting time

    // Arguments definition

    int FLAGS_k, FLAGS_m, FLAGS_time, FLAGS_seed, FLAGS_threads, FLAGS_AGrounds;
    string FLAGS_inst, FLAGS_sol, FLAGS_humanInput, FLAGS_valueInit, FLAGS_sampleInit, FLAGS_slice, FLAGS_neighbors, FLAGS_sticks, FLAGS_initH, FLAGS_pointCorrect;
    double FLAGS_c, FLAGS_beta, FLAGS_stepsize, FLAGS_sparsity, FLAGS_mipgap, FLAGS_AlloyTol, FLAGS_MatchShiftTol, FLAGS_MatchSigmaTol, FLAGS_noiseStd;
    bool FLAGS_shiftInfo, FLAGS_rec, FLAGS_Gibbs, FLAGS_oneVersion, FLAGS_Alloy, FLAGS_Connect, FLAGS_addNoise;

    try {
//        instArg, solArg, kArg, mArg, timeArg, cArg, seedArg, betaArg, shiftInfoArg, humanInputArg, recArg, initialArg, sampleInit, stepsizeArg, sparsityArg, mipgapArg, GibbsArg;
        
        ValueArg<string> instArg("", "inst", "Input instance file", true, "", "string", cmd);
        ValueArg<string> solArg("", "sol", "The output file name of the solution", true, "output.txt", "string", cmd);
        ValueArg<int> kArg("", "k", "The number of phases", true, 8, "int", cmd);
        ValueArg<int> mArg("", "m", "The number of possible different shifts", true, 10, "int", cmd);
        ValueArg<int> timeArg("", "time", "The maximum time(seconds) spent to train the model that you could accept", false, 10000, "int", cmd);
        ValueArg<double> cArg("", "c", "Related to termination criterion: In one iteration, if (old_cost-new_cost)<c*old_cost, then the loop terminates", false, 0.00008, "double", cmd);
        ValueArg<int> seedArg("", "seed", "random seed for the random number generator.(The default value -1 means time(0))", false, -1, "int", cmd);
        ValueArg<double> betaArg("", "beta", "The weighting coefficient of the sparsity term", false, 1.0, "double", cmd);
        SwitchArg shiftInfoArg("", "shiftInfo", "Whether you want to create a text file that contains information of shifts of each sample point", cmd, false);
        ValueArg<string> humanInputArg("", "humanInput", "Human Input txt file", false, "", "string", cmd);
        SwitchArg recArg("", "rec", "Whether output reconstructed signals", cmd, false);
        ValueArg<string> valueInitArg("", "valueInit", "Initialization file containing seeds for phases and phase freezing", false, "", "string", cmd);
        ValueArg<string> sampleInitArg("", "sampleInit", "The filename containing initialzation from single-phase sample points", false, "", "string", cmd);
        ValueArg<double> stepsizeArg("", "stepsize", "Initial stepsize. default shift value 0 means the user just want to use the std stepsize", false, 0.00000, "dobule", cmd);
        ValueArg<double> sparsityArg("", "sparsity", "The overall sparsity coefficient", false, 1.0, "double", cmd);
        ValueArg<double> mipgapArg("", "mipgap", "mipgap for MIP", true, 0.1, "double", cmd);
        SwitchArg GibbsArg("", "Gibbs", "whether to enforce Gibbs phase rule", cmd, false);
//        ValueArg<int> threadsArg("", "threads", "the number of threads to use", true, 12, "int", cmd);
        ValueArg<string> sliceArg("", "slice", "specify the slice constraint", false, "", "string", cmd);
        SwitchArg oneVersionArg("", "oneVersion", "whether to enforce one-used-shifted-version constraint", cmd, false);
        ValueArg<string> neighborsArg("", "neighbors", "the file telling the neighbors of every sample point", true, "", "string", cmd);
        ValueArg<double> AlloyTolArg("", "AlloyTol", "the tolerance for calculating shifts", false, 0.05, "double", cmd);
        ValueArg<string> sticksArg("", "sticks", "the file containing all the stick patterns", false, "", "string", cmd);
        ValueArg<int> AGroundsArg("", "AGrounds", "the rounds of AgileFD-Gibbs loop", false, 3, "int", cmd);
        SwitchArg AlloyAgentArg("", "AlloyAgent", "whether to enforce Alloy rule or not", cmd, false);
        SwitchArg ConnectAgentArg("", "ConnectAgent", "whether to enforce Connectivity rule or not", cmd, false);
        ValueArg<double> MatchShiftTolArg("", "MatchShiftTol", "shift tolerance when matching with ICSD patterns", false, 0.1, "double", cmd);
        ValueArg<double> MatchSigmaTolArg("", "MatchSigmaTol", "sigma boundary when matching with ICSD patterns", false, 0.25, "double", cmd);
        ValueArg<string> initHArg("", "initH", "the file for H initialization", false, "", "string", cmd);
        ValueArg<string> pointCorrectArg("", "pointCorrect", "correct or further refine a single point without changing anything else", false, "", "string", cmd);
        SwitchArg addNoiseArg("", "addNoise", "whether to add noise or not", cmd, false);
        ValueArg<double> noiseStdArg("", "noiseStd", "noise std dev", false, 0.5, "double", cmd);


    cmd.parse(argc, argv);
    
    FLAGS_inst = instArg.getValue();
    FLAGS_sol = solArg.getValue();
    FLAGS_k = kArg.getValue();
    FLAGS_m = mArg.getValue();
    FLAGS_time = timeArg.getValue();
    FLAGS_c = cArg.getValue();
    FLAGS_seed = seedArg.getValue();
    FLAGS_beta = betaArg.getValue();
    FLAGS_shiftInfo = shiftInfoArg.getValue();
    FLAGS_humanInput = humanInputArg.getValue();
    FLAGS_rec = recArg.getValue();
    FLAGS_valueInit = valueInitArg.getValue();
    FLAGS_sampleInit = sampleInitArg.getValue();
    FLAGS_stepsize = stepsizeArg.getValue();
    FLAGS_sparsity = sparsityArg.getValue();
    FLAGS_mipgap = mipgapArg.getValue();
    FLAGS_Gibbs = GibbsArg.getValue();
//    FLAGS_threads = threadsArg.getValue();
    FLAGS_slice = sliceArg.getValue();
    FLAGS_oneVersion = oneVersionArg.getValue();
    FLAGS_neighbors = neighborsArg.getValue();
    FLAGS_AlloyTol = AlloyTolArg.getValue();
    FLAGS_sticks = sticksArg.getValue();
    FLAGS_AGrounds = AGroundsArg.getValue();
    FLAGS_Alloy = AlloyAgentArg.getValue();
    FLAGS_Connect = ConnectAgentArg.getValue();
    FLAGS_MatchShiftTol = MatchShiftTolArg.getValue();
    FLAGS_MatchSigmaTol = MatchSigmaTolArg.getValue();
    FLAGS_initH = initHArg.getValue();
    FLAGS_pointCorrect = pointCorrectArg.getValue();
    FLAGS_addNoise = addNoiseArg.getValue();
    FLAGS_noiseStd = noiseStdArg.getValue();

    } catch (ArgException &e) {
        cerr << "error: " << e.error() << "for arg" << e.argId() << endl;
    }

	if (FLAGS_seed == -1) FLAGS_seed = time(0);

    cout << "inst: " << FLAGS_inst << endl;
    cout << "sol: " << FLAGS_sol << endl;
    printf("k: %d\n", FLAGS_k);
    printf("m: %d\n", FLAGS_m);
    printf("time: %d\n", FLAGS_time);
    printf("c: %.10lf\n", FLAGS_c);
    printf("seed: %d\n", FLAGS_seed);
    printf("beta: %lf\n", FLAGS_beta);
    cout << "shiftInfo: " << (FLAGS_shiftInfo?"true":"false") << endl;
    cout << "humanInput: " << FLAGS_humanInput << endl;
    cout << "rec: " << (FLAGS_rec?"true":"false") << endl;
    cout << "initial: " << FLAGS_valueInit << endl;
    cout << "sampleInit: " << FLAGS_sampleInit << endl;
    printf("stepsize: %lf\n", FLAGS_stepsize);
    printf("sparsity: %lf\n", FLAGS_sparsity);
    printf("mipgap: %lf\n", FLAGS_mipgap);
    cout << "Gibbs: " << (FLAGS_Gibbs?"true":"false") << endl;
    //cout << "threads: " << FLAGS_threads << endl;
    cout << "slice: " << FLAGS_slice << endl;
    cout << "oneVersion: " << (FLAGS_oneVersion?"true":"false") << endl;
    cout << "neighbors: " << FLAGS_neighbors << endl;
    printf("AlloyTol: %lf\n", FLAGS_AlloyTol);
    cout << "sticks: " << FLAGS_sticks << endl;
    printf("AGrounds: %d\n", FLAGS_AGrounds);
    cout << "Alloy: " << (FLAGS_Alloy?"true":"false") << endl;
    cout << "Connect: " << (FLAGS_Connect?"true":"false") << endl;
    printf("MatchShiftTol: %lf\n", FLAGS_MatchShiftTol);
    printf("MatchSigmaTol: %lf\n", FLAGS_MatchSigmaTol);
    cout << "initH: " << FLAGS_initH << endl;
    cout << "pointCorrect" << FLAGS_pointCorrect << endl;
    cout << "addNoise:" << (FLAGS_addNoise?"true":"false") << endl;
    cout << "noiseStd" << FLAGS_noiseStd << endl;

	
	char paramBuffer[1500];


    sprintf(paramBuffer, "inst: %s, k: %d, m: %d, time: %d, c: %.10lf, Gibbs: %d, AlloyTol: %lf, Alloy: %d, Connect: %d", FLAGS_inst.c_str(), FLAGS_k, FLAGS_m, FLAGS_time, FLAGS_c, (int)FLAGS_Gibbs, FLAGS_AlloyTol, (int)FLAGS_Alloy, (int)FLAGS_Connect);

    // setting running parameters
    AgileCombiFDParameters param(FLAGS_m, FLAGS_k, FLAGS_MatchShiftTol, FLAGS_MatchSigmaTol);
    AgileCombiFD iafd(param);
    if (FLAGS_seed == -1) srand(time(0));
    else srand(FLAGS_seed);


    // read instance file, more parameters
    Instance inst = read_inst((char*)(FLAGS_inst.c_str()), FLAGS_addNoise, FLAGS_noiseStd);
    N = inst.N;
    int K = FLAGS_k, M = param.M;
    vector<set<int> > samples;
    
    // value initialization
    vector<int> freeze;
    vector<vector<double> > seeds;
    bool whetherInit;
	vector<double> initQ;
    initial_from_value((char*)(FLAGS_valueInit.c_str()), freeze, seeds, whetherInit, samples, N, initQ);

    // sample initialization
    bool whetherSampleInit;
    vector<int> sampleFreeze;
    vector<int> sampleSeeds;
    initial_from_sample((char*)(FLAGS_sampleInit.c_str()), sampleFreeze, sampleSeeds, whetherSampleInit, samples, N);

    printf("sampleFreeze:\n");
    printVector(sampleFreeze);
    printf("sampleSeeds:\n");
    printVector(sampleSeeds);

    mat beta = ones<mat>(K, N); // sparsity

    // more parameters
    beta = beta * FLAGS_sparsity;
    SeedFreeze sf(whetherInit, initQ, freeze, seeds, whetherSampleInit, sampleFreeze, sampleSeeds, samples); // put all the parameters into sf

    H.clear();
    bool randinitH = true; // whether initialize H randomly

    init_H((char*)FLAGS_initH.c_str(), H, N, K, M, randinitH, (char*)FLAGS_pointCorrect.c_str()); // initialize H tensor using given file
 
    vector<double> Q, Qlogsteps; // qvalues array
    Mat<int> CC = ones<Mat<int> >(K, N); // boolean matrix: CC(k,n) indicates whether phase k is activated in sample n
    vector<int> shiftsNo(N, 0);

    // just for point correction
    if (Correct_point((char*)FLAGS_pointCorrect.c_str(), H, sf, N, K, M, (char*)FLAGS_initH.c_str(), (char*)FLAGS_neighbors.c_str(), FLAGS_mipgap, FLAGS_AlloyTol, inst, FLAGS_oneVersion, FLAGS_Alloy, Q, Qlogsteps, D, shiftsNo, W)) {

        AgileCombiFDParameters paramPrint(FLAGS_m, FLAGS_k, FLAGS_MatchShiftTol, FLAGS_MatchSigmaTol, shiftsNo);
        paramPrint.cvgRates = iafd.getParam().cvgRates;
        paramPrint.convergenceRate = iafd.getParam().convergenceRate;
        AgileCombiFD iafdPrint(paramPrint);

        vector<int> visitMarks;
        iafdPrint.Solve(inst, FLAGS_c, FLAGS_time, (char*)FLAGS_sol.c_str(), beta, FLAGS_shiftInfo, FLAGS_rec, sf, FLAGS_stepsize, W, H, D, false, init_time, 800, true, Q, Qlogsteps, /*FLAGS_Connect*/false, CC, paramBuffer, visitMarks, (char*)FLAGS_sticks.c_str());
        fprintf(stderr, "Point-correction is done.\n");

        return 0;
    }

    int rounds = FLAGS_AGrounds; // # of AgileFD-Gibbs rounds
    int batch = N/rounds+1; // # of samples to enforce Gibbs rule in each round

	PostProcessor pptmp;

    // AgileFD-Gibbs loop

    for (int iter = 0; iter < rounds; iter++) {

        if (iter > 0) randinitH = false; // after the first round, H shouldn't be randomly intialized again
        iafd.Solve(inst, FLAGS_c, FLAGS_time/rounds+1, (char*)FLAGS_sol.c_str(), beta, FLAGS_shiftInfo, FLAGS_rec, sf, FLAGS_stepsize, W, H, D, randinitH, init_time, 600, false, Q, Qlogsteps, false, CC, paramBuffer, pptmp.visitMarks, (char*)FLAGS_sticks.c_str());

        if (iter == rounds-1) break;

        // terminate if Gibbs not required
        if (FLAGS_Gibbs == false) {
            break;
        }

        
        fprintf(stderr, "Solution without Gibbs rule added is done.\n\n\n\n");
        fprintf(stderr, "Start to add Gibbs phase rule......\n");
        fprintf(stderr, "Start to generate args...\n");
        double sumcost = 0.0;


        fprintf(stderr, "preparation for MIP is ready!\n");
        time_t time0 = time(NULL);


        //create postprocessor to enforce Gibbs rule

        PostProcessor pp((char*)FLAGS_neighbors.c_str(), H,N,K, M, FLAGS_mipgap, FLAGS_AlloyTol, 0.05, Q, Qlogsteps);
        vector<int> shiftsNo(N, 0);

        #pragma omp parallel
        {
            #pragma omp for
            for (int i = iter*batch; i < min(iter*batch+batch, N); i++) {
                pp.addGibbs(i, inst.comp_xyz, H, W, D, shiftsNo, /*oneVersion*/false, /*alloy*/false, FLAGS_AlloyTol, 1e-6, false, CC);
            }
        }

        // update H
        vector<int>::iterator it;

        for (int m = 0; m < M; m++) {
            for (int i = iter*batch; i < min(iter*batch+batch, N); i++) {
                for (int k = 0; k < K; k++) {
                    H[m](k,i) = max(pp.args[i].h[m](k), 0.0);
                }
            }
        }
        
        /*
        for (int i = iter*batch; i < min(iter*batch+batch, N); i++) {
            fprintf(stderr, "H%d=", i);
            for (int k = 0; k < K; k++) {
                double tsum = 0.0;
                for (int m = 0; m < M; m++) tsum += H[m](k,i);
                fprintf(stderr, "%lf, ", tsum);
            }
            fprintf(stderr, "\n");
        }
        return 0;
        */

        time_t time1 = time(NULL); // current time
        fprintf(stderr, "MIP time cost: %d seconds\n", time1-time0);
        fprintf(stderr, "Finish!!!\n\n\n\n");
    }

    PostProcessor pp((char*)FLAGS_neighbors.c_str(), H,N,K, M, FLAGS_mipgap, FLAGS_AlloyTol, 0.05, Q, Qlogsteps);
    
    shiftsNo.clear();
    for (int n = 0; n < N; n++) shiftsNo.push_back(0);


    // Extended Gibbs rule
    if (FLAGS_Gibbs) {
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < N; i++) {
                pp.addGibbs(i, inst.comp_xyz, H, W, D, shiftsNo, FLAGS_oneVersion, FLAGS_Alloy, FLAGS_AlloyTol, 1e-6, false, CC);
            }
        }

        for (int m = 0; m < M; m++) {
            for (int i = 0; i < N; i++) {
                for (int k = 0; k < K; k++) {
                    H[m](k,i) = max(pp.args[i].h[m](k), 0.0);
                }
            }
        }
    }
    

    // enforce connectivity
    
    //FLAGS_Connect = false;
    if (FLAGS_Connect) {
        fprintf(stderr, "\n\n----------------------\n");
        fprintf(stderr, "enforce connectivity!\n");
        fprintf(stderr, "------------------------\n\n");
        
        
        // if single phase region connectivity is broken after enforcing phase field connectivity, we enforce it again. Thus here is a loop.
        while (1) {

            pp.connectphase(H, N, K, M, CC);
            pp.connect(H, N, K, M, CC);
            
            for (int n = 0; n < N; n++) shiftsNo[n] = 0;
            
            #pragma omp parallel
            {
                #pragma omp for
                for (int i = 0; i < N; i++)
                    pp.addGibbs(i, inst.comp_xyz, H, W, D, shiftsNo, FLAGS_oneVersion, FLAGS_Alloy, FLAGS_AlloyTol, 1e-6, true, CC);
            }

            for (int m = 0; m < M; m++)
                for (int i = 0; i < N; i++)
                    for (int k = 0; k < K; k++)
                        H[m](k,i) = max(pp.args[i].h[m](k), 0.0);

            if (pp.checkCC(H, N, K, M)) break;
        }
    }
    
    fprintf(stderr, "\n---------------------\nstart to print\n---------------------\n\n");
	

    // Print solutions
    
    AgileCombiFDParameters paramPrint(FLAGS_m, FLAGS_k, FLAGS_MatchShiftTol, FLAGS_MatchSigmaTol, shiftsNo);
    paramPrint.cvgRates = iafd.getParam().cvgRates;
    paramPrint.convergenceRate = iafd.getParam().convergenceRate;
    AgileCombiFD iafdPrint(paramPrint);

    
    iafdPrint.Solve(inst, FLAGS_c, FLAGS_time, (char*)FLAGS_sol.c_str(), beta, FLAGS_shiftInfo, FLAGS_rec, sf, FLAGS_stepsize, W, H, D, false, init_time, 800, true, Q, Qlogsteps, FLAGS_Connect, CC, paramBuffer, pp.visitMarks, (char*)FLAGS_sticks.c_str());
    fprintf(stderr, "Solution with Gibbs rule added is done.\n");

    return 0;
}
