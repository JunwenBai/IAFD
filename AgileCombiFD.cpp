#include "AgileCombiFD.h"
#define Max_N 999

using namespace std;
using namespace arma;


// this function is seldom used for now
void AgileCombiFD::printShiftsInfo(vector<vector<PhaseMixture> > phase_mixtures, vector<Coordinate> xyz, int M, vector<double> Qlogsteps, vector<mat> H, char* sol_file) { // print shifts information
    char str[1000];
    int length = strlen(sol_file);
    // generate the output file name of "shifts_info" output file
    for (int i = 0; i < length; i++) {
        if (sol_file[i] == '.' && sol_file[i+1] == 't' && sol_file[i+2] == 'x' && sol_file[i+3] == 't') {str[i]='\0';break;}
        str[i] = sol_file[i];
    }
    strcat(str, "_Shifts_Information.txt");
    // open the output file
    freopen(str, "w", stdout);

    int ndim = xyz.size();

    int K = H[0].n_rows; // K: the number of phases
    int N = H[0].n_cols; // N: the number of sample points
    printf("M=%d\n", M); // M: the nubmer of shifted versions
    printf("K=%d\n", K);
    printf("N=%d\n", N);
    if (ndim == 3) {
        printf("\n// (x = %s, y = %s, z = %s)   %s,%s,%s are elements from instance file,", xyz[0].s.c_str(), xyz[1].s.c_str(), xyz[2].s.c_str(), xyz[0].s.c_str(), xyz[1].s.c_str(), xyz[2].s.c_str()); // xyz[i].s is the name of the i-th element
        printf(" or (x, y)   x, y are coordinates of decomposition space\n");
    } else {
        printf("\n// (x, y)   x, y are coordinates of decomposition space\n");
    }
 
    // ternary system coordinates
    for (int i = 0; i < ndim; i++) {
        if (i == 0) printf("x=");
        if (i == 1) printf("y=");
        if (i == 2) printf("z=");
        for (int j = 0; j < N; j++) {
            printf("%lf", xyz[i].v[j]);
            if (j == N-1) printf("\n");
            else printf(",");
        }
    }

    // print Qvalues
    printf("\n// Qvalues\n");
    printf("q=");
    int L = Qlogsteps.size();
    for (int i = 0; i < L; i++) {
        printf("%lf", Qlogsteps[i]);
        if (i == L-1) printf("\n");
        else printf(",");
    }

    // print H tensor
    printf("\n// Activation matrix H :  H[m](i,j) means the activation value for phase i at sample point j with shift m(or say, the m-th layer of the tensor)\n");
    for (int m = 0; m < M; m++) {
        for (int k = 0; k < K; k++) {
            printf("H[%d](%d,*)=", m, k);
            for (int n = 0; n < N; n++) {
                printf("%lf", H[m](k,n));
                if (n == N-1) printf("\n");
                else printf(",");
            }
        }
    }
    printf("\n");
    for (int n = 0; n < N; n++) {
        printf("avgH(%d, *)=", n);
        for (int k = 0; k < K; k++) {
            double sum = 0.0, sumw = 0.0;
            for (int m = 0; m < M; m++) {
                sum += H[m](k, n) * m;
                sumw += H[m](k, n);
            }
            double v;
            if (sumw > eps) v = sum/sumw;
            else v = 0.0;
            printf("%lf", v);
            if (k == K-1) printf("\n");
            else printf(",");
        }
    }

    fclose(stdout);
}

// constructor
AgileCombiFD::AgileCombiFD(AgileCombiFDParameters _param) {
    param = _param;
}

mat AgileCombiFD::InsertRows(mat W, int rowcount, double val, bool top = true) { // if top==true, row0 is inserted into the beginning. Otherwise, row0 is inserted into the end.
    rowvec row0(W.n_cols); // one row filled with zeros of length K
    row0.fill(val); // row0=[val,val,...,val]. Usually, val == 0.
    for (int r = 0; r < rowcount; r++) {
        if (top) {
            W.insert_rows(0, row0); // insert row0 to the beginning 
        } else {
            W.insert_rows(W.n_rows, row0); // insert row0 to the end
        }
    }
    return W;
}

// reconstruct signals
mat AgileCombiFD::Rec(mat W, vector<mat> H, int M) { // Phi is an array containing shift values w.r.t # of positions to shift(not the real shift values)
    int L = W.n_rows; // the length of one XRD pattern
    int K = W.n_cols; // # of phases
    int N = H[0].n_cols; // # of sample points
    M = H.size(); // # of shifted version
//    assert((int)Phi.size() == (int)M); // check whether Phi is valid or not

    mat Rec = zeros<mat>(L, N); // reconstructed signals (reconstructed D matrix)

    for (int i = 0; i < M; i++) {
        int p = i; // shift down p positions
        mat WW = W.rows(0,L-p-1);
        WW = InsertRows(WW, p, 0.0); // top p lines are filled with zeros

        Rec += WW * H[i]; // one shfited version of W multiplies the corresponding layer of H tensor
    }
    return Rec;
}

// calculate KL-divergence
double AgileCombiFD::KLdiv(mat D, mat Rec, bool show = false) { // "show" indicates whether to show some intermediate results (mainly for debugging)
    int L = D.n_rows;
    int N = D.n_cols;

    double ckl = 0.0; // KL divergence(cost)
    for (int j = 0; j < N; j++) {
        double loss = 0.;
        for (int i = 0; i < L; i++) {
            double Dij = D(i,j); // the original value at position (i,j)
            double Rij = Rec(i,j); // the reconstructed value at position (i,j)
            
            // print intermediate results
            if (show) printf("(%.10lf, %.10lf, %.10lf)\n", Dij, Rij, log((Dij+eps)/(Rij+eps)));
            if (show) printf("(%.10lf, %.10lf, %.10lf)\n",  (Dij+eps)/(Rij+eps), log((Dij+eps)/(Rij+eps)), Dij*log((Dij+eps)/(Rij+eps)));

            loss += Dij*log((Dij+eps)/(Rij+eps)) - Dij + Rij; // KL divergence
//            ckl += (Dij-Rij)*(Dij-Rij);
            if (show) printf("%d: %.10lf\n", i, ckl);
        }
        ckl += loss;
    }

    return ckl;
}

double AgileCombiFD::L1(mat D, mat Rec, bool show = false) { // "show" indicates whether to show some intermediate results (mainly for debugging)
    int L = D.n_rows;
    int N = D.n_cols;
    double ckl = 0.0; // KL divergence(cost)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < N; j++) {
            double Dij = D(i,j); // the original value at position (i,j)
            double Rij = Rec(i,j); // the reconstructed value at position (i,j)
            ckl += abs(Dij-Rij);
        }
    }

    return ckl;
}

double AgileCombiFD::L2(mat D, mat Rec, bool show = false) { // "show" indicates whether to show some intermediate results (mainly for debugging)
    int L = D.n_rows;
    int N = D.n_cols;
    double ckl = 0.0; // KL divergence(cost)
    for (int j = 0; j < N; j++) {
        double loss = 0.;
        for (int i = 0; i < L; i++) {
            double Dij = D(i,j); // the original value at position (i,j)
            double Rij = Rec(i,j); // the reconstructed value at position (i,j)
            loss += (Dij-Rij)*(Dij-Rij);
        }
        ckl += loss;
    }

    return sqrt(ckl);
}

mat AgileCombiFD::KLupdate(mat D, mat &W, vector<mat> &H, vector<int> Phi, mat &VRec, mat beta, vector<int> freeze, bool whetherInit, int iter_cnt) { // updating W and H
    // default meanings
    int L = W.n_rows;
    int K = W.n_cols;
    int N = H[0].n_cols;
    int M = H.size();

    double cost_old = KLdiv(D, VRec); // old cost

    vector<mat> H_old; // a copy of H tensor
    for (int i = 0; i < M; i++) {
        mat Hi = H[i];
        H_old.push_back(Hi);
    }
    mat W_old = W; // a copy of W
	
    W = NormalizeW(W);

    /*
        START to update! Please refer to the technical report for more details.
    */
    
    mat VR = D / (VRec+eps); // an auxiliary matrix

    mat O = ones<mat>(L, N); // an all-1 matrix
    
    vector<mat> H_x;
    vector<mat> H_y;
    vector<mat> GradH;


    for (int i = 0; i < M; i++) {
        int p = Phi[i];
        mat Hi_x = W.rows(0, L-p-1).t()*VR.rows(p, L-1);
        H_x.push_back(Hi_x);
        mat Hi_y = W.rows(0, L-p-1).t()*O.rows(p, L-1);
        H_y.push_back(Hi_y);

        mat Gradi = Hi_x / (Hi_y+beta);
//        mat Gradi = Hi_x / (Hi_y+0.5); // do not erase
//        mat Gradi = Hi_x / (Hi_y); // do not erase
        GradH.push_back(Gradi);
    }

    double nuH = 1.0;
    double accel = 1.0;
    int pindex = 0;
    for (int i = 0; i < M; i++) {
        H[pindex] = H_old[pindex] % (pow(GradH[pindex], nuH));
        pindex++;
    }
    VRec = Rec(W, H, M);
    double cost = KLdiv(D, VRec);

    cost_old = cost;

    // update W
    VR = D / (VRec + eps); // auxilliary matrix
//    VR = D;
//    O = VRec+eps;

    mat W_x = zeros<mat>(L, K); // auxiliary matrix (numerator in the updating rule)
    mat W_y = zeros<mat>(L, K); // auxiliary matrix (denominator)


    double h[50];
    mat H1, W1;
    for (int i = 0; i < M; i++) {
        int p = Phi[i];

        // numerator
        mat Wxp = VR.rows(p, L-1) * H[i].t();
        Wxp = InsertRows(Wxp, p, 0.0, false);
        mat Wyp = O.rows(p, L-1) * H[i].t();
        Wyp = InsertRows(Wyp, p, 0.0, false); 
        
        H1 = zeros<mat>(L, K);
        H1 = Wyp % W;
        for (int k = 0; k < K; k++) {
            h[k] = 0;
            for (int j = 0; j < L; j++) {
                h[k] += H1(j, k);
            }
        }
        W1 = zeros<mat>(L, K);
        for (int l = 0; l < L; l++) {
            for (int k = 0; k < K; k++) {
                W1(l, k) = W(l,k) * h[k];
            }
        }
        W_x = W_x+Wxp+W1;
        
        // denominator
        W1 = zeros<mat>(L, K);
        W1 = Wxp % W;
        for (int k = 0; k < K; k++) {
            h[k] = 0;
            for (int l = 0; l < L; l++) {
                h[k] += W1(l, k);
            }
        }
        H1 = zeros<mat>(L,K);
        for (int l = 0; l < L; l++) {
            for (int k = 0; k < K; k++) {
                H1(l, k) = W(l, k) * h[k];
            }
        }

        W_y = W_y+Wyp+H1;
    }

    mat GradW = W_x / (W_y+eps);
    for (int i = 0; i < freeze.size(); i++) {
        if (whetherInit && freeze[i] == 1) { // freezing phases
            for (int l = 0; l < L; l++) GradW(l,i) = 1;
        }
    }

    W = W % (pow(GradW, nuH));
	
    VRec = Rec(W, H, M);
    cost = KLdiv(D, VRec);

    cost_old = cost;
    
    return VRec;
}

void AgileCombiFD::InitW(mat &W) { // initialize W
    int L = W.n_rows;
    int K = W.n_cols;

    for (int j = 0; j < K; j++) {
        for (int i = 0; i < L; i++)
            W(i,j) = (rand() % (Max_N+1)+1)/float(Max_N+1); // to avoid initial value 0
    }

}

void AgileCombiFD::InitH(vector<mat> &H) {
    int K = H[0].n_rows;
    int N = H[0].n_cols;
    int M = H.size();

    for (int m = 0; m < M; m++) {

        for (int k = 0; k < K; k++) {
            for (int n = 0; n < N; n++) {
                H[m](k, n) = (rand() % (Max_N+1)+1)/float(Max_N+1); // to avoid initial value 0
            }
        }

    }
}

mat AgileCombiFD::NormalizeW(mat W) { // normalize W phase by phase
    int L = W.n_rows;
    int K = W.n_cols;
    for (int j = 0; j < K; j++) {
        double W2norm = norm(W.cols(j, j), "fro"); // calculate frobenius norm of the j-th column
        for (int i = 0; i < L; i++) {
			if (W(i,j) < eps) continue;
            W(i,j) /= W2norm;
        }
    }
    return W;
}

void AgileCombiFD::NormalizeWH(mat &W, vector<mat> &H) { // normalize H
    int L = W.n_rows;
    int K = H[0].n_rows;
    int N = H[0].n_cols;
    int M = H.size();

    for (int k = 0; k < K; k++) {
        double maxsumH = 0.0;
        for (int n = 0; n < N; n++) {
            double sumH = 0.0;
            for (int m = 0; m < M; m++) {
                sumH += H[m](k, n);
            }
            maxsumH = max(maxsumH, sumH);
        }

        // normalize maxsumH to be 1
        for (int m = 0; m < M; m++) {
            for (int i = 0; i < N; i++) {
                double valH = H[m](k,i);
                if (maxsumH < eps) H[m](k,i) = 0;
                else H[m](k,i) = valH / maxsumH;
            }
        }

        // to maintain the resulting matrix of W*H, multiply W(i,k) by some maxsumH
        for (int i = 0; i < L; i++) {
            double valW = W(i,k);
            W(i,k) = valW * maxsumH;
        }
    }
}

// print reconstructed signals to a file
void AgileCombiFD::printRec(mat D, vector<double> Q, mat Dlog, vector<double> Qlogsteps, mat DRec, char* sol_file) {
    char str[1000];
    int length = strlen(sol_file);
    for (int i = 0; i < length; i++) {
        if (sol_file[i] == '.' && sol_file[i+1] == 't' && sol_file[i+2] == 'x' && sol_file[i+3] == 't') {str[i] = '\0'; break;} // remove suffix
        str[i] = sol_file[i];
    }
    strcat(str, "_Ori_Rec.txt"); // add new suffix and file type
    freopen(str, "w", stdout);

    int L = Dlog.n_rows; // the length of one XRD pattern
    int N = Dlog.n_cols; // # of sample points
    vector<double> xrd;
    printf("N=%d\n", N);
    printf("L=%d\n", L);
    printf("Q=");
    for (int i = 0; i < L; i++) {
        printf("%lf", Qlogsteps[i]); // qvalues
        if (i == L-1) printf("\n");
        else printf(",");
    }

    for (int n = 0; n < N; n++) {
        printf("O%d=", n);
        for (int j = 0; j < L; j++) { // the original XRD pattern at sample point n
            printf("%lf", Dlog(j, n));
            if (j == L-1) printf("\n");
            else printf(",");
        }
    }
    
    for (int n = 0; n < N; n++) {
        printf("R%d=", n);
        for (int j = 0; j < L; j++) { // the reconstructed XRD pattern at sample point n
            printf("%lf", DRec(j, n));
            if (j == L-1) printf("\n");
            else printf(",");
        }
    }
    
    fclose(stdout);
}

// calculate the L1 loss of two respective columns of two matrices
double AgileCombiFD::diff(mat W, int k, mat Dlog, int n) {
    int L = W.n_rows;
    k = 186, n = 186;
    double cost = 0.0;
    for (int l = 0; l < L; l++) {
        cost += abs(W(l,k)-Dlog(l,n));
    }
    return cost;
}

void AgileCombiFD::Solve(Instance data, double convergenceRate, int time_limit, char* sol_file, mat beta, bool showShift, bool showRec, SeedFreeze sf, double givenstepsize, mat &W, vector<mat> &H, mat &Dlog, bool randini, time_t init_time, int iter_limit, bool forprint, vector<double> &Q, vector<double> &Qlogsteps, bool whetherCC, Mat<int> &CC, char* paramBuffer, vector<int> visitMarks, char* sticksfile) {
    time_t now_time; // current time
    time_t local_init_time;

    local_init_time = time(NULL);
    this->param.convergenceRate = convergenceRate;
    
    // default meanings
    int N = data.N;
    int M = param.M;
    int K = param.K;
    int L = data.L;

    // D matrix: the collection of XRD patterns at each sample point
    mat D = data.XrdMat;

    // convert Q into log space
    double Qlogmin = log(data.qvalues[0]);
    double Qlogmax = log(data.qvalues[L-1]);
    double Qstepsize = (Qlogmax-Qlogmin)/(L-1);
    fprintf(stderr, "std Qstepsize = %.50lf\n", Qstepsize);

    int len;
    if (givenstepsize > eps) { // reset stepsize using a given stepsize
        len = (int)((Qlogmax-Qlogmin)/log(1+givenstepsize)+1.0);
        Qstepsize = (Qlogmax-Qlogmin)/(len-1);
    } else {
        len = L;
    }

    fprintf(stderr, "Qstepsize = %.50lf\n", Qstepsize);
    fprintf(stderr, "Length of qvalues array:\nOld len = %d, New len = %d\n", L, len);
 

    Q.clear();Qlogsteps.clear();
//    vector<double> Q;
    vector<double> Qlog;
//    vector<double> Qlogsteps;

    for (int i = 0; i < L; i++) {
        Q.push_back(data.qvalues[i]);
    }
    for (int i = 0; i < len; i++) {
        Qlog.push_back(Qlogmin+i*Qstepsize);
        Qlogsteps.push_back(exp(Qlogmin+i*Qstepsize)); // Qlogsteps is a geometric sequence such that shifting Qlogsteps 1 position to the right means multiplying each qvalue by a constant
    }

    Dlog = zeros<mat>(len, N);

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

        for (int j = 0; j < len; j++) {
            Dlog(j, i) = s(Qlogsteps[j]);
            if (Dlog(j, i) < 0.0) Dlog(j, i) = 0.0; // ensure non-negativity
        }
    }

    mat DRec;
    vector<int> Phi;
	//Init Phi
    for (int i = 0; i < M; i++) {
        Phi.push_back(i); // generate shift values with regard to how many positions to shift
    }

	int iter_cnt;

    // if this function is just used for printint, then skip the following part.
    if (!forprint) {  ///////////

        // Init matrix W
        W = zeros<mat>(len,K);
        InitW(W); // random initialization

        // convert seeding values into log space
        int rows = sf.valueSeeds.size();
        mat Wseed = zeros<mat>(len,rows);
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

            for (int l = 0; l < len; l++) {
                Wseed(l,r) = s(Qlogsteps[l]);
                if (Wseed(l,r) < 0.0) Wseed(l,r) = 0.0;
            }
        }
        
        vector<int> freeze;
        bool whetherInit = false;

        // if the user specifies to seed with certain values,
        if (sf.whetherValueInit) {
            int rows = Wseed.n_cols;
            for (int r = 0; r < min(rows, K); r++) {
                for (int l = 0; l < len; l++) {
                    W(l,r) = Wseed(l,r);
                }
            }
            freeze = sf.valueFreeze;
            whetherInit = sf.whetherValueInit;
        }


        // if the user specifies to seed from certain samples
        if (sf.whetherSampleInit) {
            int rows = sf.sampleSeeds.size();
            int nn = min(rows, K);
            for (int r = 0; r < nn; r++) {
                for (int l = 0; l < len; l++) {
                    W(l,r) = Dlog(l,sf.sampleSeeds[r]);
                }
            }
            freeze = sf.sampleFreeze;
            whetherInit = sf.whetherSampleInit;
        }
       

        W = NormalizeW(W); // step 1: W becomes W_tilde
        

        // Init matrices H
        if (randini == true) { // randini indicates whether you want a random initialization or not
            // if yes,
            H.clear();
            for (int m = 0; m < M; m++) {
                mat Hm = zeros<mat>(K,N);
                H.push_back(Hm);
            }
            InitH(H); // initialize H randomly
        } else {
            // if no, we initialize H based on the given H
            
            
            for (int n = 0; n < N; n++) {
                int cnt = 0;
                for (int k = 0; k < K; k++) {
                    double tmp = 0.0;
                    for (int m = 0; m < M; m++) {
                        tmp += H[m](k,n);
                    }
                    if (tmp > eps) { // this means phase k plays a role at sample point n
                        for (int m = 0; m < M; m++)
                            if (H[m](k,n) < eps) 
                                H[m](k,n) = (1e-3)*(rand() % (Max_N+1)+1.0)/(float)(Max_N+1);
                    } else {
                        for (int m = 0; m < M; m++) H[m](k,n) = 0.0;
                    }
                }
            }
            

        }

        if (whetherInit) {
            int cols = Wseed.n_cols;

            int kk = sf.samples.size();
            for (int k = 0; k < min(kk, K); k++) {
                for (int n = 0; n < N; n++) {
                    if (sf.samples[k].find(n) != sf.samples[k].end()) continue;
                    for (int m = 0; m < M; m++) {
                        H[m](k,n) = 0.0;
                    }
                }
            }
        }

        
        if (sf.whetherSampleInit) { // initialize some values in H based on sampleSeeds
            int rows = sf.sampleSeeds.size();
            for (int r = 0; r < min(rows, K); r++) {
                for (int k = 0; k < K; k++) {
                    if (k == r) {
                        continue;
                    }
                    for (int m = 0; m < M; m++) H[m](k, sf.sampleSeeds[r]) = 0.0;
                }
            }
        }

        DRec = Rec(W, H, M); // reconstructed D matrix
        double cost = KLdiv(Dlog, DRec); // calculate cost

        fprintf(stderr, "cost: %lf\n", cost);

        now_time = time(NULL); // current time
        int time_cost = now_time - local_init_time; // duration time
        fprintf(stderr, "Before KLupdate, we have spent %d seconds.\n", time_cost);

        double finalConvRate;
        iter_cnt = 0;
        while (true) {
            DRec = KLupdate(Dlog, W, H, Phi, DRec, beta, freeze, whetherInit, iter_cnt); // up-to-date reconstructed D matrix
            
            double new_cost = KLdiv(Dlog, DRec); // new cost
            if (iter_cnt % 25 == 0)
                fprintf(stderr, "%d: old: %lf, new: %lf", iter_cnt, cost, new_cost);
            now_time = time(NULL);
            time_cost = now_time - local_init_time;

            // based on cost and new_cost, we determine whether to stop or not
            // termination criterion: 1. reach the convergence point  2. running time expires
            if (abs(cost-new_cost) < cost*convergenceRate+eps || time_cost > time_limit /*|| iter_cnt >= iter_limit*/) {
                fprintf(stderr, "\n\ncost-newcost: %lf\n", cost-new_cost);
                finalConvRate = (double)(cost-new_cost)/cost;
                fprintf(stderr, "Final Convergence Rate: %.12lf\n", finalConvRate);
                param.cvgRates.push_back(finalConvRate);

                //fprintf(stderr, "iter_limit: %d\n", iter_limit);
                fprintf(stderr, "%d: old: %lf, new: %lf\n", iter_cnt, cost, new_cost);
                break;
            }
            else {
                cost = new_cost; // new cost becomes old cost
                if (iter_cnt % 50 == 0) fprintf(stderr, "   Total time cost until now: %d s.", time_cost);
                if (iter_cnt % 25 == 0) fprintf(stderr, "\n");
            }
            iter_cnt++;
        }

        NormalizeWH(W, H); // normalize H

    } //////////

    DRec = Rec(W, H, M);

    if (showRec) { // whether to print reconstructed signals
        printRec(D, Q, Dlog, Qlogsteps, DRec, sol_file);
    }



    vector<Phase> phases; // XRD patterns of phases
    vector<vector<PhaseMixture> > phase_mixtures; // shifts, proportion

    for (int k = 0; k < K; k++) {
        vector<double> xrd;
        for (int j = 0; j < len; j++) {
            xrd.push_back(W(j, k));
        }

        // convert qvalues from Qlogsteps array back to original Q array using interpolation
        vector<double> X(Qlogsteps.begin(), Qlogsteps.end());
        vector<double> Y(xrd.begin(), xrd.end());

//        Spline s(X, Y);
        tk::spline s;
        s.set_points(X,Y);

        Phase phase;
        for (int j = 0; j < L; j++) {
            phase.xrd.push_back(s(Q[j]));
        }

        phases.push_back(phase);

        vector<PhaseMixture> mixtures;
        for (int i = 0; i < N; i++) {
            double sumH = 0.0;
            double shiftH = 0.0;

            for (int m = 0; m < M; m++) {
                // Notice that Phi[m] == m
                sumH += H[m](k,i); // total proportion
                double shift = Qlogsteps[m]/Qlogsteps[0]; // generate the real shift values
                shiftH += shift * H[m](k,i);
            }
            if (sumH < eps) shiftH = 0.0; // phase k doesn't appear at sample point i
            else shiftH /= sumH; // weighted shift value

            PhaseMixture mixture;
            mixture.isInvolved = (sumH >= eps); // whether phase k involves at sample point i
            mixture.proportion = (sumH >= eps)?sumH:0.0; // proportion
            mixture.shift = shiftH; // shift
            mixtures.push_back(mixture);
        }
        phase_mixtures.push_back(mixtures);
    }
    
    // generate *_output.txt

	if (forprint)
		generateSolTxt(phases, phase_mixtures, N, M, K, L, len, Q, Qlogsteps, sol_file, H, D, Dlog, DRec, W, init_time, whetherCC, CC, paramBuffer, visitMarks, sticksfile, data.uuid);

    // generate *_Ori_Rec.txt (seldom used for now)
    if (showShift) {
        printShiftsInfo(phase_mixtures, data.comp_xyz, M, Qlogsteps, H, sol_file);
    }
}

void AgileCombiFD::generateSolTxt(vector<Phase> phases, vector<vector<PhaseMixture> > phase_mixtures, int N, int M, int K, int L, int len, vector<double> Q, vector<double> Qlogsteps, char* sol_file, vector<mat> H, mat D, mat Dlog, mat DRec, mat W, int init_time, bool whetherCC, Mat<int> &CC, char* paramBuffer, vector<int> visitMarks, char* sticksfile, string uuid) {
	
	mat Drec = Rec(W, H, M);
    double printLoss = KLdiv(Dlog, Drec);
    double L1Loss = L1(Dlog, Drec);
    double L2Loss = L2(Dlog, Drec);


//	DRec = zeros<mat>(L, N);
	double* sampleLoss = new double[N+1];
    double* sampleTot = new double[N+1];
    time_t now_time = time(NULL);
 
	DRec = zeros<mat>(L, N);
	for (int n = 0; n < N; n++) {
		vector<double> xrd;
		for (int l = 0; l < L; l++) xrd.push_back(Drec(l,n));
		vector<double> X(Qlogsteps.begin(), Qlogsteps.end());
		vector<double> Y(xrd.begin(), xrd.end());
		
		tk::spline s;
        s.set_points(X,Y);
		
		for (int l = 0; l < L; l++) DRec(l,n) = s(Q[l]);
	}

    double origin_kl = KLdiv(D, DRec);
    double origin_l1 = L1(D, DRec);
    double origin_l2 = L2(D, DRec);

    DRec = Rec(W, H, M);
    for (int n = 0; n < N; n++) {
        double loss = 0.0;
        double tot = 0.0;
        for (int l = 0; l < len; l++) {
            double d = Dlog(l,n), dr = DRec(l,n);
//            loss += d*log((d+eps)/(dr+eps))-d+dr;
            loss += abs(d-dr);
            tot += abs(d);
        }
        sampleLoss[n] = loss;
        sampleTot[n] = tot;
    }

	int below25 = 0, below50 = 0;
	double minRatio = 100.0;

	for (int n = 0; n < N; n++) {
		double ratio = sampleLoss[n]/sampleTot[n];
		if (1-ratio < 0.25) below25++;
		if (1-ratio < 0.50) below50++;
		minRatio = min(minRatio, 1-ratio);
	}

	char filename[100];
    char *extptr = strrchr(sol_file,'.');
    int extloc = strlen(sol_file);
    if (NULL!=extptr) extloc -= strlen(extptr);

    // stats file filename
	sprintf(filename, "%.*s_stats_KL=%.6e%s", extloc, sol_file, printLoss, extptr);

	delete [] sampleLoss;
	delete [] sampleTot;

    FILE *fsol = fopen(sol_file, "w");
    FILE *fstats = fopen(filename, "w");

	fprintf(fsol,"Description=IAFD solution for %s\n", paramBuffer);
    if (uuid.length()) fprintf(fsol,"UUID=%s\n", uuid.c_str());

	fprintf(fstats,"Description=IAFD solution for %s\n", paramBuffer);
    if (uuid.length()) fprintf(fstats,"UUID=%s\n", uuid.c_str());
    fprintf(fstats,"\n");

    fprintf(fstats,"time cost: %d\n", now_time-init_time); // running time
    fprintf(fstats,"KL=%.10lf\n", printLoss);
    fprintf(fstats,"L1=%.10lf\n", L1Loss);
    fprintf(fstats,"L2=%.10lf\n", L2Loss);
    fprintf(fstats,"\n");
     fprintf(fstats,"origin_KL=%.10lf\n", origin_kl);
    fprintf(fstats,"origin_L1=%.10lf\n", origin_l1);
    fprintf(fstats,"origin_L2=%.10lf\n", origin_l2);

    fprintf(fsol,"K=%d\n", K);

    fprintf(fstats,"Number of loops: %d\n", param.cvgRates.size());

    for (int i = 0; i < param.cvgRates.size(); i++) {
        fprintf(fstats,"%.12lf ", param.cvgRates[i]);
        if (param.cvgRates[i] > param.convergenceRate) fprintf(stderr," [WARNING!!! FAIL to reach convergence rate!!!]  ");
    }
    fprintf(fstats,"\n");

    fprintf(fsol,"\n//List of solution models\n");
    fprintf(fsol,"Params=[Q,R],[Q,B,C,S],[Q,B,H]\n");

    double xrd[L+1];

    fprintf(fsol,"\n// Phase pattern (basis)\n");
    // print Q values
    fprintf(fsol,"Q=");
    for (int i = 0; i < L; i++) {
        double q = Q[i];
        fprintf(fsol,"%.10lf", q);
        if (i == L-1) fprintf(fsol,"\n");
        else fprintf(fsol,",");
    }
    // print phases patterns
    for (int k = 0; k < K; k++) {
        fprintf(fsol,"B%d=", k+1);
        for (int i = 0; i < L; i++) {
            fprintf(fsol,"%.10lf", phases[k].xrd[i]);
            if (i == L-1) fprintf(fsol,"\n");
            else fprintf(fsol,",");
        }
    }

    // print the proportion each phase accounts for at every sample point
    fprintf(fsol,"\n// Phase concentrations at each sample\n");
    for (int  n = 0; n < N; n++) {
        double sum = 0.0;
        for (int k = 0; k < K; k++) {
            sum += phase_mixtures[k][n].proportion;
        }
        sum = 1.0;
        fprintf(fsol,"C%d=", n+1);
        for (int k = 0; k < K; k++) {
            double c;
            if (sum > 0.0) c = phase_mixtures[k][n].proportion/sum;
            else c = 0.0;

            if (whetherCC == true) if (CC(k, n) == true) if (c < eps) c = eps;

            fprintf(fsol,"%.12lf", c);
            if (k == K-1) fprintf(fsol,"\n");
            else fprintf(fsol,",");
        }
    }
    
    // print the reconstructed signals of each phase at each sample point
    fprintf(fsol,"\n// Per-phase model for each sample\n");
    double logxrd[len+1];
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < K; k++) {
            fprintf(fsol,"R%d_%d=", n+1, k+1); // Ri_j denotes the signals of phase j at sample point i
            for (int l = 0; l < len; l++) logxrd[l] = 0.0;

            for (int m = 0; m < M; m++) {
                for (int l = 0; l < len; l++) {
                    int tl = l-m;
                    if (tl < 0) continue;
                    if (tl > len-1) continue;
                    logxrd[l] += W(tl, k) * H[m](k,n);
                }
            }
            vector<double> X(Qlogsteps.begin(), Qlogsteps.end());
            vector<double> Y;
            for (int l = 0; l < len; l++) Y.push_back(logxrd[l]);

//            Spline s(X, Y);
            tk::spline s;
            s.set_points(X,Y);

            vector<double> rlist;
            for (int j = 0; j < L; j++) {
                rlist.push_back(s(Q[j]));
            }
            for (int l = 0; l < L; l++) {
                fprintf(fsol,"%.10lf", rlist[l]);
                if (l == L-1) fprintf(fsol,"\n");
                else fprintf(fsol,",");
            }
 
        }
    }

    // shifts
    fprintf(fsol,"\n// Per-phase shift for each sample\n");
    for (int n = 0; n < N; n++) {
        fprintf(fsol,"S%d=", n+1);
        for (int k = 0; k < K; k++) {
            fprintf(fsol,"%.10lf", phase_mixtures[k][n].shift);
            if (k == K-1) fprintf(fsol,"\n");
            else fprintf(fsol,",");
        }
    }

    // H tensor
    fprintf(fsol,"\n// H tensor\n");
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < K; k++) {
            fprintf(fsol,"H[*](%d,%d)=", k+1, n+1);
            for (int m = 0; m < M; m++) {
                fprintf(fsol,"%.10lf", H[m](k,n));
                if (m != M-1) fprintf(fsol,",");
                else fprintf(fsol,"\n");
            }
        }
    }

	// Extended Gibbs
    fprintf(fstats,"\n// Extended Gibbs\n");
    for (int n = 0; n < this->param.shiftsNo.size(); n++)
        fprintf(fstats,"EG%d=%d\n", n+1, 3-this->param.shiftsNo[n]);
    
    // sample contribution
    fprintf(fstats,"\n// Per-sample contribution (L1 loss)\n");
    fprintf(fstats,"L=");
	sampleLoss = new double[N+1];
	sampleTot = new double[N+1];

    for (int n = 0; n < N; n++) {
        double loss = 0.0;
		double tot = 0.0;
        for (int l = 0; l < len; l++) {
            double d = Dlog(l,n), dr = DRec(l,n);
//            loss += d*log((d+eps)/(dr+eps))-d+dr;
            loss += abs(d-dr);
			tot += abs(d);
        }
        fprintf(fstats,"%lf", loss);
		sampleLoss[n] = loss;
		sampleTot[n] = tot;
        if (n == N-1) fprintf(fstats,"\n");
        else fprintf(fstats,",");
    }

	fprintf(fstats,"L_proportion=");
    for (int n = 0; n < N; n++) {
        fprintf(fstats,"%.10lf", sampleLoss[n]/sampleTot[n]);
        if (n == N-1) fprintf(fstats,"\n");
        else fprintf(fstats,",");
    }
    
    // sample contribution
    fprintf(fstats,"\n// Per-sample contribution (L1 loss)\n");
    fprintf(fstats,"L1=");
    sampleLoss = new double[N+1];
    sampleTot = new double[N+1];

    for (int n = 0; n < N; n++) {
        double loss = 0.0;
        double tot = 0.0;
        for (int l = 0; l < len; l++) {
            double d = Dlog(l,n), dr = DRec(l,n);
//            loss += d*log((d+eps)/(dr+eps))-d+dr;
            loss += abs(d-dr);
            tot += abs(d);
        }
        fprintf(fstats,"%.10lf", loss);
        sampleLoss[n] = loss;
        sampleTot[n] = tot;
        if (n == N-1) fprintf(fstats,"\n");
        else fprintf(fstats,",");
    }    

	delete [] sampleLoss;
	delete [] sampleTot;


    // sample contribution
    fprintf(fstats,"\n// Per-sample contribution (L2 loss)\n");
    fprintf(fstats,"L2=");
    sampleLoss = new double[N+1];
    sampleTot = new double[N+1];

    for (int n = 0; n < N; n++) {
        double loss = 0.0;
        double tot = 0.0;
        for (int l = 0; l < len; l++) {
            double d = Dlog(l,n), dr = DRec(l,n);
//            loss += d*log((d+eps)/(dr+eps))-d+dr;
            loss += (d-dr)*(d-dr);
            tot += d*d;
        }
        fprintf(fstats,"%.10lf", loss);
        sampleLoss[n] = loss;
        sampleTot[n] = tot;
        if (n == N-1) fprintf(fstats,"\n");
        else fprintf(fstats,",");
    }

    delete [] sampleLoss;
    delete [] sampleTot;    

    
    // sample contribution
    fprintf(fstats,"\n// Per-sample contribution (KL div)\n");
    fprintf(fstats,"KL=");
    sampleLoss = new double[N+1];
    sampleTot = new double[N+1];

    for (int n = 0; n < N; n++) {
        double loss = 0.0;
        double tot = 0.0;
        for (int l = 0; l < len; l++) {
            double d = Dlog(l,n), dr = DRec(l,n);
            loss += d*log((d+eps)/(dr+eps))-d+dr;
        }
        fprintf(fstats,"%.10lf", loss);
        sampleLoss[n] = loss;
        sampleTot[n] = tot;
        if (n == N-1) fprintf(fstats,"\n");
        else fprintf(fstats,",");
    }

    delete [] sampleLoss;
    delete [] sampleTot;

    //fclose(stdout);
    //return;

    fprintf(fstats,"\n\n// stick matching\n\n");

    fclose(fstats);

    freopen(filename, "a", stdout);

    // invoke GaussianCompare to match phases to icdds
    GaussianCompare gc(this->param.MatchShiftTol, this->param.MatchSigmaTol);
    gc.compare(sticksfile, W, H, K, Q, Qlogsteps);

    fclose(stdout);
    
    fstats = fopen(filename, "a");

    now_time = time(NULL);
    fprintf(fstats,"\n\n\ntime cost: %d\n", now_time-init_time);
    
    fclose(fsol);
    fclose(fstats);
}
