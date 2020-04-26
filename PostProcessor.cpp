#include "PostProcessor.h"

PostProcessor::PostProcessor(int K_, int M_, double mipgap_, double tol_) {
    K = K_;
    M = M_;
    mipgap = mipgap_;
    tol = tol_;
}

PostProcessor::PostProcessor() {
}

PostProcessor::PostProcessor(char* filename,vector<mat> H,int N_, int K_, int M_, double mipgap_, double tol_, double alloyeps, vector<double> Q_, vector<double> Qlogsteps_) {
    K = K_;
    M = M_;
    N = N_;
    mipgap = mipgap_;
    tol = tol_;
    first = true;
    Q = Q_;
    Qlogsteps = Qlogsteps_;
    for (int i = 0; i < N; i++) visitMarks.push_back(-1);
   
    // calls the readnb function to populate neighbors vector
    readNeighbors(filename, neighbors, whethernb);
    
    // create args struct
    for (int i = 0; i < N_; i++) {
        args[i].K = K_;
        args[i].M = M_;
        args[i].sample = i;
        args[i].mipgap = mipgap_;
        args[i].h.clear();
        for (int j = 0; j < M_; j++) {
            colvec colh = zeros<colvec>(K_);
            args[i].h.push_back(colh);
        }
    }
    
    
    // create average shift matrix
    // this 
    avgshifts=zeros<mat>(K, N_);
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < N_; j++) {
            double sum = 0.0, sumw = 0.0;
            for (int m = 0; m < M; m++) {
                sum += H[m](i,j)*Qlogsteps[m]/Qlogsteps[0];
                sumw += H[m](i,j);
            }
            if (sumw > alloyeps) avgshifts(i,j) = (sum)/(sumw);
            else avgshifts(i,j) = 1.0;
        }
    }
    

}


int PostProcessor::nodeViol(int sample){
    // calculates the shift and present phases
    // compares with neighbors and gibbs phase rule
    // returns number of degrees of freedom used

    // calculate number of phases
    int numPhases = 0;
    for (int k = 0; k < K; k++) {
        double tmpsum = 0.0;
        double sum = 0.0;
        for (int m = 0; m < M; m++) {
            tmpsum +=  args[sample].h[m](k);
        }
        if (tmpsum > eps) {
            numPhases++;
        }
    }

    // looks how many phases are shifting rel neighbors
    for(int k=0; k<K; k++){
        double myShift = this->avgshifts(k, sample);
        for(int nb =0; nb<this->neighbors[sample].size(); nb++){
            if(abs(myShift - this->avgshifts(k, this->neighbors[sample][nb]) )> 0.4){
               // numPhases++;
                continue; // if one neighbor is shifting, that phase is shifting
            }
        }
    }

    return numPhases;
}


void PostProcessor::delCC(int k){
    // deletes the smallest connected component
    // finds cc through BFS

    fprintf(stderr, "Writing out CCs \n");

    //vector of nodes that are investigated
    vector <int> visited;
    vector<int> tovisit;
    vector<double> amountCCs;
    vector<vector<int> > allCCs;

    //outer loop that finds all CCs
    for(int sample=0;sample<N; sample++){

        // check if already visited
        if(find(visited.begin(), visited.end(), sample) != visited.end()) {
            continue;
        }

        visited.push_back(sample);
        double amount=0;

        for(int m=0; m<M; m++){
            amount +=  args[sample].h[m](k);
        }
        // if no phase present, skip
        if(amount < 0.01){
            continue;
        }

        fprintf(stderr, "found CC starting with %d, printing nodes ", sample);
        vector<int> nodesCC;
        nodesCC.clear();
        nodesCC.push_back(sample);    


        // do BFS from starting node
        tovisit = this->neighbors[sample];
        while(tovisit.size()>0){

            // pop first element
            int node;
            node = tovisit[0];
            tovisit.erase(tovisit.begin());

            // if element already visited, continue
            if(find(visited.begin(), visited.end(), node) != visited.end()) {
                continue;
            }

            // if not visited earlier, add to visited
            visited.push_back(node);

            // look how much of phase we have
            double tempamount=0;
            for(int m=0; m<M; m++){
                tempamount +=  args[node].h[m](k);
            }

            // if no phase present, skip
            if(tempamount < 0.01){
                continue;
            }

            fprintf(stderr, "%d ", node);
            nodesCC.push_back(node);

            // if element present, add neighbors and amount
            amount += tempamount;
            for(int nb=0; nb<this->neighbors[node].size(); nb++){
                tovisit.push_back(this->neighbors[node][nb]);
            }
        }

        // push back results and print
        amountCCs.push_back(amount);
        allCCs.push_back(nodesCC);
        fprintf(stderr, "size of the CC is: %f \n", amount);        

    }// end find CC loop

    // if only one component, break
    if (amountCCs.size() < 2) return;

    // find smallest component
    int smallestInd =0;
    double smallestVal = 100;

    for(int i=0; i<amountCCs.size(); i++){
        if(smallestVal > amountCCs[i]){
            smallestInd = i;
            smallestVal = amountCCs[i];
        }
    }

    // delete smallest component
    vector<int> toDelete = allCCs[smallestInd];
    for(int node=0; node<toDelete.size(); node++){
        for(int m=0; m<M; m++){
            args[toDelete[node]].h[m](k) = 0;
        }
    }

}//end function

void PostProcessor::readNeighbors(char* filename, vector<vector<int> > &nb, bool &whethernb) {

    if (strlen(filename) == 0) { // users do not use this functionality
        whethernb = false;
        return;
    }
    FILE * file = freopen(filename, "r", stdin);
    if (file == NULL) { // no such file
        whethernb = false;
        return;
    }

    whethernb = true;
    nb.clear();
    char c;

    int sample;
    for (int i = 0; i < this->N; i++) {
        vector<int> onenb;
        onenb.clear();

        scanf("%d", &sample);
        c = getchar();
        while (c == ',') {
            int v;
            scanf("%d", &v);
            onenb.push_back(v);
            c = getchar();
        }
        nb.push_back(onenb);
    }

    fclose(stdin);
}

void PostProcessor::govisit(int sample, double &sumv, vector<int> &visit, int cnt, vector<mat> &H, int k, map<set<int>, int> &combmap, map<int, set<int> > &combdict, vector<int> &mark, set<int> &losers) {

    int N = H[0].n_cols;
    int M = H.size();
    double tmpsum = 0.0;
    for (int m = 0; m < M; m++)
        for (set<int>::iterator it = combdict[k].begin(); it != combdict[k].end(); it++)
            tmpsum += H[m](*it, sample);
    if (tmpsum < eps) return;
    
    visit[sample] = cnt;
    sumv += tmpsum;

    for (vector<int>::iterator it = neighbors[sample].begin(); it != neighbors[sample].end(); it++) {
        if (mark[*it] == k && visit[*it] == -1) {
            govisit(*it, sumv, visit, cnt, H, k, combmap, combdict, mark, losers);
        }
    }

}

// visit single phase region
void PostProcessor::govisitphase(int sample, double &sumv, vector<int> &visit, int cnt, vector<mat> &H, int k) {
    int N = H[0].n_cols;
    int M = H.size();

    // double check whether phase k is activated or not
    double tmpsum = 0.0;
    for (int m = 0; m < M; m++)
        tmpsum += H[m](k, sample);
    if (tmpsum < eps) return;
    
    visit[sample] = cnt;
    sumv += tmpsum;

    for (vector<int>::iterator it = neighbors[sample].begin(); it != neighbors[sample].end(); it++) {
        if (visit[*it] == -1) {
            govisitphase(*it, sumv, visit, cnt, H, k);
        }
    }
}

// check connectivity
bool PostProcessor::checkCC(vector<mat> &H, int N, int K, int M) {
    for (int k = 0; k < K; k++) {
        vector<int> visit;
        for (int n = 0; n < N; n++) visit.push_back(-1);
        
        int visitCnt = 0;
        for (int n = 0; n < N; n++) {
            double sum = 0.0;
            for (int m = 0; m < M; m++) sum += H[m](k,n);
            if (sum < eps) continue;

            if (visit[n] == -1) {
                double sumv = 0.0;
                govisitphase(n, sumv, visit, visitCnt, H, k);
                visitCnt++;
                if (visitCnt > 1) return false;
            }
        }
    }
    return true;
}

// enforce single phase region connectivity
void PostProcessor::connectphase(vector<mat> &H, int N, int K, int M, Mat<int> &CC) {

    CC = zeros<Mat<int> >(K, N);
    set<int> losers;
	vector<int> mark;
	for (int i = 0; i < N; i++) mark.push_back(-1);

    map<set<int>, int> combmap;
    map<int, set<int> > combdict;
	int combcnt=0;


    for (int k = 0; k < K; k++) {
        vector<int> visit;
        for (int n = 0; n < N; n++) visit.push_back(-1);

        vector<double> weights;
        weights.clear();
        int cnt = 0;
        double finalsumv = 0.0;
        int chosencnt = 0;

        for (int n = 0; n < N; n++) {
            double sum = 0.0;
            for (int m = 0; m < M; m++) sum += H[m](k,n);
            if (sum < eps) continue;

            if (visit[n] == -1) {
                double sumv = 0.0;
                govisitphase(n, sumv, visit, cnt, H, k);
                if (sumv > finalsumv) {
                    finalsumv = sumv;
                    chosencnt = cnt;
                }

                cnt++;
            }
        }
        for (int n = 0; n < N; n++) {

            double sum = 0.0;
            for (int m = 0; m < M; m++) sum += H[m](k,n);
            if (sum < eps) continue;

            if (visit[n] != chosencnt) {
                for (int m = 0; m < M; m++) H[m](k,n) = 0.0;
                CC(k,n) = 0;
                losers.insert(n);
            } else {
                CC(k,n) = 1;
            }
        }

    }
}


void PostProcessor::connect(vector<mat> &H, int N, int K, int M, Mat<int> &CC) {

    //CC = zeros<Mat<int> >(K, N); // connected components
    map<set<int>, int> combmap; // map phase set to phase-set index
    map<int, set<int> > combdict; // map phase-set index to phase set
    vector<int> mark; // record the index of each sample
    for (int i = 0; i < N; i++) mark.push_back(-1);

    int combcnt = 0;
    set<int> losers; // empty samples

    for (int n = 0; n < N; n++) {
        set<int> comb;
        for (int k = 0; k < K; k++) {
            double sum = 0.0;
            for (int m = 0; m < M; m++) {
                sum += H[m](k,n);
            }
            if (sum > eps && CC(k, n)) {
                comb.insert(k);
            }
        }
        if (comb.size() == 0) {
            losers.insert(n);
            continue;
        }
        if (combmap.find(comb) == combmap.end()) {
            combmap[comb] = combcnt;
            combdict[combcnt] = comb;
            combcnt++;
        }
        mark[n] = combmap[comb];
    }

//    /*
    fprintf(stderr, "NUMBER of losers: %d. Notice 1st !!!!!\n", losers.size());

    // iterate all the phase sets
    for (int i = 0; i < combcnt; i++) {
        vector<int> visit;
        for (int n = 0; n < N; n++) visit.push_back(-1);
        vector<double> weights;
        weights.clear();
        int cnt = 0;
        double finalsumv = 0.0;
        int chosencnt = 0;

        for (int n = 0; n < N; n++) {
            if (visit[n] == -1 && mark[n] == i) { // this sample is not visited and its phase-set index is the same as the current phase set
                double sumv = 0.0;
                govisit(n, sumv, visit, cnt, H, i, combmap, combdict, mark, losers); // dfs from sample n
                if (sumv > finalsumv) {
                    finalsumv = sumv;
                    chosencnt = cnt;
                }
                cnt++;
            }
        }
        
        for (int n = 0; n < N; n++) {
            if (mark[n] == i) {
                if (visit[n] != chosencnt) losers.insert(n); // this sample is not in the main component. Abandoned!
                else { // main connected component
                    for (set<int>::iterator it = combdict[i].begin(); it != combdict[i].end(); it++) {
                        CC(*it, n) = true;
                    }
                }
            }
        }
    }
//    */

    map<int, int> validnb;

    fprintf(stderr, "NUMBER of losers: %d. Notice 2nd!!!!!\n", losers.size());

    // calculate the number of valid neighbors of each 'loser'
    for (set<int>::iterator loser = losers.begin(); loser != losers.end(); loser++) {
        int validn = 0;
        for (vector<int>::iterator it = neighbors[*loser].begin(); it != neighbors[*loser].end(); it++) {
            if (losers.find(*it) == losers.end()) validn += 1;
        }
        validnb[*loser] = validn;
    }

    // gradually fill 'losers' (empty samples)
    while (losers.size()) {
        int largest = 0, chosenone = -1;

        // find the 'loser' which has the most number of valid neighbors
        for (set<int>::iterator loser = losers.begin(); loser != losers.end(); loser++) {
            if (validnb[*loser] > largest) {
                largest = validnb[*loser];
                chosenone = *loser;
            }
        }
        
        // calculate the frequency of different phase sets among this sample's valid neighbors
        map<int, int> mapcnt;
        for (vector<int>::iterator it = neighbors[chosenone].begin(); it != neighbors[chosenone].end(); it++) {
            if (losers.find(*it) == losers.end()) {
                int type = mark[*it];
                if (mapcnt.find(type) == mapcnt.end()) mapcnt[type] = 1;
                else mapcnt[type] += 1;
            } else {
                validnb[*it] += 1;
            }
        }

        // find the most popular phase set
        int maxtype = -1, maxnum = 0;
        for (map<int, int>::iterator it = mapcnt.begin(); it != mapcnt.end(); it++) {
            if (it->second > maxnum) {
                maxnum = it->second;
                maxtype = it->first;
            }
        }

        // fill this empty sample
        mark[chosenone] = maxtype;
        for (set<int>::iterator it = combdict[maxtype].begin(); it != combdict[maxtype].end(); it++) {
            CC(*it, chosenone) = true;
        }
        losers.erase(chosenone);
    }

}


void PostProcessor::addGibbs(int sample, vector<Coordinate> &xyz, vector<mat> &H, mat& W, mat& D,  vector<int> &shiftsNo, bool onev, bool alloy, double epsilon, double alloyeps, bool whetherCC, Mat<int> &CC) { // use mip to enforce Gibbs phase rule

    // alloyeps denotes how much shift to put in avgshift matrix
    // epsilon denotes how much slack we can have wrt neighbors

    bool permitprint = (sample % 10 == 0);

    vector<vector<int> > nb = neighbors;
    mat avgshifts = this->avgshifts;
   
    IloEnv env; // set up the environment
    int L = W.n_rows; // the length of one XRD pattern


    try {
        
        IloModel model(env);
        IloIntVarArray I(env);
        IloIntVarArray Im(env);
        IloIntVarArray Ishifts(env);

        for (int i = 0; i < K; i++) {
            I.add(IloIntVar(env,0,1)); // K indicator variables indicating whether to select the k-th phase or not
        }
        for (int i = 0; i < M*K; i++) {
            Im.add(IloIntVar(env,0,1));
        }
    

        // decision variable for shifting
        Ishifts.add(IloIntVar(env,0,1));

        IloNumVarArray Hc(env);
        for (int m = 0; m < M; m++)
        for (int k = 0; k < K; k++) {
            if (whetherCC == true) {
                if (CC(k,sample) == 0) {Hc.add(IloNumVar(env, 0.0, 1e-30));continue;}
                else {Hc.add(IloNumVar(env, 1.01/M*eps, IloInfinity));continue;}
            }
            
            
            double tmp = 0.0;
            for (int i = 0; i < M; i++) tmp += H[i](k, sample);
            if (tmp < 1e-8) Hc.add(IloNumVar(env, 0.0, 1e-30));
            else Hc.add(IloNumVar(env, 0.0, IloInfinity)); // K phases, each of which has M shifted versions
            
            //Hc.add(IloNumVar(env, 0.0, IloInfinity));
        }

        IloNumVarArray t(env); // the L1 loss of the differences between reconstructed signals (of length L) and original signals (one column of D matrix)
        for (int i = 0; i < L; i++) t.add(IloNumVar(env, 0.0, IloInfinity));
       

        IloExpr obj(env); // objective function
     
        bool Dneg = false;

        for (int l = 0; l < L; l++) { // t[l] = |reconstructed_signals[l]-original_signals[l]|
                IloExpr e(env); // reconstructed signals
                for (int m = 0; m < M; m++) {
                    for (int k = 0; k < K; k++) {
                        if (l-m < 0) continue;
                        e += W(l-m,k)*Hc[m*K+k];
                    }
                }

                // t[l] = |reconstructed_signals[l]-original_signals[l]|
            
                model.add(t[l]-e>=-D(l,sample));
                model.add(t[l]+e>=D(l,sample));

                obj += t[l];
//                obj += (e-D(l,sample))*(e-D(l,sample));
        }


        // ----------extended Gibbs------------
        IloExpr sumIshifts(env);

        if (alloy) {

            if (permitprint) fprintf(stderr, "epsilon: %lf\n", epsilon);

            vector<int>::iterator it;
            for (int k = 0; k < K; k++) {
                IloExpr sum(env), sumw(env);
                for (int m = 0; m < M; m++) {
                    sum += Hc[m*K+k] * Qlogsteps[m]/Qlogsteps[0];
                    sumw += Hc[m*K+k];
                }

                for (it=nb[sample].begin(); it<nb[sample].end(); it++)
                if (avgshifts(k, *it) > 0.0) {
                    model.add(sum-avgshifts(k, *it)*sumw <= epsilon*sumw+M/2.0*Ishifts[0]+M/2.0*(1-I[k]));
                    model.add(sum-avgshifts(k, *it)*sumw >= -epsilon*sumw-M/2.0*Ishifts[0]-M/2.0*(1-I[k]));
                }
            }
            sumIshifts += Ishifts[0];

        }
        // ---------------end------------------

        if (onev) {
        // ---------------new------------------

            if (permitprint) fprintf(stderr, "enforcing one-shifted-version constraint!!!!\n");
            for (int k = 0; k < K; k++) {
                IloExpr e(env);
                for (int m = 0; m < M; m++) {
                    e += Im[m*K+k];
                }
                model.add(e<=1);
            }
            for (int k = 0; k < K; k++) {
				
                model.add(Im[0*K+k]+Im[1*K+k]>=Hc[0*K+k]);
                model.add(Im[(M-1)*K+k]+Im[(M-2)*K+k]>=Hc[(M-1)*K+k]);
                for (int m = 1; m < M-1; m++) {
                    model.add(Im[(m-1)*K+k]+Im[m*K+k]+Im[(m+1)*K+k]>=Hc[m*K+k]);
                }
				/*
				model.add(Im[0*K+k]>=Hc[0*K+k]);
				for (int m = 1; m < M; m++)
					model.add(Im[(m-1)*K+k]+Im[m*K+k]>=Hc[m*K+k]);
				*/
            }
			
            for (int k = 0; k < K; k++) {
                for (int m = 1; m < M; m++) model.add(0.9*Hc[m*K+k]-Hc[(m-1)*K+k]>=Im[m*K+k]-1);
                for (int m = 0; m < M-1; m++) model.add(0.9*Hc[m*K+k]-Hc[(m+1)*K+k]>=Im[m*K+k]-1);
            }
			

        // ---------------end------------------
        }

        for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) {
                    model.add(I[k]-Hc[m*K+k]>=0); // if I[k] == 0, the coefficient of every shifted version of the k-th phase should be 0. Otherwise, there is no boundaries for the coefficients
                }
        }
        

        IloExpr e(env);
        for (int k = 0; k < K; k++) {
            e += I[k]; // calculate the number of used phases
        }

        // enforce Gibbs phase rule

        double x = xyz[0].v[sample], y = xyz[1].v[sample], z = xyz[2].v[sample];
        
        int numlim = 3;
        
        
        if (abs(x+y+z-1) < eps) {
            if (x < eps) numlim--;
            if (y < eps) numlim--;
            if (z < eps) numlim--;
        }
        
        
        if (alloy) model.add(e+sumIshifts<=numlim);
        else model.add(e<=numlim);

        model.add(e>=1);

        if (alloy) obj += 0.001 * Ishifts[0];
        
        model.add(IloMinimize(env, obj)); // minimize the objective function

        
        // -------------------------------------

        IloCplex cplex(model);

        // -------------------------------------

        cplex.setParam(IloCplex::EpGap, mipgap); // set the mipgap
        cplex.setParam(IloCplex::Threads, 1); // set # of threadd
        cplex.setParam(IloCplex::MIPDisplay,0); // set verbose
        cplex.setParam(IloCplex::SimDisplay,0);
        cplex.setParam(IloCplex::BarDisplay,0);
        cplex.setParam(IloCplex::NetDisplay,0);
        cplex.setParam(IloCplex::Param::ParamDisplay,0);
        

        if (!cplex.solve()) { // start to solve the MIP
            env.error() << "Failed to optimize LP." << endl;
            throw(-1);
        }

//        sumcost += cplex.getObjValue(); // the minimum value of the objective function
        if (permitprint)
            fprintf(stderr, "Solution value = %lf\n", cplex.getObjValue());


        // get values of the variables
        IloNumArray Hvalues(env);
        cplex.getValues(Hvalues, Hc);
        
        vector<int> intvalues;
        vector<int> Imvalues;
        for (int i = 0; i < K; i++) {
            intvalues.push_back(cplex.getValue(I[i])); // indicators
        }


        int cnt = 0;

        for (int k = 0; k < K; k++) {
            double tmpsum = 0.0;
            double sum = 0.0;
            for (int m = 0; m < M; m++) {
                args[sample].h[m](k) = Hvalues[m*K+k];

                if (whetherCC == true) {
                    if (CC(k,sample) == true)
                        if (args[sample].h[m](k) < 1.01/M*eps) args[sample].h[m](k) = 1.0/M*eps;
                }

                sum += Hvalues[m*K+k]*Qlogsteps[m]/Qlogsteps[0];
                tmpsum += Hvalues[m*K+k];
            }
            if (tmpsum > alloyeps) {
                cnt++;
                if (permitprint) {
                    fprintf(stderr, "sample %d -- phase %d: %lf || ", sample, k, sum/tmpsum);
                    for (vector<int>::iterator it=nb[sample].begin(); it!=nb[sample].end(); it++)
                        fprintf(stderr, "%d: %lf(%lf) ", *it, avgshifts(k, *it), avgshifts(k, *it)-sum/tmpsum);
                    fprintf(stderr, "\n");
                }
                this->avgshifts(k, sample) = sum/tmpsum;
            } else {
                if (permitprint) {
                    fprintf(stderr, "phase %d: no\n", k);
                }
                this->avgshifts(k, sample) = 1.0;
            }
        }
        if (permitprint)
            fprintf(stderr, "sample %d: %d\nsumIshifts: %d\n", sample, cnt, (int)cplex.getValue(sumIshifts));
        shiftsNo[sample] = (int)(cplex.getValue(sumIshifts)+0.5);

    }
    catch (IloException &e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();      

}

void PostProcessor::correctPoint(int sample, vector<Coordinate> &xyz, vector<mat> &H, mat& W, mat& D, vector<int> &shiftsNo, bool onev, bool alloy, double epsilon, double alloyeps, set<int> phaseSet) { // use mip to enforce Gibbs phase rule

    // alloyeps denotes how much shift to put in avgshift matrix
    // epsilon denotes how much slack we can have wrt neighbors

    bool permitprint = (sample % 10 == 0);

    vector<vector<int> > nb = neighbors;
    mat avgshifts = this->avgshifts;
   
    IloEnv env; // set up the environment
    int L = W.n_rows; // the length of one XRD pattern


    try {
        
        IloModel model(env);
        IloIntVarArray I(env);
        IloIntVarArray Im(env);
        IloIntVarArray Ishifts(env);

        for (int i = 0; i < K; i++) {
            I.add(IloIntVar(env, 0, 1));
        }
        for (int m = 0; m < M; m++)
            for (int k = 0; k < K; k++)
                Im.add(IloIntVar(env, 0, 1));

        // decision variable for shifting
        Ishifts.add(IloIntVar(env,0,1));

        IloNumVarArray Hc(env);

        //fprintf(stderr, "new sample: %d\n", sample);

        for (int m = 0; m < M; m++)
        for (int k = 0; k < K; k++) {
            if (phaseSet.find(k) != phaseSet.end()) {
                //fprintf(stderr, "k: %d\n", k);
                Hc.add(IloNumVar(env, 1.0/M*eps, IloInfinity));
            }
            else {
                Hc.add(IloNumVar(env, 0.0, 1e-30));
            }
        }

        IloNumVarArray t(env); // the L1 loss of the differences between reconstructed signals (of length L) and original signals (one column of D matrix)
        for (int i = 0; i < L; i++) t.add(IloNumVar(env, 0.0, IloInfinity));
       

        IloExpr obj(env); // objective function
     
        bool Dneg = false;

        for (int l = 0; l < L; l++) { // t[l] = |reconstructed_signals[l]-original_signals[l]|
                IloExpr e(env); // reconstructed signals
                for (int m = 0; m < M; m++) {
                    for (int k = 0; k < K; k++) {
                        if (l-m < 0) continue;
                        e += W(l-m,k)*Hc[m*K+k];
                    }
                }

                // t[l] = |reconstructed_signals[l]-original_signals[l]|
            
                model.add(t[l]-e>=-D(l,sample));
                model.add(t[l]+e>=D(l,sample));

                obj += t[l];
//                obj += (e-D(l,sample))*(e-D(l,sample));
        }


        // ----------extended Gibbs------------
        IloExpr sumIshifts(env);

        if (alloy) {

            if (permitprint) fprintf(stderr, "epsilon: %lf\n", epsilon);

            vector<int>::iterator it;
            for (int k = 0; k < K; k++) {
                IloExpr sum(env), sumw(env);
                for (int m = 0; m < M; m++) {
                    sum += Hc[m*K+k] * Qlogsteps[m]/Qlogsteps[0];
                    sumw += Hc[m*K+k];
                }

                for (it=nb[sample].begin(); it<nb[sample].end(); it++)
                if (avgshifts(k, *it) > 0.0) {
                    model.add(sum-avgshifts(k, *it)*sumw <= epsilon*sumw+M*Ishifts[0]+M*(1-I[k]));
                    model.add(sum-avgshifts(k, *it)*sumw >= -epsilon*sumw-M*Ishifts[0]-M*(1-I[k]));
                }
            }
            sumIshifts += Ishifts[0];

        }
        // ---------------end------------------

        if (onev) {
        // ---------------new------------------

            if (permitprint) fprintf(stderr, "enforcing one-shifted-version constraint!!!!\n");
            for (int k = 0; k < K; k++) {
                IloExpr e(env);
                for (int m = 0; m < M; m++) {
                    e += Im[m*K+k];
                }
                model.add(e<=1);
            }
            for (int k = 0; k < K; k++) {
				
                model.add(Im[0*K+k]+Im[1*K+k]>=Hc[0*K+k]);
                model.add(Im[(M-1)*K+k]+Im[(M-2)*K+k]>=Hc[(M-1)*K+k]);
                for (int m = 1; m < M-1; m++) {
                    model.add(Im[(m-1)*K+k]+Im[m*K+k]+Im[(m+1)*K+k]>=Hc[m*K+k]);
                }
				/*
				model.add(Im[0*K+k]>=Hc[0*K+k]);
				for (int m = 1; m < M; m++)
					model.add(Im[(m-1)*K+k]+Im[m*K+k]>=Hc[m*K+k]);
				*/
            }
			
            for (int k = 0; k < K; k++) {
                for (int m = 1; m < M; m++) model.add(0.9*Hc[m*K+k]-Hc[(m-1)*K+k]>=Im[m*K+k]-1);
                for (int m = 0; m < M-1; m++) model.add(0.9*Hc[m*K+k]-Hc[(m+1)*K+k]>=Im[m*K+k]-1);
            }
			

        // ---------------end------------------
        }

        for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) {
                    model.add(I[k]-Hc[m*K+k]>=0); // if I[k] == 0, the coefficient of every shifted version of the k-th phase should be 0. Otherwise, there is no boundaries for the coefficients
                }
        }
        

        IloExpr e(env);
        for (int k = 0; k < K; k++) {
            e += I[k]; // calculate the number of used phases
        }

        // enforce Gibbs phase rule

        double x = xyz[0].v[sample], y = xyz[1].v[sample], z = xyz[2].v[sample];
        
        int numlim = 3;
        
        
        if (abs(x+y+z-1) < eps) {
            if (x < eps) numlim--;
            if (y < eps) numlim--;
            if (z < eps) numlim--;
        }
        
        
        if (alloy) model.add(e+sumIshifts<=numlim);
        else model.add(e<=numlim);

        model.add(e>=1);

        if (alloy) obj += 0.001 * Ishifts[0];
        
        model.add(IloMinimize(env, obj)); // minimize the objective function

        
        // -------------------------------------

        IloCplex cplex(model);

        // -------------------------------------

        cplex.setParam(IloCplex::EpGap, mipgap); // set the mipgap
        cplex.setParam(IloCplex::Threads, 1); // set # of thread
        cplex.setParam(cplex.MIPDisplay,0); // set verbose
        cplex.setParam(cplex.SimDisplay,0);
        

        if (!cplex.solve()) { // start to solve the MIP
            env.error() << "Failed to optimize LP." << endl;
            throw(-1);
        }

//        sumcost += cplex.getObjValue(); // the minimum value of the objective function
        if (permitprint)
            fprintf(stderr, "Solution value = %lf\n", cplex.getObjValue());


        // get values of the variables
        IloNumArray Hvalues(env);
        cplex.getValues(Hvalues, Hc);
        
        vector<int> intvalues;
        vector<int> Imvalues;
        for (int i = 0; i < K; i++) {
            intvalues.push_back(cplex.getValue(I[i])); // indicators
        }


        int cnt = 0;

        for (int k = 0; k < K; k++) {
            double tmpsum = 0.0;
            double sum = 0.0;
            for (int m = 0; m < M; m++) {
                H[m](k, sample) = Hvalues[m*K+k];

                if (H[m](k, sample) < 1.0/M*eps && phaseSet.find(k) != phaseSet.end()) H[m](k, sample) = 1.0/M*eps;

                sum += Hvalues[m*K+k]*Qlogsteps[m]/Qlogsteps[0];
                tmpsum += Hvalues[m*K+k];
            }
            //fprintf(stderr, "%lf, ", tmpsum); // ~~~~~~~~~~~
            if (tmpsum > alloyeps) {
                cnt++;
                if (permitprint) {
                    fprintf(stderr, "sample %d -- phase %d: %lf || ", sample, k, sum/tmpsum);
                    for (vector<int>::iterator it=nb[sample].begin(); it!=nb[sample].end(); it++)
                        fprintf(stderr, "%d: %lf(%lf) ", *it, avgshifts(k, *it), avgshifts(k, *it)-sum/tmpsum);
                    fprintf(stderr, "\n");
                }
                this->avgshifts(k, sample) = sum/tmpsum;
            } else {
                if (permitprint) {
                    fprintf(stderr, "phase %d: no\n", k);
                }
                this->avgshifts(k, sample) = 1.0;
            }
        }
        //fprintf(stderr, "\n"); // ~~~~~~~~~~~
        if (permitprint)
            fprintf(stderr, "sample %d: %d\nsumIshifts: %d\n", sample, cnt, (int)cplex.getValue(sumIshifts));
        shiftsNo[sample] = (int)(cplex.getValue(sumIshifts)+0.5);

    }
    catch (IloException &e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();      

}


