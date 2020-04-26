## Interleaved Agile Combinatorial Factor Decomposition (IAFD)

This code accompanies the following publications:

Gomes, C. P., Bai, J., Xue, Y., Bjorck, J., Rappazzo, B., Ament, S., Bernstein, R., Suram, S. K., van Dover, R. B., Gregoire, J. M. (2019). **CRYSTAL: a multi-agent AI system for automated mapping of materials' crystal structures**. MRS Communications, 1-9. DOI: 10.1557/mrc.2019.50

Bai, J., Bjorck, J., Xue, Y., Suram, S. K., Gregoire, J., & Gomes, C. (2017). **Relaxation Methods for Constrained Matrix Factorization Problems: Solving the Phase Mapping Problem in Materials Discovery.** Fourteenth International Conference on Integration of Artificial Intelligence and Operations Research Techniques in Constraint Programming (CPAIOR), 104-112. DOI: 10.1007/978-3-319-59776-8_9

Bibtex:

```
@article{Gomes2019,
    author = {Gomes, Carla P. and Bai, Junwen and Xue, Yexiang and Bj{\"{o}}rck, Johan and Rappazzo, Brendan and Ament, Sebastian and Bernstein, Richard and Kong, Shufeng and Suram, Santosh K. and van Dover, R. Bruce and Gregoire, John M.},
    doi = {10.1557/mrc.2019.50},
    issn = {2159-6859},
    journal = {MRS Communications},
    month = {apr},
    pages = {1--9},
    title = {{CRYSTAL: a multi-agent AI system for automated mapping of materials' crystal structures}},
    url = {https://www.cambridge.org/core/product/identifier/S2159685919000508/type/journal{\_}article},
    year = {2019}
}
```

```
@article{Bai2017,
    author = {Bai, Junwen and Bjorck, Johan and Xue, Yexiang and Suram, Santosh K. and Gregoire, John and Gomes, Carla},
    doi = {10.1007/978-3-319-59776-8_9},
    journal = {Fourteenth International Conference on Integration of Artificial Intelligence and Operations Research Techniques in Constraint Programming (CPAIOR)},
    pages = {104--112},
    title = {{Relaxation Methods for Constrained Matrix Factorization Problems: Solving the Phase Mapping Problem in Materials Discovery}},
    url = {http://link.springer.com/10.1007/978-3-319-59776-8{\_}9},
    year = {2017}
}
```



### Introduction and Prerequisites

IAFD is an algorithm for demixing X-ray diffraction (XRD) patterns into individual bases/phases, which correspond to distinct compounds/materials, under physical constraints. IAFD has been deployed in Caltech, Toyota Research Institute, etc.

IAFD has two external library dependencies which must be installed:

1) Armadillo, see: http://arma.sourceforge.net
2) ILOG/CPLEX Optimization Studio: https://www.ibm.com/developerworks/downloads/ws/ilogcplex/

IAFD also uses TCLAP (http://tclap.sourceforge.net/) to process command-line arguments and Tino Kluge's spline implementation (http://kluge.in-chemnitz.de/opensource/spline/). These are header-only includes and are provided in the source distribution.

### Compilation

An example Makefile for GNU make and g++ is included, but will need to be modified to reflect the installation location of CPLEX, and the BLAS library (like openblas) used with Armadillo. Additional changes will be needed for alternate compilers.

After installing all the prerequisites, one should be able to compile the code with the following command:

```bash
make
```

If the compilation fails, try to `make clean` first.

### Usage

```bash
./iafd [-h] [OPTIONS] --inst INSTANCE_FILENAME --m M --k K --sol SOLUTION_FILENAME 
```

Options:

```bash
   --pointCorrect <string>
     correct or further refine a single point without changing anything
     else
   --initH <string>
     the file for H initialization
   --MatchSigmaTol <double>
     sigma boundary when matching with ICSD patterns
   --MatchShiftTol <double>
     shift tolerance when matching with ICSD patterns
   --ConnectAgent
     whether to enforce Connectivity rule or not
   --AlloyAgent
     whether to enforce Alloy rule or not
   --AGrounds <int>
     the rounds of AgileFD-Gibbs loop
   --sticks <string>
     the file containing all the stick patterns
   --AlloyTol <double>
     the tolerance for calculating shifts
   --neighbors <string>
     (required)  the file telling the neighbors of every sample point
   --oneVersion
     whether to enforce one-used-shifted-version constraint
   --slice <string>
     specify the slice constraint
   --Gibbs
     whether to enforce Gibbs phase rule
   --mipgap <double>
     (required)  mipgap for MIP
   --sparsity <double>
     The overall sparsity coefficient
   --stepsize <dobule>
     Initial stepsize. default shift value 0 means the user just want to
     use the std stepsize
   --sampleInit <string>
     The filename containing initialzation from single-phase sample points
   --valueInit <string>
     Initialization file containing seeds for phases and phase freezing
   --rec
     Whether output reconstructed signals
   --humanInput <string>
     Human Input txt file
   --shiftInfo
     Whether you want to create a text file that contains information of
     shifts of each sample point
   --beta <double>
     The weighting coefficient of the sparsity term
   --seed <int>
     random seed for the random number generator.(The default value -1
     means time(0))
   --c <double>
     Related to termination criterion: In one iteration, if
     (old_cost-new_cost)<c*old_cost, then the loop terminates
   --time <int>
     The maximum time(seconds) spent to train the model that you could
     accept
   --m <int>
     (required)  The number of possible different shifts
   --k <int>
     (required)  The number of phases
   --sol <string>
     (required)  The output file name of the solution
   --inst <string>
     (required)  Input instance file
   --addNoise
     whether to add noise or not
   --noiseStd <double>
     define the standard deviation if add noise
   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.
   --version
     Displays version information and exits.
   -h,  --help
     Displays usage information and exits.
```



### Example Command

```bash
./iafd --inst input/Ta-Rh-Pd/Ta-Rh-Pd_inst.txt --m 30 --k 6 --time 10 --sol output/Ta-Rh-Pd_output.txt --c 1e-5 --beta 1.0 --mipgap 0.1 --Gibbs --sparsity 0.1 --neighbors input/Ta-Rh-Pd/Ta-Rh-Pd_edges.txt --AGrounds 3 --MatchShiftTol 0.1 --MatchSigmaTol 2.0 --ConnectAgent --AlloyAgent --AlloyTol 0.003 --sticks input/Ta-Rh-Pd/sticks/sticks.txt
```



### Instance File Format

For more information, see this publication and associated datasets:

Le Bras, R., Bernstein, R., Gregoire, J. M., Suram, S. K., Gomes, C. P., Selman, B., & van Dover, R. B. (2014). A Computational Challenge Problem in Materials Discovery: Synthetic Problem Generator and Real-World Datasets. In Twenty-Eighth International Conference on Artificial Intelligence (AAAI'14).

Example:

    Description=Human readable description of the instance, including origin and any preprocessing 
    UUID=unique-identifier
    Format_Version=1.0
    
    //Number of elements
    M=3
    //Element labels
    Elements=Ta,Rh,Pd
    //Sample count
    N=197
    //Coordinate labels for coordinate systems, e.g. substrate deposition and composition
    Deposition=X,Y
    Composition=Ta,Rh,Pd
    
    // Coordinate values data: lists of length N
    X=-4.7857,-6.81E-05,4.7857,......
    Y=-33.5,-33.5,-33.5,.......
    Ta=0.8072182,0.71555525,0.5822375,......
    Rh=0.10984136,0.11910693,0.1169005,......
    Pd=0.08294043,0.16533777,0.30086198,......
    
    //Q values: data domain: length L
    Q=16.000000,16.150000,16.300000,16.450000,16.600000,......
         
    //Intensity measurements: length L for each of N samples
    I1=3.508038,6.712338,5.096438,......
    ......
    I197=0.000000,18.186840,19.199840,......

### Appendix

#### value Initialization File Format

FLAG: --valueInit

A value initialization file specifies the basis patterns to use for initialization, as well as related configuration options. The format is as follows:

    // Comments can be included as complete lines beginning with two slashes
    
    // Q values corresponding to the basis vectors, in the same units as in the instance file (e.g. nm^-1 or A^-1). 
    // Basis vectors are resampled, so these values do not need to exactly match the ones in the instance file.
    
        Q=1,1.1,1.2,... 
    
    // Basis patterns are only shifted to the right (positive shift) in IAFD, 
    // so initial or frozen basis patterns should be specified as far
    // to the left as expected. The V parameter is a multiplicative shift which 
    // is applied to the Q vector and affects all basis patterns in the same way.
    // The following effectively shifts all basis patterns 1% to the left.
    
        V=0.99
    
    // B1...BK specify the basis patterns you wish to initialize. The indices must be 
    // sequential, but fewer than K can be specified if desired.
    
        B1=0.3123,0.545234,......
        B2=0.4324,0.454243,......
        B3=0.42345,1.42344,......
        ...
    
    // Whether to freeze each basis pattern: 0=seed, 1=freeze
    
        F1=1
        F2=0
        F3=1
        ...
    
    // Lists sample indices in which each phase is allowed to appear. For example, 
    // phase 1 could appear at sample points 1,2,54,65,76......, etc. If S is not
    // specified for an initialized basis, it can appear at any sample by default.
    
        S1=1,2,54,65,76,......
        S2=1,2,4,7,111,......
        S3=67,89,111,......
    ...

#### Sample Initialization File Format

FLAG: --sampleInit

Sample initialization is an alternative to value initialization, where the desired basis 
patterns are taken from samples in the instance file. To initialize phase 1 using sample 
point 123 and phase 2 from sample point 56:

    B1=123
    B2=56

Q should not be specified because the patterns come from the instance file. However, 
V, F, and S are set the same as in a value initialization file.

#### H tensor Initialization File Format

FLAG: --initH

Usually, users could use SOLUTION(s) generated by IAFD as the initialization file for initializing H tensor.

If users want to use an independent initialization file for intializing H tensor, please follow the following format:

```
H[*](1,1)=0.000000,0.000000,0.000000,0.321602,0.357336,0.168312,0.000000,0.000000,0.000000,0.000000
......
H[*](9,1)=0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
H[*](1,2)=0.000000,0.000000,0.000000,0.299825,0.333138,0.110080,0.000000,0.000000,0.000000,0.000000
......
H[*](9,2)=0.000000,0.000000,0.000000,0.056948,0.063276,0.023566,0.000000,0.000000,0.000000,0.000000
......
```

The m-th value in "H[*](k,n)=......" denotes the activation of the m-th shifted version of phase k at sample n.

#### ICDD sticks

FLAG: --sticks

Matching basis patterns with icdd sticks is embedded into the code. But the sticks should be specified in advance. Users can use FLAG "--sticks" to tell the path of the sticks file.

The file should follow the format:

```
Q0=6.467400,8.588200,9.204500,12.101500,12.934900,13.392000,16.587100,17.176300,17.400900,17.765100,...
P0=0.126100,0.200200,0.035000,0.032000,0.102100,0.068100,0.010000,0.006000,0.121100,0.211200,1.00000,...
Q1=17.020800,17.102500,21.050000,21.116100,21.590100,24.770400,24.890200,26.424900,29.183400,30.0547,...
P1=0.203200,0.202200,1.000000,0.983000,0.191200,0.264300,0.329300,0.119100,0.001000,0.084100,0.08410,...
...
...
```

Each pair (Qi,Pi) specifies one ICDD pattern (q values and corresponding intensities).

#### Neighbors

FLAG: --neighbors

This flag is used to specify neighbors:

```
0,1,7,8,9,18
1,0,9,2
2,9,10,3,1
3,2,4,10,11,12,21
4,13,3,12,5,23
5,36,6,13,24,4
6,5,14,15,24,25
7,16,17,28,0,8
......
```

Each line tells the neighbors of one sample. For example, "0,1,7,8,9,18" means sample 0 has neighbors 1,7,8,9,18. Note that all the edges are bidirectional. Thus sample 1 also has one neighbor sample 0.

#### Solution File Format

For more information, see this publication and associated datasets:

Le Bras, R., Bernstein, R., Gregoire, J. M., Suram, S. K., Gomes, C. P., Selman, B., & van Dover, R. B. (2014).  A Computational Challenge Problem in Materials Discovery: Synthetic Problem Generator and Real-World Datasets. In Twenty-Eighth International Conference on Artificial Intelligence (AAAI'14).

Example:

```
  Description=Human-readable description of solution method, parameters, and corresponding instance
  UUID=unique-identifier-of-instance
  Format_Version=1.0

  // Number of phases
  K=5

  //List of solution models: each entry lists the variable prefixes associated with a particular mode
  //[Q,R,C] or [Q,R] is required in order to provide an algorithm-agnostic representation (Q assumed to
  //match instance file if not listed).
  Params=[Q,R,C],[Q,B,C,S],[Q,B,H]

  //Values for the listed parameters at each sample; representations for most parameters can be algorithm-specific
  //Q values typically match the instance file
  Q=......

  //Basis vector representation
  B1=......
  ......
  B5=......

  // Phase concentrations at each sample
  C1=......
  ......
  C197=......

  //Representation of each phase as reconstructed at each sample
  R1_1=......
  ......
  R1_5=......
  ......
  R_197_5=......

  //Per-phase multiplicative (scalar) shift at each sample,
  S1=......
  ......
  S197=......

  //H[*](k,n) Tensor (specific to IAFD - S is computed as a weighted aveage using these as weights)
  H[*](1,1)=......
  ......
  H[*](5,197)=......
```



### Contact me

If you have any questions or suggestions, please feel free to contact jb2467@cornell.edu or rab38@cornell.edu.

