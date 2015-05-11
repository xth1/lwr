#ifndef TESTER_HEADER
#define TESTER_HEADER
#include "HeaderDefault.h"
#include "Solver.h"
#include "Permutation.h"
#include "Util.h"
#include "DivideConquer.h"
#include "InstanceGroup.h"
#include "ExactSolution.h"
#include "MessReduction.h"
#include "SequentialDataAnalyser.h"
#include "Util.h"
#include "Min.h"

#define mss map<string,string>

#define pdd pair<double, double>
#define BIN_DATA_EXACT_SOLUTION "/home/thiago/projetos/LWR_zanoni/lwr/data/all_1_10.dat"

typedef struct {
    vector<Solver *> solvers;
    vector<InstanceGroup> instanceGroups;
    string outputPath;
} TestSet;

typedef struct {
    Permutation p;
    ll cost;
    ll num_rev;
} Line;

typedef struct {
    ll begin;
    ll end;
    ll n;
    string name;
    vector<string> nameSolver;
} MinResult;

enum CSVMode {CREATE, ADD};

class Tester
{
    vector<Solver *> solvers;
    vector<MinResult> minResults;
    map<string,ofstream *> mapLog;
    mss param;
    map<string,Solver *> mapSolvers;
    vector<TestSet> testSets;

    map<string, string> solverLabel;



    //path of the Directory of the test-set file
    string origPath;
    string testFileName;

public:

    Tester(vector<TestSet> testSet,mss param = mss());
    Tester(string testFileName,mss param = mss());
    Tester() { };


    virtual bool test();

    virtual bool testZerOneSeq(string ins_file_name, string out_file_name);

    static void buildCSV(vector<Table> tables, string file_name, CSVMode mode);
    static void buildVerticalCSV(vector<Table> tables, string path);
    static void buildCSVByField(vector<Table> tables, string path, string field);
    static void read_log_line(FILE *fin,int n, Line &line);

    void computeMinFile(ll begin, ll end, ll n,vector<string> path, string out_file);
private:
    string solver_name(string sv, string label);
    Solver *getSolver(string sn, string label, parameter solver_param);
    /*
     * Procces path
     */
    string doPath(string path);
    void addMinLog(vector<string> &outputPath);
    void execReversals(vector<ll> &seq, vector<Reversal> r);
    bool cmpSeq(vector<ll> s1,vector<ll> s2);
    void doLog(Permutation &p, Result &r, ofstream &fout);
    vector<TestSet> readTestFile(string fileName);
    bool test(TestSet &testSet);
    Table testInstanceGroup(InstanceGroup &ig, vector<Solver *> &solvers,SequentialDataAnalyser &sda);
    void evaluateDiameter(Table &table, ll n);
    void initMapLog(vector<Solver *> &solvers, vector<string> &outputPath);

    // To binary sequences...
    vector<string> readInstancesFileNames(string list_file_name);
    string extractFileNameFromPath(string path);

    void generateTableResults(ll begin, ll end, ll n,vector<string> path, vector<string> names, vector<string> fields, string out_path);


};

#endif
