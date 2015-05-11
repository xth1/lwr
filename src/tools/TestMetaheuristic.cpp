#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>

#include "ExactUnitary.h"
#include "HeaderDefault.h"
#include "GRIMMInterface.h"
#include "HeuristicChooseReversal.h"
#include "ImproveReversals.h"
#include "GRASPHeuristic.h"
#include "Parameter.h"
#include "NonDeterministicSolution.h"

#define NUM_IT_GRIMM_UNSIGNED 1000

#define psd pair<string, double>

/* Return the parameters for unsigned permutation to be sorted
 * by MessReduction Solver (That is a meta-heuristic) configured to use
 * GRASP heuristic to choose reversals.
 */

parameter getUnsignedGraspParameters()
{
    /* Permutation parameters */
    parameter param_perm;
    param_perm["reversal"] = true;
    param_perm["lwr"] = true;
    param_perm["signed"] = false;
    param_perm["MessType"] = 1;
    param_perm["extended"] = true;

    return param_perm;
}
/* Return default parameters for signed permutations */
parameter getSignedParameters()
{
    parameter param_perm;
    param_perm["signed"] = true;
    param_perm["lwr"] = true;
    param_perm["reversal"] = true;
    param_perm["MessType"] = 1;
    return param_perm;
}

parameter getDefaultParameters() {
    parameter param_b;

    // Inserido por Ulisses
    param_b["IMP.COST_FUNCTION"] = 0;
    param_b["HCR.COST_FUNCTION"] = 0;
    param_b["BUILD_SOLUTION.COST_FUNCTION"] = 0;

    param_b["IMP.NUM_ROUND"] = 2000;
    param_b["IMP.NUM_IT"] = 1;
    param_b["IMP.ROLETTE_PERCENT_DISCARD"] = 25;
    param_b["IMP.USE_MINUS_MIN"] = 1;
    param_b["IMP.WINDOW_ROLETTE"] = 1;
    param_b["IMP.WINDOW_SZ_BEGIN"] = 5;
    param_b["IMP.WINDOW_SZ_INC"] = 1;
    param_b["IMP.WINDOW_SZ_END"] = 20;
    param_b["IMP.DYNAMIC_SIZE"] = 1;
    param_b["HCR.GRASP"] = 1;
    param_b["HCR.NUM_SELECTED_REV"] = 5;
    param_b["HCR.SQUARE_OF_SCORE"] = 1;
    param_b["BUILD_SOLUTION.GRIMMRevs"] = 1;
    param_b["IMP.EXACT_UNITATY_LS"] = 1;
    param_b["HCR.ENTROPY"] = 1;

    return param_b;
}

vector<string> split(string s,char sep)
{
    vector<string> vs;
    size_t i;
    int ini = 0;
    for (i = 0; i < s.size(); i++) {
        if (s[i] == sep) {
            vs.push_back(s.substr(ini,i+1));
            ini = i+1;
        }
    }
    vs.push_back(s.substr(ini,s.size()));
    return vs;
}

vector<ll> read_seq(char *s)
{

    vector<string> vs = Util::split(string(s), ',');

    vector<ll> v;
//  cout<<"seq : "<<endl;
    for (size_t i = 0; i < vs.size(); i++) {
        ll x;
        sscanf(vs[i].c_str(), "%lld", &x);
        v.push_back(x);
    }

    return v;
}

unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;
    a=a-c;
    a=a^(c >> 13);
    b=b-c;
    b=b-a;
    b=b^(a << 8);
    c=c-a;
    c=c-b;
    c=c^(b >> 13);
    a=a-b;
    a=a-c;
    a=a^(c >> 12);
    b=b-c;
    b=b-a;
    b=b^(a << 16);
    c=c-a;
    c=c-b;
    c=c^(b >> 5);
    a=a-b;
    a=a-c;
    a=a^(c >> 3);
    b=b-c;
    b=b-a;
    b=b^(a << 10);
    c=c-a;
    c=c-b;
    c=c^(b >> 15);
    return c;
}

pair<string,string> getParameterPair(string s)
{
    vector<string> vs;
    /* Get the type of parameter */
    vs = Util::split(s, '.');
    assert(vs.size() == 2);
    return make_pair(vs[0], vs[1]);
}

typedef unsigned long long timestamp_t;

/* Get passed time in microseconds */
static timestamp_t get_timestamp ()
{
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t) now.tv_sec * 1000000;
}

/* Get the interval time (in seconds) */

static double get_interval_time(timestamp_t t0, timestamp_t t1)
{
    double time_sec = (t1 - t0) / 1000000.0L;
    return time_sec;
}

vector<HeuristicChooseReversal *> build_hcr_seq(
    string UNUSED(intervalMethod),
    parameter &improveP,
    parameter &hcrP,
    parameter &exactP,
    string heuristicReversalsName
)
{
    /* Create Solvers */
    vector<HeuristicChooseReversal *> hcr_seq;
    HeuristicChooseReversal *hcr;

    string heuristicName;
    // Set the unsigned parater for the Exact unitary algorithm
    exactP["UNSIGNED"] = improveP["UNSIGNED"];

    // Push the heuristic
    hcr = HeuristicChooseReversal::findHeuristic(heuristicReversalsName);
    hcr->setParameter(hcrP);
    hcr_seq.push_back(hcr);
    // Push the alternative heuristic
    hcr = HeuristicChooseReversal::findHeuristic("MinLength");
    hcr_seq.push_back(hcr);

    hcr->setParameter(hcrP);

    return hcr_seq;
}


void cleanResources()
{
    HeuristicChooseReversal::destroyHeuristics();
}

Result tryImprove(
    Result &seed_result,
    ImproveReversals * ir,
    Permutation perm_orig,
    int cost_function
)
{

    // show result
    bool valid = Solver::validateResult(perm_orig, seed_result, cout);
    assert(valid == true);

    Result r_improved = ir->sort(perm_orig, seed_result);

    r_improved = Solver::calcTotalCost(r_improved, cost_function);

    return r_improved;
}

int main(int argc, char *argv[])
{
    vector<HeuristicChooseReversal *> hcr_seq;

    /* Set random seed */
    unsigned long seed = mix(clock(), time(NULL), getpid());

    /* Set random seed */
    srand(seed);

    /* Exact Unitaty Parameters */
    ofstream fout;

    /* Time results */
    timestamp_t time_begin;
    timestamp_t time_end;
    time_begin = get_timestamp();

    /* Log result file */
    ostream *logOut = NULL;
    ostream *logTime = NULL;

    /* Check essential parameters */
    // TODO: start to use parameter names instead of parameter sequence, ex: -p 1,2,4...
    if (argc < 7) {
        printf("Usage:\n");
        printf("lwr_test_metaheuristic permutation number_of_iteration inversion_limit frame_limit is_signed cost_function {-r=[[...]]}.\n\n");
        printf("Parameters description:\n");
        printf(" 1. permutation: the input permutation.\n\n");
        printf(" 2. number_of_iteration: number of iterations that will\n");
            printf("be performed by the Metaheuristic.\n\n");
        printf(" 3. inversion_limit: limit (upper bound) for the number of inversions\n");
            printf("that will be selected at each iteration of the heuristic for building solutions.\n\n");
        printf(" 4. frame_limit: percentege (in the interval (0,100]) of frames that \n");
            printf("will be selected at each iteration of the Metaheuristic.\n\n");
        printf(" 5. is_signed: 1 if the input permutation is signed, 0 otherwise.\n\n");
        printf(" 6. cost_function: 1 for linear function, 2 for quadratic function and 3 logaritmic function.\n\n");
        printf(" 7. -r=[[...]]: Optional parameter to give a seed solution as a input parameter\n");

        return 0;
    }

    /* Seed result */
    Result seed_result;

    /* Permutation parameters */
    parameter perm_param;

    /* Exact Unitary parameters */
    parameter exactP;

    /* Improve parameters */

    parameter improveP;

    /* ChooseHeuristic parameters*/
    parameter hcrP;

    /* list of solvers for choose reversals */
    vector<Solver *> solver_hcr;

    string intervalMethod = "WINDOW";
    string heuristicReversalsName = "BenefitLoss";

    /* -------------------Read the parameters  From terminal----------------- */

    int p_read = 1;

    // Read the permutation's sequence
    vector<ll> perm_seq = read_seq(argv[p_read++]);

    int num_it;
    sscanf(argv[p_read++], "%d", &num_it);

    int inv_limit;
    sscanf(argv[p_read++], "%d", &inv_limit);

    double frame_limit;
    sscanf(argv[p_read++], "%lf", &frame_limit);

    int is_signed;
    sscanf(argv[p_read++], "%d", &is_signed);

    int cost_function;
    sscanf(argv[p_read++], "%d", &cost_function);

    /* Construct the permutation */
    if (is_signed) {
        perm_param = getSignedParameters();
    } else {
        perm_param = getUnsignedGraspParameters();
    }
    Permutation perm = Permutation(perm_param, perm_seq);
    // makes a copy from the original input permutation and set the normal mode
    Permutation perm_orig = perm;

    // Look for a input solution
    Result input_seed_result;
    bool has_input_seed_result = false;
    if (p_read < argc) {
        string solution = string(argv[p_read++]);
        if (solution.substr(0, 2) == "-r") {
            // in this caseheuristicName is supposed to be like -r=[[...]]
            string str_result = solution.substr(3);
	    // Alterado por Ulisses
            input_seed_result = Util::getResultFromStr(str_result.c_str(), cost_function);
            has_input_seed_result = true;

            Permutation perm_copy = perm;
            bool valid = Solver::validateResult(perm_copy, input_seed_result, cout);
            if (!valid) {
                cout<<"Invalid input solution."<<endl;
                exit(1);
            }
        }
    }


    /* -------------Check if the permutation is already sorted--------------- */
    if (perm.isSorted(1, perm.size())) {
        Result::printResult(seed_result, cout);
        cout<<endl;
        return 0;
    }

    /* ----------------------- Set Default Parameters ----------------------- */
 
    parameter default_parameters = getDefaultParameters();

    /* Set the output for results logging */
    if (p_read < argc) {
        logOut = new ofstream(argv[p_read++], std::ofstream::out|ofstream::app);
    }
    if (p_read < argc) {
        logTime = new ofstream(argv[p_read++], std::ofstream::out|ofstream::app);
    }

    /* Look for Parameters */
    string type;
    parameter::iterator it;
    for (it = default_parameters.begin(); it != default_parameters.end(); it++) {
        pair<string, string> p = getParameterPair(it->F);
        string type = p.F;
        string name = p.S;
        double value = it->S;
        // get Num selected rev
        if (type == "HCR") {
            hcrP[name] = value;
        } else if (type == "IMP") {
            improveP[name] = value;
        } else if (type == "BUILD_SOLUTION") {
            exactP[name] = value;
        }
        p_read++;
    }
    /* Look for the unsigned mode */
    improveP["UNSIGNED"] = !is_signed;

    /* Create log file */
    //ofstream log_out("log_test_signed.txt",ios_base::app);

    /* ----------------------Set Parameters from terminal---------------------- */
    hcrP["COST_FUNCTION"]  = cost_function;
    improveP["COST_FUNCTION"]  = cost_function;
    exactP["COST_FUNCTION"]  = cost_function;

    improveP["NUM_ROUND"] = num_it;
    exactP["NUM_SELECTED_REV"] = inv_limit;

    if (frame_limit <= 0 || frame_limit > 100) {
        cerr<<"Error: frame_limit must be in the interval (0, 100]"<<endl;
        exit(1);
    }
    improveP["ROLETTE_PERCENT_DISCARD"] = (100.0 - frame_limit);


    /*---------------- Build the Choosen reversal Heuristic -------------------*/
    hcr_seq = build_hcr_seq(
      intervalMethod,
      improveP,
      hcrP,
      exactP,
      heuristicReversalsName
    );
    if (improveP["EXACT_UNITATY_LS"]) {
        solver_hcr.push_back(new ExactUnitary(exactP, hcr_seq));
    } else if (improveP["GRASP_HEURISTIC_LS"]) {
        solver_hcr.push_back(new GRASPHeuristic(exactP, hcr_seq));
    }


    /*---------------------- Build The seed solution -------------------------*/
    if (has_input_seed_result) {
        seed_result = input_seed_result;
    } else if (is_signed){
      // Alterado por Ulisses
      seed_result = sortByGRIMM(perm, true, cost_function);
    } else {
      seed_result = sortByGRIMMUnsigned(perm, NUM_IT_GRIMM_UNSIGNED, true, cost_function);
    }
        

    /*--------------------- Build Improvement Metaheuristic ------------------*/
    ImproveReversals * ir = new ImproveReversals(
        solver_hcr,
        intervalMethod,
        improveP,
        logOut
    );

    /*-------------------------- Try improve the result ----------------------*/
    Result r_improved = tryImprove(
        seed_result,
        ir,
        perm_orig,
	cost_function
     );


    /*----------------------- Validate result --------------------------------*/
    bool valid = Solver::validateResult(perm_orig, r_improved, cout);

    Result seedR = seed_result;
    //seedR.reversals = seed_result;
    seedR = Solver::calcTotalCost(seedR, cost_function);


    /*------------------------- Print the result -----------------------------*/
    if (valid) {
        Result::printResult(r_improved, cout);
        cout<<endl;
    } else {
        cout<<"ERROR";
        return 1;
    }

    /*------------------------- Print time -----------------------------------*/
    time_end = get_timestamp();
    double time_sec = get_interval_time(time_begin, time_end);

    if (logTime != NULL)
        *(logTime)<<time_sec<<endl;

    /* Clean all resources */
    for (size_t i = 0; i < solver_hcr.size(); i++) {
        delete solver_hcr[i];
    }
    solver_hcr.clear();
    if (logOut != NULL)
        delete logOut;
    delete ir;
    cleanResources();

    return valid;
}
