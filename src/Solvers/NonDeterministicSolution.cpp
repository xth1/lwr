#include "NonDeterministicSolution.h"
#define TRACE(x...)
#define PRINT(x...) TRACE(printf(x));fflush(stdout)
#define WATCH(x) TRACE(cout << #x" = " << x << "\n")


NonDeterministicSolution::NonDeterministicSolution(
    Parameter local_param,
    parameter solver_param,
    parameter hcrP
)
{
    this->local_param = local_param;
    this->solver_param = solver_param;
    this->hcrP = hcrP;

    // heuristic to select reversals in the Bergeron algorithm
    string hc_name;
    if (local_param.hasString("HEURISTIC_CHOOSE_REVERSAL")) {
        hc_name = local_param.getString("HEURISTIC_CHOOSE_REVERSAL");
    } else {
        cerr<<"Heuristic to choose reversal not found: "<<hc_name<<endl;
        exit(1);
    }

    HeuristicChooseReversal *hc = HeuristicChooseReversal::findHeuristic(
                                      hc_name
                                  );

    WATCH(hc_name);
    WATCH(hc);

    // set the parameters
    hc->setParameter(hcrP);

    vector<HeuristicChooseReversal *> hc_seq;
    hc_seq.push_back(hc);

    // Set statistics
    this->solver_param["STATISTICS"] = 1;

    // create solver object
    this->solver = new ExactUnitary(
        this->solver_param,
        hc_seq
    );
    PRINT("FINISHED");
}

Result NonDeterministicSolution::sort(Permutation &p)
{
    PRINT("sorting...");
    int num_it = (int) this->local_param.getNum("NUM_IT");
    WATCH(num_it);
    // Run the solver for several iterations and pick the min-cost result
    ll min_cost = INF;
    Result min_result;
    for (int k = 0; k < num_it; k++) {
        Result r;
        Permutation p_copy = p;
        WATCH(k);
        r = this->solver->sort(p_copy);
        WATCH(r.totalCost);
        if (r.totalCost < min_cost) {
            min_cost = r.totalCost;
            min_result = r;
        }
    }
    return min_result;
}

vector<Result> NonDeterministicSolution::sortAllIterations(
    Permutation &p,
    bool output_steps,
    ostream &fout
)
{
    PRINT("sorting...");
    int num_it = (int) this->local_param.getNum("NUM_IT");
    bool just_diff_sol = false;
    if (this->local_param.hasNum("JUST_DIFF_SOL")) {
        just_diff_sol = (bool)this->local_param.getNum("JUST_DIFF_SOL");
    }
    WATCH(num_it);
    WATCH(just_diff_sol);

    set<Result, ResultComparer> diff_sols;
    // Run the solver for several iterations and pick the min-cost result
    vector<Result> vr(num_it);

    int k = 0;

    if (output_steps) {
        fout<<num_it<<endl;
    }

    ll min_cost = INF;
    Result min_result;
    while (k < num_it) {
        // time computation
        ll begin_time = Util::getTimestamp();

        Result r;
        Permutation p_copy = p;
        r = this->solver->sort(p_copy);
        if (just_diff_sol) {
            if (diff_sols.find(r) == diff_sols.end()) {
                vr[k++] = r;
                diff_sols.insert(r);
            }
            WATCH(k);
        } else {
            vr[k++] = r;
        }
        WATCH(r.totalCost);
        WATCH(k);
        if (r.totalCost < min_cost) {
            min_cost = r.totalCost;
            min_result = r;
        }
        // Insert the current result into set


        ll end_time = Util::getTimestamp();
        double interval_time = Util::calcIntervalTime(begin_time, end_time);

        if (output_steps) {
            // Print the min cost solution
            //this->printResult(min_result, fout);
            //Result::printResult(min_result, fout);
            //(diff_sols.end())->print(fout);
            fout<<min_cost<<": "<<interval_time<<endl;
        }
    }

    this->statistics = solver->getStatistics();
    parameter::iterator it;
    for (it = this->statistics.begin(); it != this->statistics.end(); it++) {
        cout<<it->F<<": "<<it->S<<endl;
    }

    // Compute statistcs
    /*double avg_num_rev = this->statistics["sum_revs"] /
                         this->statistics["num_choose_execs"];
    double avg_num_rev_opt = this->statistics["sum_opt_revs"] /
                             this->statistics["num_choose_execs"];
    cout<<"avg_rev:"<<avg_num_rev<<endl;
    cout<<"avg_rev_opt:"<<avg_num_rev_opt<<endl;
    */

    return vr;
}

/* Function to sort permutations in non-increasing way */
bool resultCostCMP(const Result &a, const Result &b)
{
    if (a.totalCost < b.totalCost)
        return true;
    return false;
}

/* @return the K best reversals: one to each start... */
vector<Result> NonDeterministicSolution::getKBestResults(Permutation &p)
{
    int K = this->local_param.getNum("NUM_STARTS");
    vector<Result> vr = this->sortAllIterations(p);
    std::sort(vr.begin(), vr.end(), resultCostCMP);
    vector<Result> kvr(K);

    // copy the best K results...
    for (int i = 0; i < K; i++) {
        WATCH(vr[i].totalCost);
        kvr[i] = vr[i];
    }

    return kvr;
}