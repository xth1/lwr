#ifndef NON_DETERMINISTIC_SOLUTION
#define NON_DETERMINISTIC_SOLUTION 1

#include "Solver.h"
#include "Permutation.h"
#include "Parameter.h"
#include "ExactUnitary.h"
#include "HeuristicChooseReversal.h"
#include "Util.h"
/* This class build non-deterministic solutions for signed LWR-problem.
 * This resuts are intended to used as seed-solution for a multi-start
 * MetaHeuristic
 */
class NonDeterministicSolution: public Solver
{
public:
    Parameter local_param;
    parameter solver_param;
    parameter hcrP;
    parameter statistics;
    ExactUnitary *solver;

    NonDeterministicSolution(
        Parameter local_param,
        parameter solver_param,
        parameter hcrP
    );
    virtual Result sort(Permutation &p);
    vector<Result> sortAllIterations(
        Permutation &p,
        bool output_steps = false,
        ostream &fout = cout
    );
    vector<Result> getKBestResults(Permutation &p);
    virtual ~NonDeterministicSolution() {}
    // void printResult(Result r, ostream &out);
};

#endif