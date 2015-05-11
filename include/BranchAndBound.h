#ifndef BRANCHANDBOUND_HEADER
#define BRANCHANDBOUND_HEADER 1
#include "HeaderDefault.h"
#include "Permutation.h"
#include "MessReduction.h"
#include "ExactSolution.h"

class BranchAndBound
{
private:
    ExactSolution *exactSolution;
    bool hasExactSolution;
public:

    BranchAndBound();
    BranchAndBound(ExactSolution *exactSolution);

    ll sort(Permutation &p);

    ll maxLengthIndependentRev(vector<Reversal> R);

    ll _sort(Permutation &perm, ll &lower, ll &upper, ll &cost,
             Reversal prev_rev);

    ll getLowerBound(Permutation &p);
    ll computeLowerBound(Permutation &p);

};
#endif
