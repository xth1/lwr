#include "Solver.h"
#include "HeuristicChooseReversal.h"

class GRASPHeuristic: public Solver
{
    vector<HeuristicChooseReversal *> hcr_seq;
    parameter param;
public:
    GRASPHeuristic(
        parameter param,
        vector<HeuristicChooseReversal *> hcr_seq
    )
    {
        this->param = param;
        this->hcr_seq = hcr_seq;
    }
    virtual Result sort(Permutation &p, ll upper_bound);
};
