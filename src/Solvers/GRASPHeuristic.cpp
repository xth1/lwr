#include "GRASPHeuristic.h"
#define TRACE(x...) x
#define PRINT(x...) TRACE(printf(x));fflush(stdout)
#define WATCH(x) TRACE(cout << #x" = " << x << "\n")

Result GRASPHeuristic::sort(Permutation &perm, ll upper_bound)
{
  int cost_function = this->param["COST_FUNCTION"];
    Result result;
    result.totalCost = 0;
    perm = perm.toUnsigned();
    WATCH(perm);

    while (perm.isSorted(1, perm.size()) == false) {
        vector<Reversal> revs;
        vector<Strip> strips = perm.getAllStrips();

        /* Creates reversals list where the reversals points (begin or end) are
         *the begin of a breapoint or a breakpoint's internal neighbor point.
         *Each breakpoint [a,b] has at most two internal neighbor points:
         * a + 1 and b - 1.
         *It prevents  generating unitary and repeated reversals.
        */
        for (size_t i = 0; i < strips.size(); i++) {
            Strip si = strips[i];
            /* Reversal [si.begin + 1, si.end - 1] */
            if (si.begin + 1 < si.end - 1) {
                revs.push_back(Reversal(
                                   si.begin + 1,
                                   si.end - 1,
				   0,
				   cost_function
                               ));
            }
            for (size_t j = i + 1; j < strips.size(); j++) {
                Strip sj = strips[j];
                /* Reversals starting at si.begin */
                //1st
                if (sj.begin + 1 < sj.end) {
                    revs.push_back(Reversal(
                                       si.begin,
                                       sj.begin + 1,
				       0,
				       cost_function
                                   ));
                }
                //2nd
                if (sj.end - 1 > sj.begin && sj.end - 1 != sj.begin + 1) {
                    revs.push_back(Reversal(
                                       si.begin,
                                       sj.end - 1,
				       0,
				       cost_function
                                   ));
                }
                //3rd
                revs.push_back(Reversal(
                                   si.begin,
                                   sj.end,
				   0,
				   cost_function
                               ));
                /* Reversals starting at si.begin + 1 */
                if (si.begin + 1 < si.end) {
                    //4th
                    if (sj.begin + 1 < sj.end) {
                        revs.push_back(Reversal(
                                           si.begin + 1,
                                           sj.begin + 1,
					   0,
					   cost_function
                                       ));
                    }
                    //5th
                    if (sj.end - 1 > sj.begin && sj.end - 1 != sj.begin + 1) {
                        revs.push_back(Reversal(
                                           si.begin + 1,
                                           sj.end - 1,
					   0,
					   cost_function
                                       ));
                    }
                    //6th
                    revs.push_back(Reversal(
                                       si.begin + 1,
                                       sj.end, 
				       0,
				       cost_function
                                   ));
                }
                /* Reversals starting at si.end - 1 */
                if (si.end - 1 > si.begin && si.end - 1 != si.begin + 1) {
                    //7th
                    if (sj.begin + 1 < sj.end) {
                        revs.push_back(Reversal(
                                           si.end - 1,
                                           sj.begin + 1, 
					   0,
					   cost_function
                                       ));
                    }
                    //8th
                    if (sj.end - 1 > sj.begin && sj.end - 1 != sj.begin + 1) {
                        revs.push_back(Reversal(
                                           si.end - 1,
                                           sj.end - 1,
					   0,
					   cost_function
                                       ));
                    }
                    //9th
                    revs.push_back(Reversal(
                                       si.end - 1,
                                       sj.end,
				       0,
				       cost_function
                                   ));
                }
            }
        }

        /* Create the nops vector (just to preserve the API) */
        vector<ll> nops(revs.size(), 1);

        int best_rev_id = -1;
        for (size_t k = 0; k < hcr_seq.size(); k++) {
            best_rev_id = this->hcr_seq[k]->choose(perm, revs, nops);
            if (best_rev_id != -1)
                break;
        }

        if (best_rev_id == -1) {
            cerr<<"Error ar GRASPHeuristic: reversal not found"<<endl;
            exit(1);
        }

        Reversal best_rev = revs[best_rev_id];
        WATCH(best_rev);
        WATCH(perm);
        result.totalCost += best_rev.cost;
        result.reversals.push_back(best_rev);

        if (result.totalCost > upper_bound) {
            PRINT("OUT: %lld %lld\n", result.totalCost, upper_bound);
            result.totalCost = INF;
            return result;
        }

        perm.reverse(best_rev);
    }
    return result;
}
