#ifndef IMPROVE_REVERSALS_HEADER
#define IMPROVE_REVERSALS_HEADER 1
#include "HeaderDefault.h"
#include "Permutation.h"
#include "Solver.h"
#include "Util.h"
#include "RandomGenerator.h"

#define NUM_IT  10
#define NUM_ROUND 100
#define MIN_PEAK_LENGTH 3

#define SAMPLING_RATE 50
#define ppll pair<pii, ll>
/* Parameters
 * NUM_ROUND
 * NUM_IT
 * NUM_ANCHOR_POINTS
 */
class ImproveReversals: public Solver
{
    vector<Solver *> heuristics;
    string intervalMethod;

    /* This file Used to output results along iterations */
    ostream *logOut;
    //Parameter param;
public:
    ImproveReversals(
        vector<Solver *> &heuristics,
        string chooseInterval = "GRASP",
        parameter param = parameter(),
        ostream *logOut = &cout
    );

    bool hasOutputLog();
    pii findSubsequenceByGRASP(Result &R, ll min_length, ll max_length);
    pii chooseInterval(Permutation &p, Result &R, int length);
    int tryImprove(Permutation &p, Result &R, pii interval);
    void mergeResults(Result &R, Result &nR, pii itv);
    pii findDenseSubsequence(Result &R);

    /* Kinds of sorting */
    Result sortByGraspInterval(Permutation &p, Result &R);
    Result sortByRandInterval(Permutation &p, Result &R);
    Result sortByEntropyPeaks(Permutation &p, Result &R);

    /* Anchor method */
    void generateAnchorPoints(
        Result &R,
        int anchor_points,
        vector<int> &points
    );
    Result sortByAnchorPoints(Permutation &perm,
                              Result &R,
                              int anchor_points
                             );
    /* Window method */
    pii selectMaxCostInterval(vector<ppll> window);
    pii selectIntervalByRandom(vector<ppll> window, double &prob_chosen);
    pii selectIntervalByRolette(
        vector<ppll> window,
        double &prob_chosen,
        ll &itv_cost
    );
    Result sortByWindow(Permutation &perm, Result &R, int window_size);
    Result sortByValueWindow(Permutation &perm, Result &R, ll threshold);

    virtual Result sort(Permutation &p, Result &result);

    /* Entropy Peak */
    void findEntropyPeaks(Permutation &perm, Result &R, vector<pii> &peaks);
    virtual ~ImproveReversals() {}
};

#endif
