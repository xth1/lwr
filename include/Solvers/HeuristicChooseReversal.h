#ifndef HCR_HEADER
#define HCR_HEADER
#include "HeaderDefault.h"
#include "Permutation.h"
#include "RandomGenerator.h"
#define pdi pair<double, int>
/*
 * Prototype class
 */

struct pdiComparer {
    bool operator() (pdi a, pdi b)
    {
        return a.F > b.F;
    }
};

class HeuristicChooseReversal
{
    /*
     * This method get three parameters:
     * perm: current permutation, in extended way: [0 .. n].
     * revs: array of revesals
     * nops: nops[i] is the number of oriented pairs of
     *  reversal i.
     *
     * Returns index of choosed reversal
     */
protected:
    parameter param;
public:
    virtual int choose(
        Permutation perm,
        vector<Reversal> revs,
        vector<ll> nops,
        parameter local_param = parameter()
    ) = 0;

    /* Gets name of Heuristic */
    virtual string getName();

    void setParameter(string name, double value);
    void setParameter(parameter &param);
    double getParameter(string name);
    bool hasParameter(string name);

    static vector<HeuristicChooseReversal *> heuristics;
    /* Flag to indicate whether the heuristics objects has been built or not*/
    static bool is_heuristics_ready;
    /* List the current implemented heuristics */
    static void BuildHeuristics();

    static void destroyHeuristics();
    /* Find a Heuristic by name */
    static HeuristicChooseReversal * findHeuristic(string name);

    virtual ~HeuristicChooseReversal() {} //do nothing yet
};

class Roulette
{
public:

    static int chooseByRoulette(vector<pdi> &score, int num_rev);
};

/*
 * Abstract class for Grasp Algorithm
 */
class HCRGRASP: public HeuristicChooseReversal
{
public:
    virtual int choose(
        Permutation perm,
        vector<Reversal> revs,
        vector<ll> nops,
        parameter local_param = parameter()
    );
    virtual double computeBenefit(
        const Permutation &perm,
        Reversal rev,
        parameter local_param = parameter()
    ) = 0;
    virtual string getName();
};

/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum length.
 */
class HCRMinLength : public HCRGRASP
{
public:
    virtual int choose(
        Permutation perm,
        vector<Reversal> revs,
        vector<ll> nops,
        parameter local_param = parameter()
    );
    virtual string getName();
    virtual double computeBenefit(
        const Permutation &perm,
        Reversal rev,
        parameter local_param = parameter()
    );
};

/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum length.
 */
class HCRRandom : public HCRGRASP
{
public:
    virtual int choose(
        Permutation perm,
        vector<Reversal> revs,
        vector<ll> nops,
        parameter local_param = parameter()
    );
    virtual string getName();
    virtual double computeBenefit(
        const Permutation &perm,
        Reversal rev,
        parameter local_param = parameter()
    );
};

/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum length.
 */
class HCRMaxLength : public HeuristicChooseReversal
{
public:
    virtual int choose(
        Permutation perm,
        vector<Reversal> revs,
        vector<ll> nops,
        parameter local_param = parameter()
    );
    virtual string getName();
};

/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum Entropy Loss.
 */
class HCRBenefitLoss : public HCRGRASP
{
public:
    double calcLoss(const Permutation &perm, const Reversal rev);
    virtual int choose(
        Permutation perm,
        vector<Reversal> revs,
        vector<ll> nops,
        parameter local_param = parameter()
    );
    virtual double computeBenefit(
        const Permutation &perm,
        Reversal rev,
        parameter local_param = parameter()
    );
    virtual string getName();
};

/*
 * Choose lexicographically minimun reversal
 */
class HCRLexMin : public HeuristicChooseReversal
{
public:

    virtual bool isLower(const Reversal &r1, const Reversal r2);
    virtual int choose(
        Permutation perm,
        vector<Reversal> revs,
        vector<ll> nops,
        parameter local_param = parameter()
    );
    virtual string getName();
};

/*
 * Test
 */
class HCRTest : public HeuristicChooseReversal
{
public:

    virtual double calcBenefit(const Permutation &perm, const Reversal rev);
    virtual double calcHeuristic(Permutation &perm);
    virtual int choose(
        Permutation perm,
        vector<Reversal> revs,
        vector<ll> nops,
        parameter local_param = parameter()
    );
    virtual string getName();
};

class HCRGRASPByEntropy : public HCRBenefitLoss
{
public:

    double calcHeuristic(
        Permutation &perm,
        vector<Reversal> &revs,
        vector<ll> &nops, vector<double> &H
    );
    virtual int choose(
        Permutation perm,
        vector<Reversal> revs,
        vector<ll> nops,
        parameter local_param = parameter()
    );
    virtual string getName();
};

#endif
