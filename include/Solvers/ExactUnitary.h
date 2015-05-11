#ifndef EXACT_UNITARY_HEADER
#define EXACT_UNITARY_HEADER

/* Author: Thiago da Silva Arruda <thiago.xth1@gmail.com>
 * This code is a implementation of algorithm to compute
 * the exact reversal distance (unitary model) of signed permutations.
 * This algorithm was proposed by Bergeron and other authors in the paper
 * "Reversal Distance without Hurdles and Fortresses", 2004.
 */
#include "HeaderDefault.h"
#include "GRIMMInterface.h"
#include "Solver.h"
#include "HeuristicChooseReversal.h"
/* Structure to a store a component [begin..end].
 * direct means the following:
 * -true, direct component: (perm[begin], perm[end])
 * -false, reversed component: (-perm[begin, -perm[end])
 * */
struct Component {
    int begin;
    int end;
    bool direct;
    bool oriented;
    Component(int begin, int end, bool direct, bool oriented)
    {
        this->begin = begin;
        this->end = end;
        this->direct = direct;
        this->oriented = oriented;
    }
};
typedef struct Component Component;

/*
enum nodes {
  ROUND,
  SQUARE
};

typedef nodes type_node;
*/

/* Tree data struct's:
 * root: have index 0
 * square nodes: have id_component equals SQUERE (-1)
 */

#define SQUARE -1
struct Node {
    int id;
    int id_parent;
    int id_component;
    vector<int> adj;
    Node(int id, int id_parent, int id_component)
    {
        this->id = id;
        this->id_parent = id_parent;
        this->id_component = id_component;
    }
};
typedef struct Node Node;

typedef vector<Node> Tree;

// other way
//typedef vector<int,vector<int> > Tree;

class ExactUnitary : public Solver
{

private:

    vector<bool> visited;
    parameter param;

    /* Sequence of Heuristics for choose reversals:
     * Used to choose a
     *
     */
    vector<HeuristicChooseReversal *> hcr_seq;

public:

    ExactUnitary() { }
    ExactUnitary(
        parameter param,
        vector<HeuristicChooseReversal *> hcr_seq
    );

    /* Given a (extended) permutation  perm, this method build a array M.
     * Where M[i] is the nearest element of perm that preceds perm[i]
     * and is greather than perm[i]*/
    void makeNearestGreatherArray(Permutation &perm, vector<ll> &M);

    /* Given a (extended) permutation  perm, this method build a array M.
     * Where M[i] is the nearest element of perm that preceds perm[i]
     * and is smaller than perm[i]*/
    void makeNearestSmallerArray(Permutation &perm, vector<ll> &m);

    /* Given a (extended) permutation  perm, this method find all
     * components of perm.
     * This algorithm was proposed by Bergeron:
     * =>Direct Components
     * The interval (s..i) is a valid direct component if and only if:
     *  - both perm.sign(i) and perm.sign(s) are positive
     *  - all elements between perm.at(s) and perm.at(i) in perm are greather
     *    than perm.at(s) and smaller than perm.at(i), the latter being
     *    equivalent to the simple test M[i] = M[s].
     *  - no element between perm.at(s) and perm.at(i) is missing, i.e,
     *    i - s = perm.at(i) - perm.at(s).
     * =>Reversed Components
     *  Similiar to direct case.
     * */
    void findComponents(Permutation &perm, vector<Component> &components);

    void buildComponentsTree(Permutation &perm, vector<Component> &components, Tree &T);
    void buildBranchTree(vector<Component> &components, Tree &T);

    bool dfs(vector<Component> &components, Tree &T, int v);

    Result sort(Permutation &perm);

    Result sortUnsigned(Permutation &perm);
    Result sortByBergeron(Permutation &perm);


    /* Functions to:
     * A very elementary presentation of the Hannenhalli–Pevzner theory
     */
    void getAllRev(Permutation &perm, vector<Reversal> &revs_orig,
                   vector<Reversal> &revs_reduced,
                   vector<ll> &nops);
    void getOrientedRev(Permutation &perm, vector<Reversal> &revs_orig,
                        vector<Reversal> &revs_reduced,
                        vector<ll> &nops);
    void getGRIMMRev(Permutation &perm, vector<Reversal> &revs_orig,
                     vector<Reversal> &revs_reduced,
                     vector<ll> &nops);
    void getSiepelRev(Permutation &perm, vector<Reversal> &revs_orig,
                     vector<Reversal> &revs_reduced,
                     vector<ll> &nops);

    ll getBreakpointAndNeghborReversals();
    Reversal getUnsignedRev(Permutation &perm);

    void getBreakpointRevs(
        Permutation &perm,
        vector<Reversal> &revs_orig,
        vector<Reversal> &revs_reduced,
        vector<ll> &nops
    );

    ll getNumOrientedPairs(Permutation &perm);
    Reversal getRevWithMaxOrientedPairs(Permutation &perm);
    Reversal mergeHurdles(pii h1, pii h2);
    Reversal cutHurdle(Permutation &perm, pii h1);

    int getNumHurdlesWhenUnsigned(Permutation &perm);
    pii findSimpleHurdle(Permutation &perm, vector<pii> &hurdles);
    void makeConsecutiveMap(Permutation &perm, vector<int> &begin,
                            vector<int> &end);
    void findFramedIntervals(Permutation &perm, vector<pair<int, int> > &hurdles);

    virtual ~ExactUnitary() {};
};

#endif
