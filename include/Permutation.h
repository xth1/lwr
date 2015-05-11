#ifndef PERMUTATION_HEADER
#define PERMUTATION_HEADER 1
#include "HeaderDefault.h"
#include "Moviment.h"

#define N_BITS 4 // para representação de permutações em ll
#define MASK 15LL

typedef map<string,double> parameter;


enum StripDirection {NEGATIVE, POSITIVE};
struct Strip {
    ll begin;
    ll end;
    StripDirection dir; // Direction of the strip: POSITIVE or NEGATIVE
    Strip() {}
    Strip(ll b, ll e, StripDirection d)
    {
        begin = b;
        end = e;
        dir = d;
    }
};

typedef struct Strip Strip;

class Permutation
{
    //TODO: this structures are 0 indexed...
    vector<ll> seq;
    vector<ll> pos;

    /* Data structures to store information about permutation strips:
     * For P = 3,1,2,-5,-4,
     *  -strips = [[0,0], [3,3], [1,2], [-5, -4]]
     *  -strips_pos = [[0,0],[1,1], [2,3], [4, 5]]
     */
    vector<pll > strips;
    vector<pll > strips_pos;

    // Current permuation parameters
    parameter param;
    // Permuation size
    int N;

    // flag regarding the Permution reduction state
    bool is_reduced;

public:
    typedef pair<Permutation, Reversal> perm_rev;
    Permutation() {};
    Permutation(parameter param, ll permI64);
    Permutation(parameter param);
    Permutation(parameter param, vector<ll> seq);
    Permutation(vector<ll> seq);

    static Permutation makeIdentity(parameter param, int n);

    static Permutation convertFromTo(Permutation &init, Permutation &end);
    //Reversal Model
    bool isReversalValid(int i, int j);
    bool reverse(Reversal rev);

    Permutation toUnsigned();
    //Methods for Almost Simetric Reversals
    int slice(int k);
    int sign(int k) const;
    int signAt(int k);
    int position(int k) const;
    bool isASRValidReversal(int i,int j);
    void updatePositions();
    bool isSorted(ll begin, ll end, bool asUnsigned = false);

    /* Methods for reduce permutation */
    void expand();
    void reduce();
    void setNormalMode();
    void setReducedExtendedMode();
    void setExtendedMode();
    Reversal convertReversal(Reversal rev);
    Reversal convertReversalToReduced(Reversal rev);

    //Auxiliar methods
    bool indexInBounds(int i) const;

    //General Methods
    bool isValid();
    void initSeq(vector<ll> seq);

    /* Return the absolute value at position k */
    ll at(int k) const;
    /* Return the value at position k, not the absolute value */
    ll get(int k) const;
    ll operator[](int k) const;

    int size() const;
    parameter getParameter() const;
    void setParameter(string param_name, double val);
    double getParameter(string param_name) const;

    static vector<Permutation> makeAllPermutations(ll n);
    static void makeAllPermutationsToFile(ll n, string fileName);
    static void makeRandPermutationsToFile(ll n,ll q, string fileName);

    static ll numPerm(ll n);
    void computePositiveSum(vector<ll> &sum) const;
    ll getNumOrientedPairs();
    ll messLevel() const;
    ll inversionLevel();
    ll numBreakPoints(bool asUnsigned = false);
    ll getNumDecreasingStrips();


    void  getBreakpointReversals(
				 vector<Permutation::perm_rev> &out_reversals, 
				 int cost_function
    ) const;
    /*void  getPreservingReversals(
    ) const;*/
    bool next();

    vector<Strip> getAllStrips() const;

    //Representation
    ll getI64();
    vector<ll> getBinary();
    pll getMedian(ll ini, ll fim);

    // vector<perm_rev> makeNeighborhood(int cost_function);
    //Operators
    int cmp(const Permutation &p) const;
    bool operator<(const Permutation &p) const;
    bool operator>(const Permutation &p) const;
    bool operator==(const Permutation &p) const;
    bool operator!=(const Permutation &p) const;

    //bool operator()

    //Input and Output
    friend ostream &operator<<(ostream &stream, Permutation p);
    friend istream &operator>>(istream &stream, Permutation &p);

    void printReversal(ostream &stream, Reversal r);

    void print(ostream &out, char sep) const;
    string toString() const;
    /* Mixed code */
    void changeSignsRandomly();

private:
    void set(int k, ll v);
    //Representation
    vector<ll> llToVector(ll k);

    bool invert(int i, int j);
};

struct PermutationComparer {
    bool operator() (
        const Permutation &a,
        const Permutation &b
    ) const
    {
        return a.cmp(b) == -1;
    }
};


#endif
