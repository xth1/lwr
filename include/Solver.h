#ifndef SOLVER_HEADER
#define SOLVER_HEADER 1
#include "HeaderDefault.h"
#include "Permutation.h"
class Result
{
public:
    vector<Reversal> reversals;
    ll totalCost;
    Result()
    {
        totalCost = 0;
    }
    void removeUnitaryReversals();
    friend ostream &operator<<(ostream &out, Result r);
    static void printResult(Result r, ostream &out);
    void print(ostream &out) const;
};

struct ResultComparer {
    /* Return true if a < b */
    bool operator() (
        const Result &a,
        const Result &b
    ) const
    {
        int sza = a.reversals.size();
        int szb = b.reversals.size();
        int min_sz = min(sza, szb);
        for (int i = 0; i < min_sz; i++) {
            Reversal ra = a.reversals[i];
            Reversal rb = b.reversals[i];
            if (ra < rb)
                return true;
            if (!(ra == rb))
                return false;
        }
        if (sza >= szb)
            return false;
        return true;
    }
};

struct ResultCostComparer {
    /* Return true if a < b */
    bool operator() (
        const Result &a,
        const Result &b
    ) const
    {
        return a.totalCost < b.totalCost;
    }
};


class Solver
{
protected:
    parameter param;
    parameter statistics;

public:
    Solver() {};
    Solver(parameter param);

    static bool validateResult(
        Permutation p,
        Result r,
        ostream &stream = cout,
        bool asUnsigned = false
    );
    //Ulisses modificou
    static Result calcTotalCost(Result &r, int cost_function);
    //Result calcTotalCost(Result &r);

    static void printResult(Result &r, ostream &out);
    //TODO: change to sort(const Permutation p)
    virtual Result sort(Permutation &p);
    virtual Result sort(Permutation &p, Result &result);
    virtual Result sort(Permutation &p, ll upper_bound);
    virtual Result distFromTo(Permutation &p1, Permutation &p2);
    virtual string getName();
    virtual void setName(string name);
    parameter getStatistics();
    virtual ~Solver() {}
protected:
    string name;
};

#endif
