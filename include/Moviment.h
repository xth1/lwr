#ifndef MOVIMENT_HEADER
#define MOVIMENT_HEADER 1
#include "HeaderDefault.h"


class Moviment
{
    //virtual class
    virtual double getCost();
};

class Reversal: public Moviment
{
public:
    ll begin;
    ll end;
    ll cost;
    // Alterado por Ulisses toda vez que aparecer cost_function
    int cost_function;
    void calcCost();
    Reversal() {};
    Reversal(ll begin,ll end,ll cost, int cost_function);
    // Reversal(ll begin,ll end, int cost_function);
    virtual double getCost();

    friend ostream &operator<<(ostream &stream, Reversal r);
    friend istream &operator>>(istream &stream, Reversal &r);
    bool operator==(Reversal &rev) const;
    bool operator<(Reversal &rev) const;
    bool operator!=(Reversal &rev) const;
};


#endif
