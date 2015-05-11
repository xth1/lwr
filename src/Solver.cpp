#include "Solver.h"
#include<cmath>

/* Result methods */

void Result::print(ostream &out) const
{
    Result::printResult(*this, out);
}

void Result::removeUnitaryReversals()
{
    vector<Reversal> new_revs;
    // reset the result's total cost
    this->totalCost = 0;
    for (size_t i = 0; i < this->reversals.size(); i++) {
        Reversal rev = this->reversals[i];
        assert(rev.end >= rev.begin);
        // check if the current reversal is not unitary
        if (rev.end > rev.begin) {
            // add up the reversal
            new_revs.push_back(rev);
            // update the total cost
            this->totalCost += rev.cost;
        }
    }
    // set up the new_revs as the sequence of reversals
    this->reversals = new_revs;
}

void Result::printResult(Result r, ostream &out)
{
    out<<r.totalCost<<" "<<r.reversals.size();
    int sz = r.reversals.size();

    out<<" [";
    if (sz > 0) {
        for (int i = 0; i < sz - 1; i++)
            out<<"["<<r.reversals[i].begin<<","<<r.reversals[i].end<<"],";
        out<<"["<<r.reversals[sz-1].begin<<","<<r.reversals[sz-1].end<<"]";
    }
    out<<"]";
}

ostream &operator<<(ostream &out, Result r)
{
    r.print(out);
    return out;
}

Solver::Solver(parameter param)
{
    this->param = param;
    this->setName("Error: Don't have name");
}

Result Solver::sort(Permutation &p)
{
    cerr<<"Fatal error: Solver::sort(Permutation &p) has not been implemented"<<endl;
    cerr<<"\tWith parameters ("<<p<<endl;
    Result r;
    exit(1);
    return r;
}

Result Solver::sort(Permutation &p, Result &result)
{
    cerr<<"Fatal error: Solver::sort(Permutation &p, Result &result) has not been implemented"<<endl;
    cerr<<"\tWith parameters ("<<p<<" "<<result<<")"<<endl;
    Result r;
    exit(1);
    return r;
}

Result Solver::sort(Permutation &p, ll upper_bound)
{
    cerr<<"Fatal error: Solver::sort(Permutation &p, ll upper_bound = INF, ll lower_bound = 0) has not been implemented"<<endl;
    cerr<<"\tWith parameters ("<<p<<" "<<upper_bound<<")"<<endl;
    Result r;
    exit(1);
    return r;
}

Result Solver::distFromTo(Permutation &p1, Permutation &p2)
{
    cerr<<"Fatal error: Solver::distFromTo(Permutation &p1, Permutation &p2) has not been implemented"<<endl;
    cerr<<"\tWith parameters ("<<p1<<" "<<p2<<")"<<endl;
    Result r;
    exit(1);
    return r;
}
/* Check if a result is valid regarding the permutation p.
 * If the parameter asUnsigned is seted, only is necessary that  the elements of
 * the resulting permutation (after performing the reversals)
 * be correctly placed in order to be validated.
 */
bool Solver::validateResult(
    Permutation p,
    Result r,
    ostream &stream,
    bool asUnsigned
)
{
    vector<Reversal> vr = r.reversals;

    int i, sz;
    sz = vr.size();

    bool valid = true;

    if (p.isValid() == false) {
        stream<< "Permutation "<<p	<<" is invalid"<<endl;
        valid = false;
    }

    Permutation pp = p;
    Permutation id = pp.makeIdentity(p.getParameter(), p.size());

    ll ct = 0;
    for (i = 0; i < sz; i++) {
        if (vr[i].begin < 1 || vr[i].end > p.size()) {
            stream<<"At permutation "<<pp<<"["<<p<<"]"<<endl;
            stream<<"\tReversal ["<<vr[i].begin<<" , "<<vr[i].end<<"]"<<" is invalid."<<endl;
            valid = false;
        }
        pp.reverse(vr[i]);
        ct += vr[i].cost;
    }

    assert(pp.isValid());
    //cout<<"PERM AFTER REVERSALS "<<pp<<endl;
    if (!pp.isSorted(1, pp.size(), asUnsigned)) {
        stream<<"At permutation"<<pp<<"["<<p<<"]"<<endl;
        stream<<"\tThis is not a sorting sequence of reversals:"<<endl;
        for (i = 0; i <sz; i++) {
            stream<<"\t\t"<<vr[i].begin<<" "<<vr[i].end<<" "<<vr[i].cost<<endl;
        }
        valid = false;
    }

    if (ct != r.totalCost) {
        stream<<"At permutation "<<pp<<"["<<p<<"]"<<endl;
        stream<<"\tResult cost: "<<r.totalCost<<", real cost: "<<ct<<endl;
        valid = false;
    }

    return valid;
}

Result Solver::calcTotalCost(Result &r, int cost_function)
{
    //vector<Reversal> vr = r.reversals;
    ll ct = 0;

    Result nr;
    for (size_t i = 0; i < r.reversals.size(); i++) {
        Reversal rr;
        rr.begin = r.reversals[i].begin;
        rr.end   = r.reversals[i].end;
	rr.cost_function = cost_function;
	rr.calcCost();

        nr.reversals.push_back(rr);
        ct += rr.cost;
    }

    nr.totalCost = ct;

    return nr;
}

parameter Solver::getStatistics()
{
    return this->statistics;
}

string Solver::getName()
{
    return this->name;
}

void Solver::setName(string name)
{
    cerr<<"NAME : "<<name;
    this->name = name;
}


/*ll Solver::sort(){
	cout<<"Do Nothing..."<<endl;
	return NIL;
}*/
