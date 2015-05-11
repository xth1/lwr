#include "HeaderDefault.h"
#include "HeuristicChooseReversal.h"

#define TRACE(x...)
#define PRINT(x...) TRACE(printf(x));fflush(stdout)
#define WATCH(x) TRACE(cout << #x" = " << x << "\n")


/* Initialize static variable */
bool HeuristicChooseReversal::is_heuristics_ready = false;
vector<HeuristicChooseReversal *> HeuristicChooseReversal::heuristics = vector<HeuristicChooseReversal *>();

/*======================HeuristicChooseReversal=======================*/


void HeuristicChooseReversal::setParameter(string name, double value)
{
    this->param[name] = value;
}

void HeuristicChooseReversal::setParameter(parameter &param)
{
    this->param = param;

    parameter::iterator it;
    /*
    cout<<" Choose heuristic parameters" <<endl;
    for (it = param.begin(); it != param.end(); it++) {
      cout<<it->F<<" "<<it->S<<endl;
    }*/
}

double HeuristicChooseReversal::getParameter(string name)
{
    return this->param[name];
}

bool HeuristicChooseReversal::hasParameter(string name)
{
    if (this->param.find(name) != this->param.end())
        return true;
    return false;
}

/* Gets name of Heuristic */
string HeuristicChooseReversal::getName()
{
    return "UnImplemented";
}

/* List the current implmented heuristics */
void HeuristicChooseReversal::BuildHeuristics()
{
    if (HeuristicChooseReversal::is_heuristics_ready == true)
        return;
    HeuristicChooseReversal::is_heuristics_ready = true;
    /* Add here objects of new Heuristics */
    HeuristicChooseReversal::heuristics.push_back(new HCRMinLength());
    HeuristicChooseReversal::heuristics.push_back(new HCRMaxLength());
    HeuristicChooseReversal::heuristics.push_back(new HCRBenefitLoss());
    HeuristicChooseReversal::heuristics.push_back(new HCRLexMin());
    HeuristicChooseReversal::heuristics.push_back(new HCRTest());
    HeuristicChooseReversal::heuristics.push_back(new HCRGRASPByEntropy());
    HeuristicChooseReversal::heuristics.push_back(new HCRRandom());
}

void HeuristicChooseReversal::destroyHeuristics()
{
    if (HeuristicChooseReversal::is_heuristics_ready == false)
        return;
    for (size_t i = 0; i < HeuristicChooseReversal::heuristics.size(); i++) {
        delete HeuristicChooseReversal::heuristics[i];
    }
    HeuristicChooseReversal::heuristics.clear();
}

HeuristicChooseReversal * HeuristicChooseReversal::findHeuristic(string name)
{
    if (HeuristicChooseReversal::is_heuristics_ready == false)
        HeuristicChooseReversal::BuildHeuristics();

    for (size_t i = 0; i < HeuristicChooseReversal::heuristics.size(); i++) {
        if (HeuristicChooseReversal::heuristics[i]->getName() == name)
            return HeuristicChooseReversal::heuristics[i];
    }
    cerr<<"Heuristic "<<name<<" not found"<<endl;
    return NULL;
}

/*==========================Roulette========================*/
/* Choose a reversal between the first num_rev by Roulette Whell
 * Selection algorithm.
 */
int Roulette::chooseByRoulette(vector<pdi> &score, int num_rev)
{
    int i;
    double sum = 0;

    //Check if the set of reversals is empty or not
    if (score.size() == 0)
        return -1;
    // Calc total sum
    for (i = 0; i < num_rev; i++) {
        sum += score[i].F;
    }

    // choose romdonly a val in interval [0..sum];
    RandomGenerator random_generator;
    double random_value = random_generator.getFloatPoint(0, sum);

    // current sum
    double curr_sum = 0;
    int chosen_id = -1;

    // find id
    for (i = 0; i < num_rev; i++) {
        curr_sum += score[i].F;
        chosen_id = score[i].S;
        if (curr_sum >= random_value)
            break;
    }

    return chosen_id;
}

/*============================GRASP============================*/
int HCRGRASP::choose(
    Permutation perm,
    vector<Reversal> revs,
    vector<ll> nops,
    parameter local_param
)
{
    int sz = nops.size();
    int i;
    /* Number of reversals that will be selected for run roulette algorithm */
    int num_selected_rev = revs.size();
    ll max_nop = 0;

    if (HeuristicChooseReversal::getParameter("NUM_SELECTED_REV")) {
        num_selected_rev = HeuristicChooseReversal::getParameter("NUM_SELECTED_REV");
    }

    /* Check if just optimal reversals will be used */
    if (HeuristicChooseReversal::getParameter("JUST_OPTIMAL_REV")) {
        max_nop = *max_element(nops.begin(), nops.end());
    }

    // benefit of reversals
    vector<double> B(sz);

    // score of reversals
    vector< pdi > score;

    // min negative benefit
    double min_neg_ben = 0;

    for (i = 0; i < sz; i++) {
        B[i] = this->computeBenefit(perm, revs[i], local_param);
        WATCH(revs[i]);
        WATCH(B[i]);
        if (B[i] < min_neg_ben)
            min_neg_ben = B[i];
    }

    // calc scores
    min_neg_ben = abs(min_neg_ben);

    for (i = 0; i < sz; i++) {
        pdi sc;

        // Check it's a optimal reversal
        if (HeuristicChooseReversal::getParameter("JUST_OPTIMAL_REV")
                && nops[i] != max_nop
           ) {
            continue;
        }

        if (HeuristicChooseReversal::getParameter("AVOID_NON_POSITIVE_LOSS")
                && B[i] <= 0
           ) {
            continue;
        }

        /* All scores will be positive ( >= 1)*/
        sc.F = min_neg_ben + B[i] + 1;

        if (HeuristicChooseReversal::getParameter("SQUARE_OF_SCORE")) {
            sc.F *= sc.F;
        }

        if (HeuristicChooseReversal::hasParameter("COEFFICIENT")) {
            double c = HeuristicChooseReversal::getParameter("COEFFICIENT");
            //PRINT("COEFFICIENT %lf\n", c);
            sc.F *= c;
        }

        if (HeuristicChooseReversal::getParameter("EXPONENTIAL_BASE_2")) {
            //PRINT("EXPONENTIAL_BASE_2");
            sc.F = pow(2.0, sc.F);
        }

        WATCH(sc.F);

        sc.S = i;
        score.push_back(sc);
    }

    //cout<<"NSR: "<<num_selected_rev<<" NR "<<score.size()<<endl;
    sort(score.begin(), score.end(), pdiComparer());

    // num of reversals that will be used
    int num_rev = min(num_selected_rev, (int)score.size());

    // Choose a reversal by Roulette Whell Selection algorithm
    return Roulette::chooseByRoulette(score, num_rev);
}

string HCRGRASP::getName()
{
    return "GRASP";
}

/*=============================MinLength==============================*/

/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum length.
 */
int HCRMinLength::choose(
    Permutation perm,
    vector<Reversal> revs,
    vector<ll> nops,
    parameter local_param
)
{
    size_t i;
    int min_length = INF;
    int id_best = -1;
    ll max_nop = *max_element(nops.begin(), nops.end());

    if (HeuristicChooseReversal::getParameter("GRASP")) {
        ll max_length = 0;
        for (size_t i = 0; i < revs.size(); i++) {
            max_length = max(max_length, revs[i].cost);
        }
        local_param["max_rev_cost"] = max_length;
        return HCRGRASP::choose(perm, revs, nops, local_param);
    }

    for (i = 0; i < revs.size(); i++) {
        if (HeuristicChooseReversal::getParameter("JUST_OPTIMAL_REV")
                && nops[i] != max_nop
           ) {
            continue;
        }

        if (revs[i].cost < min_length) {
            min_length = nops[i];
            id_best = i;
        }
    }

    return id_best;
}

double HCRMinLength::computeBenefit(
    const Permutation &UNUSED(perm),
    Reversal rev,
    parameter param
)
{
    ll max_rev_cost = param["max_rev_cost"];

    double b = max_rev_cost - rev.cost + 1;

    return b;
}

string HCRMinLength::getName()
{
    return "MinLength";
}
/*=============================Random==============================*/

/*
 * Choose a reversal randomly
 */
int HCRRandom::choose(
    Permutation perm,
    vector<Reversal> revs,
    vector<ll> nops,
    parameter UNUSED(local_param)
)
{
    return HCRGRASP::choose(perm, revs, nops);
}

double HCRRandom::computeBenefit(
    const Permutation &UNUSED(perm),
    Reversal UNUSED(rev),
    parameter UNUSED(param)
)
{
    //constant probability
    return 1;
}

string HCRRandom::getName()
{
    return "Random";
}

/*=============================MaxLength==============================*/

/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum length.
 */
int HCRMaxLength::choose(
    Permutation UNUSED(perm),
    vector<Reversal> revs,
    vector<ll> nops,
    parameter UNUSED(local_param)
)
{
    size_t i;
    int id_best = 0;
    ll max_nop = nops[0];

    for (i = 1; i < revs.size(); i++) {
        if (nops[i] > max_nop ||
                (nops[i] == max_nop && revs[i].cost > revs[id_best].cost)) {
            max_nop = nops[i];
            id_best = i;
        }
    }
    return id_best;
}

string HCRMaxLength::getName()
{
    return "MaxLength";
}

/*===========================EntropyLoss==============================*/

double HCRBenefitLoss::calcLoss(const Permutation &perm, const Reversal rev)
{
    double e1;
    double e2;
    double loss;
    // Makes a copy
    Permutation new_perm = perm;
    PRINT("CALC LOSS\n");
    // Original Permutation
    e1 = 0;
    if (HeuristicChooseReversal::hasParameter("ENTROPY")) {
        PRINT("ENTROPY\n");
        e1 += (double) new_perm.messLevel();
    }

    if (HeuristicChooseReversal::hasParameter("INVERSION")) {
        PRINT("INVERSION\n");
        e1 += (double) new_perm.inversionLevel();
    }

    new_perm.reverse(rev);
    // After performs reversal
    e2 = 0;
    if (HeuristicChooseReversal::hasParameter("ENTROPY")) {
        e2 += (double) new_perm.messLevel();
    }

    if (HeuristicChooseReversal::hasParameter("INVERSION")) {
        e2 += (double) new_perm.inversionLevel();
    }

    // Compute entropy loss
    loss = (e1 - e2) / (double) rev.cost;
    return loss;
}

double HCRBenefitLoss::computeBenefit(
    const Permutation &perm,
    Reversal rev,
    parameter UNUSED(local_param)
)
{
    return this->calcLoss(perm, rev);
}

/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum length.
 */
int HCRBenefitLoss::choose(
    Permutation perm,
    vector<Reversal> revs,
    vector<ll> nops,
    parameter UNUSED(local_param)
)
{

    if (HeuristicChooseReversal::hasParameter("GRASP")) {
        return HCRGRASP::choose(perm, revs, nops);
    }

    ll max_nop = *max_element(nops.begin(), nops.end());
    double max_loss = -INF;
    int id_best = -1;
    assert(revs.size() > 0);
    assert(revs.size() == nops.size());

    PRINT("\nSTART CHOOSE------------------------------------\n");
    WATCH(max_nop);
    WATCH(revs.size());

    for (size_t i = 0; i < revs.size(); i++) {
        double loss = this->calcLoss(perm, revs[i]);
        WATCH(nops[i]);
        WATCH(revs[i]);
        WATCH(loss);
        if (HeuristicChooseReversal::hasParameter("JUST_OPTIMAL_REV")
                && nops[i] != max_nop
           ) {
            continue;
        }

        if (HeuristicChooseReversal::hasParameter("AVOID_NON_POSITIVE_LOSS")
                && loss <= 0
           ) {
            continue;
        }

        if (loss > max_loss) {
            max_nop = nops[i];
            max_loss = loss;
            id_best = i;
        }
        WATCH(id_best);
    }

    assert(id_best == - 1 || (id_best >= 0 && (size_t)id_best < revs.size()));

    return id_best;
}

string HCRBenefitLoss::getName()
{
    return "BenefitLoss";
}


/*============================LexMin============================*/

bool HCRLexMin::isLower(const Reversal &r1, const Reversal r2)
{
    if (r1.begin < r2.begin)
        return true;
    else if (r1.begin == r2.begin)
        return r1.end < r2.end;
    return false;
}

/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum length.
 */
int HCRLexMin::choose(
    Permutation UNUSED(perm),
    vector<Reversal> revs,
    vector<ll> UNUSED(nops),
    parameter UNUSED(local_param)
)
{
    size_t i;
    int id_best = 0;

    for (i = 1; i < revs.size(); i++) {
        if (this->isLower(revs[i], revs[id_best]))
            id_best = i;
    }

    return id_best;
}

string HCRLexMin::getName()
{
    return "LexMin";
}

/*============================Test============================*/
/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum length.
 */

double HCRTest::calcHeuristic(Permutation &perm)
{
    double ml = (double) perm.messLevel();
    double hc =  ml;
    return hc;
}

double HCRTest::calcBenefit(const Permutation &perm, const Reversal rev)
{
    double e1;
    double e2;
    double benefit;
    // Makes a copy
    Permutation new_perm = perm;

    ll n = new_perm.size();
    double idp = (double)min(rev.begin, n - rev.end);

    // Original Permutation
    e1 = (double) this->calcHeuristic(new_perm);
    new_perm.reverse(rev);
    // After performs reversal
    e2 = (double) this->calcHeuristic(new_perm);

    // Compute benefit
    benefit = (e1 - e2) / (((double) rev.cost) * idp);
    return benefit;
}

int HCRTest::choose(
    Permutation perm,
    vector<Reversal> revs,
    vector<ll> UNUSED(nops),
    parameter UNUSED(local_param)
)
{
    size_t i;
    int id_best = -1;
    double max_benefit = -INF;
    double benefit;

    for (i = 0; i < revs.size(); i++) {
        benefit = this->calcBenefit(perm, revs[i]);
        if (benefit > max_benefit) {
            max_benefit = benefit;
            id_best = i;
        }
    }

    if (max_benefit <= 0)
        return -1;
    return id_best;
}

string HCRTest::getName()
{
    return "Test";
}

/*============================GRASPByEntropy============================*/
/* Grasp By Positive entropy */
double HCRGRASPByEntropy::calcHeuristic(
    Permutation &perm,
    vector<Reversal> &revs,
    vector<ll> &UNUSED(nops),
    vector<double> &H
)
{
    size_t sz = revs.size();
    H = vector<double>(sz);
    Permutation p2;
    /* Sum of heuristics values */
    double H_sum = 0;
    size_t i;
    ll n = perm.size();
    double max_loss = 0;

    for (i = 0; i < revs.size(); i++) {
        p2 = perm;
        //p2.reverse(revs[i]);
        H[i] = HCRBenefitLoss::calcLoss(p2, revs[i]) * (double) n;
        max_loss = max(max_loss, H[i]);
        // cout<<revs[i]<<" : "<<H[i]<<endl;

        // just positive loss
        if (H[i] > 0)
            H_sum += H[i];
    }

    /* returns heuristics sum */
    return H_sum;
}

/*
 * Choose reversal with maximun NOP. In case of tie,
 * choose any with minimum length.
 */

int HCRGRASPByEntropy::choose(
    Permutation perm,
    vector<Reversal> revs,
    vector<ll> nops,
    parameter UNUSED(local_param)
)
{
    size_t i;
    int best_id = -1;
    double H_sum;
    vector<double> H;
    H_sum = this->calcHeuristic(perm, revs, nops, H);

    /* Don't exist reversal with positive loss */
    if (H_sum <= 0) {
        return -1;
    }

    size_t sz = revs.size();
    /* Randomly choose a value in interval [1, H_sum] */
    double random_val = rand() % (ll) round(H_sum);
// cout<<"H_sum "<<H_sum<<endl;
// cout<<"randon val "<<random_val<<endl;
    i = 0;
    /* Current sum of heuristics values */
    double curr_sum = 0;

    /* Find for a reversal i that:
     * H[0] + H[1] ... + H[i] >= random_val
     */
    for (i = 0; curr_sum < random_val && i < sz; i++) {

        /* Ignore reverals with negative loss */
        if (H[i] <= 0)
            continue;
        curr_sum += H[i];
        best_id = i;
        // cout<<"Add rev "<<i<<" "<<revs[i]<<" : "<<H[i]<<" Sum: "<<curr_sum<<endl;
    }

    //cout<<"Chosed "<<best_id<<endl;

    return best_id;
}

string HCRGRASPByEntropy::getName()
{
    return "GRASPByEntropy";
}
