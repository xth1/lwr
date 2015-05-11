#include "ImproveReversals.h"
#define DEBUG 1
#define TRACE(x...)
#define PRINT(x...) TRACE(printf(x));fflush(stdout)
#define WATCH(x) TRACE(cout << #x" = " << x << "\n")
ImproveReversals::ImproveReversals(vector<Solver *> &heuristics,
                                   string intervalMethod, parameter param,
                                   ostream *logOut)
{
    this->heuristics = heuristics;
    this->intervalMethod = intervalMethod;
    this->param = param;
    this->logOut = logOut;
}


bool ImproveReversals::hasOutputLog() {
    if (this->logOut != NULL && logOut->good())
        return true;
    return false;
}

/*
 * Remove  interval itv from R  and insert nR in the same place.
 * nR can be lower or greather than interval itv.
 */
void ImproveReversals::mergeResults(Result &R, Result &nR, pii itv)
{
    size_t i;
    size_t begin = itv.F;
    size_t end = itv.S;
    size_t nr_sz = nR.reversals.size();
    Result posR;


    /* Add  the  reversals that are before the interval */
    for (i = 0; i < begin; i++) {
        posR.reversals.push_back(R.reversals[i]);
    }

    /* Add new reversals */
    for (i = 0; i < nr_sz; i++) {
        posR.reversals.push_back(nR.reversals[i]);
    }

    /* Add the reversals that are after interval */
    for (i = end + 1; i < R.reversals.size(); i++) {
        posR.reversals.push_back(R.reversals[i]);
    }
    // Alterado por Ulisses
    R = Solver::calcTotalCost(posR, this->param["COST_FUNCTION"]);
}

void ImproveReversals::generateAnchorPoints(Result &R, int anchor_points,
        vector<int> &points)
{
    int sz = R.reversals.size();

    if (anchor_points == 3) {
        /* choose a value near the reversal sequence middle, in the
         * interval [sz/3 ... 2*(sz/3)]
         */
        int d = sz / 3;
        int m = d + rand() % (d + 1);

        //cout<<"sz "<<sz<<"r "<<r<<" d "<<d<<" m "<<m<<endl;
        points = vector<int>(3);
        points[0] = 0;
        points[1] = m;
        points[2] = sz - 1;
    } else if (anchor_points == 4) {
        /* choose a value near the reversal sequence middle, in the
         * interval [sz/3 ... 2*(sz/3)]
         */
        int d = sz / 3;
        int m = d + rand() % (d + 1);

        // left point
        int left = ceil(((double)(2 * m)) / 3.0);

        // right point
        int right = (2 * m + sz) / 3;
        //cout<<"sz "<<sz<<"r "<<r<<" d "<<d<<" m "<<m<<endl;
        //cout<<"left = "<<left<<", right = "<<right<<endl;
        points = vector<int>(4);
        points[0] = 0;
        points[1] = left;
        points[2] = right;
        points[3] = sz - 1;
    } else if (anchor_points == 5) {
        /* choose a value near the reversal sequence middle, in the
        * interval [sz/3 ... 2*(sz/3)]
        */
        int d = sz / 3;
        int m = d + rand() % (d + 1);

        // left point
        int left = ceil(((double)(m)) / 2.0);

        // right point
        int right = (m + sz) / 2;
        //cout<<"sz "<<sz<<"r "<<r<<" d "<<d<<" m "<<m<<endl;
        //cout<<"left = "<<left<<", right = "<<right<<endl;
        points = vector<int>(5);
        points[0] = 0;
        points[1] = left;
        points[2] = m;
        points[3] = right;
        points[4] = sz - 1;
    }
}

Result ImproveReversals::sortByAnchorPoints(Permutation &perm, Result &R,
        int anchor_points)
{
    vector<int> points;
    Result NR = R;
    size_t i, j;

    //cout<<"SORT BY ANCHOR"<<endl;
    /* Look for num of rounds */
    size_t num_round;
    if (this->param.find("NUM_ROUND") != this->param.end())
        num_round = (size_t) this->param["NUM_ROUND"];
    else
        num_round = NUM_ROUND;

    //cout<<"num_round: "<<num_round<<endl;
    for (i = 0; i < num_round; i++) {
        this->generateAnchorPoints(R, anchor_points, points);
        int shift = 0;
        //cout<<"ROUND "<<i<<endl;
        for (j = 0; j < points.size() - 1; j++) {
            pii interval(points[j] + shift, points[j + 1] + shift);
            // cout<<"\t"<<interval.F<<" "<<interval.S<<endl;
            int new_end = this->tryImprove(perm, NR, interval);
            shift = new_end - interval.S;
        }
    }
    return NR;
}

/* This method uses a GRASP strategy to choose a subsequence of reversals:
 * The probability of a subsequende be chosen is proportional to its
 * total cost (sum of costs).
 */
pii  ImproveReversals::findSubsequenceByGRASP(Result &R,
        ll min_length, ll max_length)
{
    pii itv;
    size_t i, j;
    size_t sz = R.reversals.size();
    vector<pii> subs;
    vector<ll> subs_cost;

    ll sum;
    ll sum_total = 0;
    for (i = 0; i < sz; i++) {
        sum = 0;
        for (j = i; j < sz; j++) {
            int length = j - i + 1;
            sum += R.reversals[j].cost;

            if (length < min_length || length > max_length)
                continue;

            subs.push_back(pii(i,j));
            subs_cost.push_back(sum * sum);
            sum_total += sum * sum;
        }
    }

    /* Randomly choose a value in interval [0, sum_total - 1] */
    assert(sum_total > 0);
    ll sum_ch = rand() % sum_total;

    /* Choose subsequence (interval) of reversals */
    ll sum_curr = 0;

    //cout<<"sum_total "<<sum_total<<" sum_ch "<<sum_ch<<endl;

    for (i = 0; i < subs.size(); i++) {
        sum_curr += subs_cost[i];
        itv = subs[i];
        //cout<<"\titv : "<<itv.F<<" "<<itv.S<<" "<<subs_cost[i]<<" curr "<<sum_curr<<endl;
        if (sum_curr >= sum_ch)
            break;
    }

    //cout<<"Choosed interval: "<<itv.F<<" "<<itv.S<<endl;
    return itv;
}

pii ImproveReversals::findDenseSubsequence(Result &R)
{
    int i;
    int sz = R.reversals.size();
    pii best_interval = pii(0,0);

    vector<ll> sum(sz + 1);

    /* Compute sub-sums */
    sum[0] = 0;
    for (i = 0; i < sz; i++) {
        sum[i + 1] = sum[i] + R.reversals[i].cost;
    }

    /* Choose begin */
    int begin = rand() % sz;

    /* Choose end */
    ll sub_sum = sum[sz] - sum[begin];
    ll ch_sum = rand() % sub_sum + 1;
    ll curr_sum = 0;
    int end = begin;
    for (i = begin; i < sz; i++) {
        curr_sum += R.reversals[i].cost;
        end = i;
        if (curr_sum >= ch_sum)
            break;
    }

    best_interval = pii(begin, end);


    /*for (i = 0; i < sz; i++) {
      double sum = 0;
      for (j = i; j < sz; j++) {
        sum += (double) R.reversals[i].cost;
        double dense = sum * ( j - i + 1);
        //cout<<"Sub ("<<i<<" , "<<j<<") : "<<dense<<endl;
        if (dense > max_dense) {
          max_dense = dense;
          best_interval = pii(i, j);
        }
      }
    }*/
    //cout<<"Best interval: "<<best_interval.F<<" "<<best_interval.S<<endl;
    return best_interval;
}

pii ImproveReversals::chooseInterval(
    Permutation &UNUSED(p),
    Result &R,
    int length
)
{
    pii interval;

    int sz = R.reversals.size();
    int begin =  rand() % (sz - length + 1);
    int end = min(begin + length, sz - 1);

    //cout<<"Sz "<<sz<<" length "<<length<<endl;
    interval = pii(begin, end);

    //cout<<"Best Interval: "<<interval.F<< " "<<interval.S<<endl;
    return interval;
}

int ImproveReversals::tryImprove(Permutation &p, Result &R, pii interval)
{
    Result nr = R;
    Result hnr;

    Result best_result;
    ll min_cost;

    size_t i, j;
    Permutation p1 = p, p2, new_p;

    /* Look for num of iterations */
    int num_it;
    if (this->param.find("NUM_IT") != this->param.end())
        num_it = this->param["NUM_IT"];
    else
        num_it = NUM_IT;

    /* Compute p1, the initial permutation of the interval */
    for (i = 0; i < (size_t)interval.F; i++) {
        p1.reverse(R.reversals[i]);
    }

    /* Compute p2, the final permutation of the interval */
    p2 = p1;
    ll interval_cost = 0;
    for (i = interval.F; i <= (size_t) interval.S; i++) {
        p2.reverse(R.reversals[i]);
        /* Compute cost of the interval */
        interval_cost += R.reversals[i].cost;
    }

    /* Gets new_p: p1 induced by p2 */
    new_p = Permutation::convertFromTo(p1, p2);
    min_cost = interval_cost;

    if (this->param.find("dump_windows") != this->param.end()) {
        fstream fout ("window_dump.txt", fstream::app | fstream::out);
        new_p.print(fout, ',');
        fout<<" ";
        for (i = interval.F; i < (size_t) interval.S; i++) {
            fout<<"["<<R.reversals[i].begin<<","<<R.reversals[i].end<<"],";
        }
        fout<<"["<<R.reversals[interval.S].begin<<","<<R.reversals[interval.S].end<<"]";
        fout<<endl;
    }

    /* Try find a better sequence of reversals to sort p1 */
    for (i = 0; i < (size_t) num_it; i++) {
        /* For each heuristic */
      // So existe uma heuristica disponivel na versao final. -> Ulisses
      //cout << "Num heuristics = " << heuristics.size() << endl;
        for (j = 0; j < heuristics.size(); j++) {
            Permutation new_p_copy = new_p;
            if (this->param["USE_UPPER_BOUND"]){

                hnr = heuristics[j]->sort(new_p_copy, min_cost);
	    } else {
	      //cout << "not using upper bound" << endl;
                hnr = heuristics[j]->sort(new_p_copy);
	    }

            if (p.getParameter("signed") == false) {
                hnr.removeUnitaryReversals();
            }
	    
	    // Alterado por Ulisses
            hnr = Solver::calcTotalCost(hnr, this->param["COST_FUNCTION"]);

            if (hnr.totalCost < min_cost) {
                min_cost = hnr.totalCost;
                best_result = hnr;
            }
        }
    }

    int end = interval.S;

    /* Check if it has an improvement */
    if (min_cost < interval_cost) {
        ll prev_cost = R.totalCost;
        this->mergeResults(R, best_result, interval);
        end = interval.F + best_result.reversals.size() - 1;

        // Assert the result's cost after merge with the improvement
        WATCH(prev_cost);
        WATCH(interval_cost);
        WATCH(best_result.totalCost);
        WATCH(R.totalCost);
        assert(prev_cost - interval_cost + best_result.totalCost == R.totalCost);
    }
    // Return end of sequence of reversals
    return end;
}

Result ImproveReversals::sort(Permutation &p, Result &R)
{
    parameter::iterator it;


    PRINT("IMPROVE REVERSALS SORT\n");
    if (this->intervalMethod == "GRASP")
        return this->sortByGraspInterval(p, R);
    else if (this->intervalMethod == "RAND")
        return this->sortByRandInterval(p, R);
    else if (this->intervalMethod == "PEAKS")
        return this->sortByEntropyPeaks(p, R);
    else if (this->intervalMethod == "ANCHOR") {
        // Default value
        int num_anchor_points = 3;
        if (this->param.find("NUM_ANCHOR_POINTS") != this->param.end())
            num_anchor_points = this->param["NUM_ANCHOR_POINTS"];
        return this->sortByAnchorPoints(p, R, num_anchor_points);
    } else if (this->intervalMethod == "WINDOW") {
        /* Default window size */

        //cout<<"Enter window: "<<this->param["WINDOW_SZ_BEGIN"]<<endl;
        int windows_size = 10;
        if (this->param.find("WINDOW_SIZE") != this->param.end() ||
                this->param["DYNAMIC_SIZE"]) {
            windows_size = this->param["WINDOW_SIZE"];
            return this->sortByWindow(p, R, windows_size);
        } else if (this->param.find("WINDOW_SZ_BEGIN") != this->param.end()) {
            //cout<<"Window increase size"<<endl;
            ll w_begin = this->param["WINDOW_SZ_BEGIN"];
            ll w_end = this->param["WINDOW_SZ_END"];
            ll w_inc = this->param["WINDOW_SZ_INC"];

            /* Check if w_inc sign is correct */
            if (Util::sign(w_end - w_begin) != Util::sign(w_inc))
                w_inc *= -1;

            /* Run for windows sizes in interval [w_begin, w_end] */
            Result nR = R;
            ll w = w_begin;

            // first trie
            nR = this->sortByWindow(p, nR, w);

            do {
                w += w_inc;
                /* Update nR results */
                nR = this->sortByWindow(p, nR, w);
            } while (w != w_end);

            /* Return result */
            return nR;
        }
    } else if (this->intervalMethod == "VALUE_WINDOW") {
        ll threshold = this->param["WINDOW_THRESHOLD"];
        return this->sortByValueWindow(p , R, threshold);
    }
    return R;
}

/*-------------------------sortByWindow-------------------------------*/

/* Select intertal that have max cost */
pii ImproveReversals::selectMaxCostInterval(vector<ppll> window)
{
    size_t i;
    ll max_cost = 0;
    pii max_interval;
    for (i = 0; i < window.size(); i++) {
        // cout<<"\t itv "<<I[i].F<<" "<<I[i].S<<" "<<C[i]<<endl;
        if (window[i].S > max_cost) {
            max_interval = window[i].F;
            max_cost = window[i].S;
        }
    }

    return max_interval;
}

bool cmp_ppll(ppll a, ppll b)
{
    return a.S > b.S;
}

/* Select a interval by a Rolette Wheel algorithm */
pii ImproveReversals::selectIntervalByRandom(vector<ppll> window,
        double &prob_chosen)
{
    int sz = window.size();
    int id = rand() % sz;
    prob_chosen = 1 / (double) window.size();

    return window[id].F;
}

/* Select a interval by a Rolette Wheel algorithm */
pii ImproveReversals::selectIntervalByRolette(
    vector<ppll> window,
    double &prob_chosen, 
    ll &itv_cost
)
{
    int i;
    double sum = 0;
    pii itv;
    int rolette_size = (int)window.size();
    bool prob_single_value = false;
    // subtract min value from all values
    bool use_minus_min = false;
    // factor of max probabilite interval that will be added in all itv
    double max_prob_factor_sum = 0;

    // Assertion concerning the input
    assert(window.size() > 0);

    PRINT("input\n");
    for (i = 0; i < rolette_size; i++) {
        PRINT("[%d %d]: %lld\n", window[i].F.F, window[i].F.S, window[i].S);
    }

    if (this->param.find("MAX_PROB_FACTOR_SUM") != this->param.end())
        max_prob_factor_sum = (double)this->param["MAX_PROB_FACTOR_SUM"];

    if (this->param.find("PROB_SIMPLE_VALUE") != this->param.end())
        prob_single_value = (bool)this->param["PROB_SIMPLE_VALUE"];

    if (this->param.find("USE_MINUS_MIN") != this->param.end())
        use_minus_min = (bool)this->param["USE_MINUS_MIN"];

    /* Look for rolette size parameter */
    if (this->param.find("ROLETTE_SIZE") != this->param.end()) {
        rolette_size = min(rolette_size, (int) this->param["ROLETTE_SIZE"]);
        //sort by cost to select intervals with greather cost
        std::sort(window.begin(), window.end(), cmp_ppll);
    } else if ((this->param.find("ROLETTE_PERCENT_DISCARD") != this->param.end())) {

        double percent = (double)this->param["ROLETTE_PERCENT_DISCARD"];
        percent /= 100.0;
        rolette_size -= percent * (double) rolette_size;
        // To avoid rollet size equals 0
        rolette_size = max(rolette_size, 1);
        //sort by cost to select intervals with greather cost
        std::sort(window.begin(), window.end(), cmp_ppll);
    }

    assert(rolette_size > 0);

    //cout<<"Rolette size: "<<rolette_size<<endl;
    vector<double> H(rolette_size);
    //cout<<"windows: "<<endl;

    double max_H = 0;
    double min_H = INF;
    for (i = 0; i < rolette_size; i++) {
        if (prob_single_value) {
            H[i] = window[i].S;
        } else {
            H[i] = window[i].S * window[i].S;
        }

        max_H = max(max_H, (double)H[i]);
        min_H = min(min_H, (double)H[i]);
        sum += H[i];
    }

    if (max_prob_factor_sum > 0) {
        sum = 0;
        double F = max_H * max_prob_factor_sum;
        for (i = 0; i < rolette_size; i++) {
            H[i] = H[i] + F;
            sum += H[i];
        }
    }

    if (use_minus_min) {
        sum = 0;
        for (i = 0; i < rolette_size; i++) {
            H[i] = H[i] - min_H + 1.0;
            sum += H[i];
        }
    }

    assert(H.size() > 0);
    assert(sum > 0);

    /* Run the rollet wheel algorithm */
    RandomGenerator random_generator;
    // Gets a value in the range [0, sum]
    double random_value = random_generator.getFloatPoint(0, sum);
    double curr_sum = 0;

    for (i = 0; i < rolette_size; i++) {
        curr_sum += H[i];
        if (curr_sum >= random_value) {
            itv = window[i].F;
            itv_cost = window[i].S;
            break;
        }
    }
    prob_chosen = (double) itv_cost / (double) sum;

    return itv;
}

/*
 * At each iteration a window (interval) is chosen to try improve
 */
Result ImproveReversals::sortByWindow(Permutation &perm, Result &R, int window_size)
{
    Result NR = R;

    size_t i, k;
    //cout<<"Sort by Window "<<endl;
    /* Gets number of rounds */
    int num_round;
    int sampling_rate;
    if (this->param.find("NUM_ROUND") != this->param.end())
        num_round = this->param["NUM_ROUND"];
    else
        num_round = NUM_ROUND;

    if (this->param.find("SAMPLING_RATE") != this->param.end())
        sampling_rate = this->param["SAMPLING_RATE"];
    else
        sampling_rate = SAMPLING_RATE;

    if (this->hasOutputLog()) {
        /* LOG: Output num_round and num_samples */
        (*logOut)<<num_round<<" "<<sampling_rate<<",";

        /* LOG: Output initial result */
        (*logOut)<<R.totalCost<<" "<<R.reversals.size();
    }

    for (i = 0; i < (size_t) num_round; i++) {
        //printf("ROUND %d\n", i);
      // cout << "Round " << i << endl; // Ulisses
        PRINT("ROUND %d\n", i);
        PRINT("Window size: %d\n ", window_size);
        //cout<<"ROUND "<<i<<endl;
        vector<ppll> window;
        // Dinamic window size
        if ( this->param["DYNAMIC_SIZE"]) {
            ll w_begin = this->param["WINDOW_SZ_BEGIN"];
            ll w_end = this->param["WINDOW_SZ_END"];
            // Choose randomly
            window_size = w_begin + rand() % (w_end - w_begin + 1);
            // Adjust the window size
            window_size = min(window_size, (int) NR.reversals.size());
        }

        assert(window_size > 0);
        // Probability of being picked of interval that was chosen
        double prob_chosen;

        // Get the windows
        ll sum_cost = 0;
        //Result::printResult(NR, cout);
        for (k = 0; k < (size_t) window_size; k++) {
            sum_cost += NR.reversals[k].cost;
            PRINT("%d: %lld\n", k, NR.reversals[k].cost);
        }
        // First window
        ppll itv_cost;

        // Others windows
        for (k = window_size - 1; k < NR.reversals.size(); k++) {
            if (k >= (size_t) window_size) {
                sum_cost += NR.reversals[k].cost;
                sum_cost -= NR.reversals[k - window_size].cost;
            }
            WATCH(sum_cost);
            itv_cost = ppll(pii(k - window_size + 1, k), sum_cost);
            window.push_back(itv_cost);
        }

        // Select a interval through a selection algorithm
        pii itv;
        if (this->param.find("WINDOW_MIN") != this->param.end()) {
            itv = this->selectMaxCostInterval(window);
        } else if (this->param.find("WINDOW_ROLETTE") != this->param.end()) {
            ll cost_itv;
	    // cout << "choosing by roulette" << endl; // Ulisses
            itv = this->selectIntervalByRolette(window, prob_chosen, cost_itv);
        } else if (this->param.find("WINDOW_RANDOM") != this->param.end()) {
            itv = this->selectIntervalByRandom(window, prob_chosen);
        } else  {
            cout<<"No window mode found"<<endl;
            exit(1);
        }


        // try improve the current result
        Permutation copy_perm = perm;
        this->tryImprove(copy_perm, NR, itv);

        // LOG: Outupt iteration results
        if (this->hasOutputLog()) {
            if ((i + 1) % (size_t) sampling_rate == 0 || (i+1) == (size_t) num_round) {
                (*logOut)<<","<<itv.F<<" "<<itv.S<<" "<<NR.totalCost<<" "<<NR.reversals.size()
                         <<" "<<prob_chosen;
            }
        }
    }
    // LOG: end of line
    if (this->hasOutputLog()) {
        (*logOut)<<endl;
    }
    return NR;
}

Result ImproveReversals::sortByValueWindow(Permutation &perm, Result &R, ll threshold)
{

    Result NR = R;

    int i;

    /* Gets number of rounds */
    int num_round;
    if (this->param.find("NUM_ROUND") != this->param.end())
        num_round = this->param["NUM_ROUND"];
    else
        num_round = NUM_ROUND;

    /* For each round */
    for (i = 0; i < num_round; i++) {
        /* Get current size of sequence of reversals */
        int N = NR.reversals.size();
        vector<ppll> window;

        double prob_chosen;

        //DEBUG
        /* cout<<"Cost "<<endl;
         for (k = 0; k < NR.reversals.size(); k++) {
           cout<<NR.reversals[k].cost<<" ";
         }
         cout<<endl;
        */
        /* Get the windows */
        ll curr_cost = 0;
        int begin;
        int end;
        ppll itv_cost;

        /* Find windows (i, j) such that cost[i..j] >= threshold and
         * cost[i .. j - 1] < threshold. 0 <= i, j < N */
        for (begin = 0, end = 0; end < N ; end++) {

            // New cost
            ll nc = NR.reversals[end].cost;

            // Increase the current cost
            curr_cost += nc;

            //cout<<" cost "<<curr_cost<<" begin "<<begin<<" end "<<end<<endl;
            if (curr_cost >= threshold) {
                /* Look for windows */
                while (curr_cost >= threshold) {
                    // Add interval
                    itv_cost = ppll(pii(begin, end), curr_cost);
                    window.push_back(itv_cost);
                    //cout<<"Add "<<itv_cost.F.F<<" "<<itv_cost.F.S<<" "<<itv_cost.S<<endl;
                    curr_cost -= NR.reversals[begin].cost;
                    begin++;
                }
            }
        }

        /* Selected window (interval) */
        pii itv;

        /* Mode for select interval (window) */
        if (this->param.find("WINDOW_MIN") != this->param.end()) {
            itv = this->selectMaxCostInterval(window);
        } else if (this->param.find("WINDOW_ROLETTE") != this->param.end()) {
            ll itv_cost;
            itv = this->selectIntervalByRolette(window, prob_chosen, itv_cost);
        } else {
            cout<<"No window mode found"<<endl;
            exit(1);
        }

        //cout<<"Window: "<<itv.F<<" "<<itv.S<<" "<<endl;
        Permutation p = perm;
        this->tryImprove(p, NR, itv);
    }

    return NR;
}

//Choose interval ramdonly
Result ImproveReversals::sortByRandInterval(Permutation &p, Result &R)
{
    int i;
    int factor = 1;
    vector<pii> peaks;
    Result nr = R;
    pii interval;
    int length = nr.reversals.size();
    findEntropyPeaks(p, R, peaks);
    for (i = 0; i < NUM_ROUND; i++) {
        //update
        length = min(length, (int)nr.reversals.size());
        int min_length = length / 3;

        if (length >= min_length)
            length -= factor;

        interval = this->chooseInterval(p, nr, length);
        tryImprove(p, nr, interval);
    }

    findEntropyPeaks(p, nr, peaks);
    return nr;
}


Result ImproveReversals::sortByGraspInterval(Permutation &p, Result &R)
{
    int i;
    Result nr = R;
    pii interval;
    int length = nr.reversals.size();
    int min_length = length / 3;
    int max_length = length;

    for (i = 0; i < NUM_ROUND; i++) {
        //update fields
        length = min(length, (int)nr.reversals.size());
        min_length = max(min_length, length / 3);
        max_length = min(max_length, length);

        interval = this->findSubsequenceByGRASP(nr, min_length, max_length);
        tryImprove(p, nr, interval);
    }
    return nr;
}


/* Entropy Peaks */
Result ImproveReversals::sortByEntropyPeaks(Permutation &p, Result &R)
{
    size_t i, j;
    Result NR = R;

    /* Gets number of rounds */
    size_t num_round;

    if (this->param.find("NUM_ROUND") != this->param.end())
        num_round = this->param["NUM_ROUND"];
    else
        num_round = NUM_ROUND;

    //cout<<"num_round: "<<num_round<<endl;
    for (i = 0; i < num_round; i++) {
        vector<pii> peaks;
        this->findEntropyPeaks(p, NR, peaks);

        // find peak with maximum length
        pii max_interval, interval;
        ll max_length = 0;
        for (j = 0; j < peaks.size(); j++) {
            ll length = peaks[j].S - peaks[j].F + 1;
            //printf("%lld: %lld , (%d, %d)\n",j , length, peaks[j].F, peaks[j].S);
            if (length > max_length) {
                max_interval = peaks[j];
                max_length = length;
            }
        }

        //cout<<"max length "<<max_length<<endl;
        //cout<<"max Interval "<<max_interval.F<<" "<<max_interval.S<<endl;
        //Given (i,j) -> (i,n)
        interval = pii(max_interval.F, max_interval.S);
        //cout<<"choosed Interval "<<interval.F<<" "<<interval.S<<endl;
        // was not found peaks..
        if (max_length == 0)
            break;
        tryImprove(p, NR, interval);
    }

    return NR;
}

void ImproveReversals::findEntropyPeaks(Permutation &perm, Result &R, vector<pii> &peaks)
{
    int sz = R.reversals.size();

    /* Entropy of sequence of reversals */
    vector<ll> E(sz + 1);
    Permutation curr_perm;

    int i;

    /* Compute entropy */
    curr_perm = perm;
    E[0] = perm.messLevel();

    /* Gets min peak length */
    int min_peak_length = 3;

    if (this->param.find("MIN_PEAK_LENGTH") != this->param.end())
        min_peak_length = this->param["MIN_PEAK_LENGTH"];
    else
        min_peak_length = MIN_PEAK_LENGTH;

    //printf("Entropy\n");
    //printf("\tE[%d] = %lld :",0, E[0]);
    //cout<<curr_perm<<endl;
    for (i = 0; i < sz; i++) {
        curr_perm.reverse(R.reversals[i]);
        E[i + 1] = curr_perm.messLevel();
        //printf("\tE[%d] = %lld :",i + 1, E[i + 1]);
        //cout<<curr_perm<<endl;
    }

    /* Find peaks */
    int begin = 0;
    for (i = 1; i <= sz; i++) {
        if (E[i] < E[begin]) {
            /* Found a peak that have size at least min_peak_length */
            if (i - begin >= min_peak_length) {
                //printf("peak: (%d, %d)\n", begin + 1, i - 1);
                /* add peak, convert to 0-indexed scale*/
                peaks.push_back(pii(begin , i - 2));
            }
            /* Update begin */
            begin = i;
        }
    }
}
