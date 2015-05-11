#include "ExactUnitary.h"
#include "pstream.h"
#include <sstream>
#include <stdlib.h>
#define TRACE(x...)
#define PRINT(x...) TRACE(printf(x));fflush(stdout)
#define WATCH(x) TRACE(cout << #x" = " << x << "\n")

ExactUnitary::ExactUnitary(
    parameter param,
    vector<HeuristicChooseReversal *> hcr_seq
)
{
    this->param = param;
    this->hcr_seq = hcr_seq;
}

void ExactUnitary::makeNearestGreatherArray(Permutation &perm, vector<ll> &M)
{
    int i;
    int n;
    n = perm.size();
    /* Alloc M */
    M = vector<ll>(n + 1);

    /* Stack that stores nearest greather element */
    stack<ll> S;

    /* Initial state */
    M[0] = n;
    S.push(n);

    for (i = 1; i <= n; i++) {
        /* Compute M[i] */
        if (perm.at(i - 1) > perm.at(i))
            S.push(perm.at(i - 1));
        else {
            /* Pop from S all entries that smaller than perm[i].
             * S never will be empty, because M[0] = n */

            while (!S.empty() && S.top() < perm.at(i)) {
                S.pop();
            }
        }

        /* S.top() is the nearest greater element of perm[i]  */
        M[i] = S.top();
    }
}

void ExactUnitary::buildBranchTree(vector<Component> &components, Tree &T)
{
    this->visited = vector<bool>(components.size() + 2);

    /* Call dfs to root node */
    dfs(components, T, 0);

}

bool ExactUnitary::dfs(vector<Component> &components, Tree &T, int v)
{
    size_t i;
    int idc = T[v].id_component;
    bool can_delete = true;
    bool is_unoriented = false;

    if (idc != SQUARE && components[idc].oriented == false)
        is_unoriented = true;

    for (i = 0; i < T[v].adj.size(); i++) {
        int w = T[v].adj[i];
        can_delete &= dfs(components, T, w);
    }

    /* Delete current node if it is not a  a leaf node or
     * all nodes in its childhood are square or oriented node*/
    if (can_delete && is_unoriented == false ) {
        printf("Delete Node %d, Component [%d, %d], oriented = %d\n",
               v, components[idc].begin, components[idc].end, components[idc].oriented);
        return true;
    } else {
        return false;
    }
}

void ExactUnitary::makeNearestSmallerArray(Permutation &perm, vector<ll> &m)
{
    int i;
    int n;
    n = perm.size();
    /* Alloc m */
    m = vector<ll>(n + 1);

    /* Stack that stores nearest smaller element */
    stack<ll> S;

    /* Initial state */
    m[0] = 0;
    S.push(n);

    for (i = 1; i <= n; i++) {
        /* Compute M[i] */
        if (perm.at(i - 1) < perm.at(i))
            S.push(perm.at(i - 1));
        else {
            /* Pop from S all entries that greather than perm[i].
             * S never will be empty, because M[0] = n */
            while (!S.empty() && S.top() > perm.at(i)) {
                S.pop();
            }
        }

        /* S.top() is the nearest smaller element of perm[i]  */
        m[i] = S.top();
    }
}

void ExactUnitary::findComponents(Permutation &perm, vector<Component> &components)
{
    int i;
    int s;
    int n = perm.size();
    vector<ll> M, m;
    stack<ll> S1, S2;

    /* Oriented mark */
    vector<bool> mark1(n + 2, false);
    vector<bool> mark2(n + 2, false);

    /* Initialize stacks */
    S1.push(0);
    S2.push(0);

    /* Compute M */
    this->makeNearestGreatherArray(perm, M);

    /* Compute m */
    this->makeNearestSmallerArray(perm, m);

    for (i = 1; i <= n; i++) {

        /* Find direct components*/
        s = S1.top();
        /* While a interval begin s is not found */
        while (perm.at(s) > perm.at(i) || M[s] < perm[i]) {
            S1.pop();
            mark1[S1.top()] = mark1[S1.top()] || mark1[s];
            s = S1.top();
        }
        /* Check if the interval [s..i] is a valid component */
        if (perm.signAt(i) == 1 && M[i] == M[s]
                && (i - s == perm.at(i) - perm.at(s))) {
            /* TODO: Check correctness of this: */
            if (i - s == 1) {
                mark1[s] = true;
            }
            components.push_back(Component(s, i, true, mark1[s]));

            mark1[S1.top()] = false;
        }

        /* Find reversed components */
        s = S2.top();
        /* While a interval begin s is not found */
        while ((perm.at(s) < perm.at(i) || m[s] > perm.at(i))
                && (s > 0)) {
            S2.pop();
            mark2[S2.top()] = mark2[S2.top()] || mark2[s];
            s = S2.top();
        }

        /* Check if the interval [s..i] is a valid component */
        if (perm.signAt(i) == -1 && m[i] == m[s]
                && (i - s == perm.at(s) - perm.at(i))) {

            /* TODO: Check correctness of this: */
            if (i - s == 1) {
                mark2[s] = true;
            }

            components.push_back(Component(s, i, false, mark2[s]));
            mark2[S2.top()] = false;
        }

        /* Update stacks */
        if (perm.signAt(i) == 1) {
            S1.push(i);
        } else {
            S2.push(i);
        }

        /* Update marks */
        mark1[S1.top()] = perm.signAt(i) == -1;
        mark2[S2.top()] = perm.signAt(i) == 1;

    }
}

void print_tree(vector<Component> &C, Tree &T)
{
    size_t i, j;

    for (i = 0; i < T.size(); i++) {
        printf("Node %lu [%d]", i, T[i].id_component);

        if (T[i].id_component != -1) {
            printf(", comp: [%d %d %d]\n",
                   C[T[i].id_component].begin,
                   C[T[i].id_component].end,
                   C[T[i].id_component].direct);
        } else
            printf("\n");
        printf("\tAdj\n");
        for (j = 0; j < T[i].adj.size(); j++) {
            printf("\t\t%d\n", T[i].adj[j]);
        }
    }

}

void  ExactUnitary::buildComponentsTree(Permutation &perm, vector<Component> &C
                                        , Tree &T)
{

    size_t i;
    size_t n = perm.size();

    /* Current id */
    int curr_id = 0;
    /* CBegin[i] stores a id of the component that have begin at
     * CBegin[i] */
    vector<int> CBegin(n + 2, -1);
    /* CBegin[i] stores a id of the component that have end at
     * CBegin[i] */
    vector<int> CEnd(n + 2, -1);

    /* Compute CBegin and CEnd */
    for (i = 0; i < C.size(); i++) {
        CBegin[C[i].begin] = i;
        CEnd[C[i].end] = i;
    }

    /* Most recent square node */
    Node q(curr_id++, SQUARE, SQUARE);
    T.push_back(q);

    /* Most recent round node */
    Node p(curr_id++, q.id, CBegin[0]);
    T.push_back(p);
    T[q.id].adj.push_back(p.id);

    for (i = 1; i < n; i++) {
        /* if there is a component starting at position i */
        if (CBegin[i] != -1) {
            /* If there is no component ending at position i */
            if (CEnd[i] == -1) {
                /* Create a new square node q as a child of p */
                q = Node(curr_id++, p.id, SQUARE);
                T.push_back(q);
                T[p.id].adj.push_back(q.id);
            }
            /* Create a new round node p (representing the new component) as
             * child of q */
            p = Node(curr_id++, q.id, CBegin[i]);
            T.push_back(p);
            T[q.id].adj.push_back(p.id);
        } else if (CEnd[i] != -1) {
            p = T[q.id_parent];
            q = T[p.id_parent];
        }
    }
    print_tree(C, T);

}

Reversal ExactUnitary::mergeHurdles(pii h1, pii h2)
{
    if (h2.F < h1.F)
        swap(h1, h2);

    Reversal rev;

    rev.begin    = min(h1.S, h2.F);
    rev.end      = max(h1.S, h2.F);
    rev.cost_function = this->param["COST_FUNCTION"];
    rev.calcCost();

    return rev;
}

Reversal ExactUnitary::cutHurdle(Permutation &perm, pii h1)
{
    Reversal rev;

    int j = h1.F;
    ll e = perm.at(j) + 1;
// cout<<"e = "<<e<<endl;
    int k = perm.position(e);

// cout<<"CUT "<<j <<" "<<k<<endl;

    rev.begin = min(j, k) + 1;
    rev.end = max(j, k) - 1;
    rev.cost_function = this->param["COST_FUNCTION"];
    rev.calcCost();

    return rev;
}

int ExactUnitary::getNumHurdlesWhenUnsigned(Permutation &perm)
{
    Reversal rev;
    vector<pii> hurdles;
    bool asUnsigned = this->param["UNSIGNED"];
    // cout<<"asUnsigned "<<asUnsigned<<endl;
    while (perm.isSorted(1, perm.size(), asUnsigned) == false) {
        rev = this->getRevWithMaxOrientedPairs(perm);

        if (rev.begin != -1) {
            perm.reverse(rev);
        } else {
            break;
        }
    }

    if (perm.isSorted(1, perm.size()), asUnsigned)
        return 0;

    this->findFramedIntervals(perm, hurdles);
    return hurdles.size();
}

pii ExactUnitary::findSimpleHurdle(Permutation &perm, vector<pii> &hurdles)
{
    size_t i;
    int orig_num_hurdles = hurdles.size();
    int min_hurdles = orig_num_hurdles;
    pii hurdle = make_pair(-1, -1);
    for (i = 0; i < hurdles.size(); i++) {
        Permutation new_perm = perm;
        Reversal rev = this->cutHurdle(new_perm, hurdles[i]);
        new_perm.reverse(rev);
        int new_num_hurdles = this->getNumHurdlesWhenUnsigned(new_perm);

        if (new_num_hurdles < min_hurdles) {
            new_num_hurdles = min_hurdles;
            hurdle = hurdles[i];
        }
    }
    return hurdle;
}

Result ExactUnitary::sort(Permutation &perm)
{
  // Muda dependendo da linha de comando.
    // Select sorting method
  // Isso aqui pega a solucao inicial
  //cout << "estou aqui??" << endl;
    if (perm.getParameter("signed")) {
      //cout << "entrando em sort by bergeron?" << endl;
        return this->sortByBergeron(perm);
    } else {
        return this->sortUnsigned(perm);
    }
}

Result ExactUnitary::sortUnsigned(Permutation &perm)
{
    Result result;
    perm.setExtendedMode();
    // The reduced permutation is used for speed improvement
    while (perm.isSorted(1, perm.size()) == false) {
        Reversal rev = this->getUnsignedRev(perm);
        perm.reverse(rev);
        WATCH(rev);
        WATCH(perm);
        // Insert the current (original version of) reversal in the result's end
        //Reversal orig_rev = perm.convertReversal(rev);
        result.reversals.push_back(rev);
    }
    // Alterado por Ulisses
    result = Solver::calcTotalCost(result, this->param["COST_FUNCTION"]);
    return result;
}

Result ExactUnitary::sortByBergeron(Permutation &perm)
{
    Reversal rev;
    Result result;
    perm.setParameter("reduced", true);
    perm.setParameter("extended", true);
    perm.reduce();

    // FOR DEBUG
    map<Permutation, bool, PermutationComparer > isRepeted;
    isRepeted[perm] = true;
    //cout<<"Sort by Bergeron "<<perm<<endl;
    PRINT("INITIAL PERM\n");
    WATCH(perm);

    while (perm.isSorted(1, perm.size()) == false) {
        WATCH(perm);

	// O nome do metodo nos leva ao erro, ele nao esta pegando a
	// reversao com o maximo numero de orientedpairs, ele estah
	// fazendo todo o processo de selecao usando a roleta nesta
	// funcao!!!
        rev = this->getRevWithMaxOrientedPairs(perm);
        if (rev.begin == -1) { /* no oriented reversal has not been found */
            vector<pair<int, int> > hurdles;
            this->findFramedIntervals(perm, hurdles);
            int num_hurdles = hurdles.size();


            if (num_hurdles % 2 == 0) {
                if (num_hurdles == 2)
                    rev = this->mergeHurdles(hurdles[0], hurdles[1]);
                else { /* Merges any non consecutive hurdles */
                    /* TODO: verify if this part of code is correct
                     * What really mean non consecutive???*/
                    rev = this->mergeHurdles(hurdles[0], hurdles[2]);
                }

            } else {
                if (this->param["GRIMMRevs"]) {
                    cerr<<string("ERROR: Hurdles are not supposed to be treat by Bergeron")
                        + string(" algorithm when GRIMM reversals are used.")<<endl;
                    exit(1);
                }

                if (num_hurdles == 1) {
                    rev = this->cutHurdle(perm, hurdles[0]);
                } else {
                    pii hurdle = this->findSimpleHurdle(perm, hurdles);
                    if (hurdle.F != -1) { /* Has a simple hurdle */
                        rev = this->cutHurdle(perm, hurdle);
                    } else if (num_hurdles == 3) {
                        /* A reversal that merge any two hurdles */
                        rev = this->mergeHurdles(hurdles[0], hurdles[1]);
                    } else {
                        /* A reversal that merge any two non-consecutive hurdles */
                        /* TODO: verify if this part of code is correct
                          * What really mean non-consecutive???*/
                        rev = this->mergeHurdles(hurdles[0], hurdles[2]);
                    }
                }
            }
        }


        /* Add original reversal */
        result.reversals.push_back(perm.convertReversal(rev));

        /* Peforms reversal */
        perm.reverse(rev);
    }
    // Alterado por Ulisses

    result = Solver::calcTotalCost(result, this->param["COST_FUNCTION"]);

    return result;
}

void ExactUnitary::getAllRev(
    Permutation &perm,
    vector<Reversal> &revs_orig,
    vector<Reversal> &revs_reduced,
    vector<ll> &nops
)
{
    int i, j;
    int n = perm.size();

    // Get number of breakpoints emulating perm as unsigned
    ll bk_num = perm.numBreakPoints(true);
    bool asUnsigned = this->param["UNSIGNED"];
    //cout<<"asUnsigned "<<asUnsigned<<endl;
    for (i = 1; i < n; i++) {
        for (j = i; j < n; j++) {
	  Reversal rev(i, j, 0, this->param["COST_FUNCTION"]);
	  rev.calcCost();
	  Reversal rev_orig = perm.convertReversal(rev);
	  Permutation p2 = perm;

	  // performs the revarsal
	  p2.reverse(rev);

	  // Get number of breakpoints emulating p2 as unsigned
	  ll new_bk_num = p2.numBreakPoints(true);

	  if (asUnsigned && (rev.cost == 1)) {
	    //cout<<rev<<endl;
	    continue;
	  }

	  /* Just add reversals that have less or equal number of breakpoints */
	  if (new_bk_num > bk_num)
	    continue;

	  //cout<<"ADD "<<rev<<endl;
	  // add rev
	  revs_reduced.push_back(rev);

	  // add orig rev
	  revs_orig.push_back(rev_orig);

	  // get and add num of oriented pairs in p2 after performing the rev
	  ll nop = p2.getNumOrientedPairs();
	  nops.push_back(nop);
        }
    }
}

void ExactUnitary::getBreakpointRevs(
    Permutation &perm,
    vector<Reversal> &revs_orig,
    vector<Reversal> &revs_reduced,
    vector<ll> &nops
)
{
    size_t i;

    // Get number of breakpoints emulating perm as unsigned
    ll bk_num = perm.numBreakPoints(true);

    PRINT("BEGIN BreakpointsRevs\n");
    PRINT("\t");
    WATCH(bk_num);

    // Get the breakpoints reversals list
    vector<Permutation::perm_rev> bkr;
    perm.getBreakpointReversals(bkr, this->param["COST_FUNCTION"]);

    //cout<<"asUnsigned "<<asUnsigned<<endl;
    for (i = 0; i < bkr.size(); i++) {
        Reversal rev = bkr[i].S;
        Reversal rev_orig = perm.convertReversal(rev);
        Permutation p2 = perm;

        PRINT("\t");
        WATCH(rev);
//    PRINT("\t"); WATCH(new_bk_num);
        PRINT("\n");

        // performs the revarsal
        p2.reverse(rev);

        // Get number of breakpoints emulating p2 as unsigned
        ll new_bk_num = p2.numBreakPoints(true);

        /* Just add reversals that have less or equal number of breakpoints */
        if (new_bk_num > bk_num)
            continue;

        //cout<<"ADD "<<rev<<endl;
        // add rev
        revs_reduced.push_back(rev);

        // add orig rev
        revs_orig.push_back(rev_orig);

        // get and add num of oriented pairs in p2 after performing the rev
        ll nop = p2.getNumOrientedPairs();
        nops.push_back(nop);
    }
    WATCH(nops.size());
    PRINT("END ALL REV\n");
}

/* Parameters:
 * @perm: a reduced perm
 *TODO: It was found that GRIMM generates some reversals that
 *increase the number of breakpoints. The simplest way to handle it
 *is just ignore these reversals, but in other hand these reversals can be
 *kinda of useful.  Would be necessary some deeper changes in this codebase
 *in order to allow breakpoint increasing reversals, because the
 *Bergeron algorithm (that is implemented here) works with reduced permutation.
 */
void ExactUnitary::getGRIMMRev(Permutation &perm, vector<Reversal> &revs_orig,
                               vector<Reversal> &revs_reduced,
                               vector<ll> &nops)
{
    // Assert the input
    //cout << perm.getParameter("reduced") << perm.getParameter("extended") << endl;
    assert(perm.getParameter("reduced") && perm.getParameter("extended"));

    // Makes a non reduced copy of the permutation
    Permutation non_reduced_perm = perm;



    non_reduced_perm.setNormalMode();

    // cout << "perm: " << non_reduced_perm.toString() << endl;


    // Get a set of reversals from grimm

    vector<Reversal> grimm_revs;

    getAllOptRev(non_reduced_perm, grimm_revs, this->param["COST_FUNCTION"]);    


    // Get the number of breakpoints of the non recuced permutation
    ll bk_num = non_reduced_perm.numBreakPoints();

    for (size_t i = 0; i < grimm_revs.size(); i++) {
        Reversal rev = grimm_revs[i];
        /* Select valid GRIMM reversals, i.e, those ones that don't increase
         the number of breakpoints
         */
        try {
            // Get the reduced reversal
            Reversal rev_reduced = perm.convertReversalToReduced(rev);

            /* Create a local copy of the non reduced permutation and perform the
             * inversion
             */
            Permutation non_reduced_perm_copy = non_reduced_perm;
            non_reduced_perm_copy.reverse(rev);
            ll new_bk_num = non_reduced_perm_copy.numBreakPoints();

            if (new_bk_num > bk_num) {
                continue;
            }

            // Store reversals
            revs_reduced.push_back(rev_reduced);
            revs_orig.push_back(rev);

            // Get and add num of oriented pairs in p2 after performing the rev
            ll nop = non_reduced_perm_copy.getNumOrientedPairs();
            nops.push_back(nop);
        } catch (string &exp) {
            /* found a reversal that increase number of breakpoints, i.e,
             * breaks any strip.
             */
            WATCH(exp);
            if (exp == "CANNOT_REDUCE") {
                continue;
            } else {
                cerr<<"UNEXPECTED EXCEPTION at ExactUnitary.cpp"<<endl;
                exit(1);
            }
        }
    }

    assert (revs_reduced.size() == revs_orig.size() && revs_orig.size() == nops.size());
    assert (grimm_revs.size() > 0 && revs_orig.size() > 0);
}

/* Ulisses: Esse cohdigo foi baseado no getGRIMMRev. Eu tentei fazer o
 menor numero de modificacoes para rodar os testes antes do
 paper. Entretanto, eu acredito sinceramente que ambas as abordagens
 precisam ser mudadas com urgencia. Estou tentando entender porque
 estamos usando a permutacao reduzida.
 */
void ExactUnitary::getSiepelRev(Permutation &perm, vector<Reversal> &revs_orig,
                               vector<Reversal> &revs_reduced,
                               vector<ll> &nops)
{
    // Assert the input
    //cout << perm.getParameter("reduced") << perm.getParameter("extended") << endl;
    assert(perm.getParameter("reduced") && perm.getParameter("extended"));

    // Makes a non reduced copy of the permutation
    Permutation non_reduced_perm = perm;
    non_reduced_perm.setNormalMode();

    // Get a set of reversals from grimm

    vector<Reversal> grimm_revs;

    // ##################################################
    // ##################################################
    // ##################################################
    string perm_string = non_reduced_perm.toString();
    // cout << "perm: " << perm_string << endl;


    redi::ipstream proc("~/run_siepel.sh " + perm_string);
    string line;

    while (getline(proc.out(), line)){
      string sub;
      istringstream iss(line);
      iss >> sub;
      int i = atoi(sub.c_str());
      iss >> sub;
      int j = atoi(sub.c_str());

      //cout << i <<  "," << j << endl;
      Reversal rev(i, j, 0, this->param["COST_FUNCTION"]);
      grimm_revs.push_back(rev);
    }
    //cout << endl;
    // ##################################################
    // ##################################################
    // ##################################################


    // getAllOptRev(non_reduced_perm, grimm_revs, this->param["COST_FUNCTION"]);



    // Get the number of breakpoints of the non recuced permutation
    ll bk_num = non_reduced_perm.numBreakPoints();

    for (size_t i = 0; i < grimm_revs.size(); i++) {
        Reversal rev = grimm_revs[i];
        /* Select valid GRIMM reversals, i.e, those ones that don't increase
         the number of breakpoints
         */
        try {
            // Get the reduced reversal
            Reversal rev_reduced = perm.convertReversalToReduced(rev);

            /* Create a local copy of the non reduced permutation and perform the
             * inversion
             */
            Permutation non_reduced_perm_copy = non_reduced_perm;
            non_reduced_perm_copy.reverse(rev);
            ll new_bk_num = non_reduced_perm_copy.numBreakPoints();

            if (new_bk_num > bk_num) {
                continue;
            }

            // Store reversals
            revs_reduced.push_back(rev_reduced);
            revs_orig.push_back(rev);

            // Get and add num of oriented pairs in p2 after performing the rev
            ll nop = non_reduced_perm_copy.getNumOrientedPairs();
            nops.push_back(nop);
        } catch (string &exp) {
            /* found a reversal that increase number of breakpoints, i.e,
             * breaks any strip.
             */
            WATCH(exp);
            if (exp == "CANNOT_REDUCE") {
                continue;
            } else {
                cerr<<"UNEXPECTED EXCEPTION at ExactUnitary.cpp"<<endl;
                exit(1);
            }
        }
    }

    assert (revs_reduced.size() == revs_orig.size() && revs_orig.size() == nops.size());
    assert (grimm_revs.size() > 0 && revs_orig.size() > 0);
}


struct ReversalCmp {
    bool operator() (Reversal a, Reversal b)
    {
        if (a.begin < b.begin)
            return true;
        else if (a.begin == b.begin) {
            return a.end < b.end;
        }
        return false;
    }
};

/* Get only oriented reversals. */
void ExactUnitary::getOrientedRev(Permutation &perm, vector<Reversal> &revs_orig,
                                  vector<Reversal> &revs_reduced,
                                  vector<ll> &nops)
{
    int i, j, k;
    int n = perm.size();
    //cout<<"curr perm "<<perm<<endl;

    /* *LAZY* solution for avoid repated reversals */
    //map<Reversal, bool, ReversalCmp> E;

    /* Choose best reversal under typeHeuristic used */
    for (i = 0; i < n; i++) {
        int p = perm.position(perm.at(i) + 1);

        /* Check if is a oriented pair */
        if (i < p && perm.signAt(i) * perm.signAt(p) == -1) {
            j = i;
            k = p;

            //cout<<"i "<<i<<" p "<<p<<endl;

            if (k < j)
                swap(j, k);

            /* basead in those cases :
             *  +(i) ... -(-i + 1) => \rho(j + 1, k)
             *  -(i) ... +(-i + 1) => \rho(j, k - 1)
             *  +(i + 1) ... -(i) => \rho(j, k - 1)
             *  -(i + 1) ... +(i) => \rho(j + 1, k)
             *
             * Where j is the lower index and k is the greater index.
             */
            if (perm.signAt(i) > 0)
                j++;
            else
                k--;

            /* Build Reversal */
            Reversal rev(j, k, 0, this->param["COST_FUNCTION"]);
	    rev.calcCost();

            /* Check if rev already was used */
            /* if (E[rev] == true)
               continue;

             cout<<"Add rev "<< rev<<endl;
             E[rev] = true;
             */
            /* Copy current permutation */
            Permutation p2 = perm;

            /* Get the original (from the expanded permutation) reversal */
            Reversal rev_orig = perm.convertReversal(rev);

            /* Peforms reversal */
            p2.reverse(rev);

            /* Get p2's number of oriented pairs */
            ll num_op = p2.getNumOrientedPairs();

            /* Add rev and nop number */
            revs_orig.push_back(rev_orig);
            revs_reduced.push_back(rev);
            nops.push_back(num_op);
        }
    }
}

Reversal ExactUnitary::getUnsignedRev(Permutation &perm)
{
    /* Gets permutation's size */
    int n = perm.size();

    /* Reversal in the original permutation, no reduced */
    Reversal best_rev(-1, -1, -1, this->param["COST_FUNCTION"]);

    /* Categorize reversals into ones that decrease the number of breakpoints
     * and ones that preserve the number of breakpoints. Reversals that increase
     * the number of breakpoints are descarted.
     * NOTE: The Oriented-pair concept does not make sense for unsigned permutations
     * however it is used here  in order to be compatible with the
     * HeuristicChooseReversal.
     * API (That was originally designed for signed permutations)
     */
    vector<Reversal> decresingBK;
    vector<ll> decresingBK_nops;
    vector<Reversal> preservingBK;
    vector<ll> preservingBK_nops;
    ll bk_num = perm.numBreakPoints();
    WATCH(bk_num);

    // Get the breakpoints reversals list
    vector<Permutation::perm_rev> bkr;
    perm.getBreakpointReversals(bkr, this->param["COST_FUNCTION"]);

    // Add decresing reversals
    for (size_t i = 0; i < bkr.size(); i++) {
        decresingBK.push_back(bkr[i].S);
        decresingBK_nops.push_back(1);
        PRINT("D [%d, %d]: { %lld, %lld} \n", bkr[i].S.begin, bkr[i].S.end, perm[bkr[i].S.begin], perm[bkr[i].S.end]);
    }

    //TODO: Swap to use all reversals generated by the combination of strips
    // and only make that computation when it be strictly necessarily, e.g,
    // when there is not any BK decreasing reversal
    // Add preserving reversals from the current strips
    vector<Strip> strips = perm.getAllStrips();
    for (size_t i = 0; i < strips.size(); i++) {
        Reversal rev;
        rev.begin = strips[i].begin;
        rev.end = strips[i].end;
	rev.cost_function = this->param["COST_FUNCTION"];
        if ((rev.begin == 1 && rev.begin == perm[rev.begin])
                || (rev.end == perm.size() - 1 && rev.end == perm[rev.end]))
            continue;
	rev.calcCost();
        assert(rev.begin >= 1 && rev.begin <= n);
        assert(rev.end >= 1 && rev.end <= n);
        preservingBK.push_back(rev);
        preservingBK_nops.push_back(1);
        // Validate reversal
        Permutation pp = perm;
        pp.reverse(rev);
        ll new_bk = pp.numBreakPoints();
        PRINT("P [%d, %d]: {%lld, %lld}, bk %lld\n", rev.begin, rev.end, perm[rev.begin], perm[rev.end], new_bk);
        assert( new_bk <= bk_num);
    }

    /*for (i = 1; i < n; i++) {
      for (j = i + 1; j <= n; j++) {
        Reversal rev(i, j, abs(j - i + 1));
        //Reversal rev_orig = perm.convertReversal(rev);
        Permutation p2 = perm;

        // performs the revarsal
        p2.reverse(rev);

        // Get number of breakpoints of the resulting permutation
        ll new_bk_num = p2.numBreakPoints();
        PRINT("[%lld, %lld]: %d ", rev.begin, rev.end, new_bk_num);
        WATCH(p2);
        PRINT("\n");
        if (new_bk_num < bk_num) {
          decresingBK.push_back(rev);
          decresingBK_nops.push_back(0);
        } else if (new_bk_num == bk_num) {
          preservingBK.push_back(rev);
          preservingBK_nops.push_back(0);
          //PRINT("P [%lld, %lld]: %d\n", rev.begin, rev.end, new_bk_num);
        }
      }
    }*/

    // Select the appropriate set of reversals
    vector<Reversal> &reversalsSet = decresingBK;
    vector<ll> &nopsSet = decresingBK_nops;
    //WATCH(decresingBK.size());
    //WATCH(preservingBK.size());

    if (decresingBK.size() == 0) {
        if (preservingBK.size() > 0) {
            reversalsSet = preservingBK;
            nopsSet = preservingBK_nops;
        } else {
            cerr<<"getUnsignedRev: factible reversals have not been found"<<endl;
            exit(1);
        }
    }
    /* A sequence of Heuristics are used here to try select a reversal, It keep
     * the choice of the first heuristic (in the sequence) that successful choose
     * a reversal. If any is not selected a fatal error will happen.
     */
    Permutation orig_perm = perm;
    int id_rev = -1;
    if (reversalsSet.size() != 0) {
        for (size_t h = 0; h < this->hcr_seq.size(); h++) {
            HeuristicChooseReversal *hcr = hcr_seq[h];
            id_rev = hcr->choose(orig_perm, reversalsSet, nopsSet);
            /* Break when a reversal is found */
            if (id_rev != -1) {
                break;
            }
        }
    }

    /* Verify if was chosen a reversal */
    if (id_rev == -1) {
        cerr<<"getUnsignedRev: Was no possible to choose a reversal!"<<endl;
        exit(1);
    }

    best_rev = reversalsSet[id_rev];
    return best_rev;
}

/*
 * Return a reversal with maximun score.
 * Score of a reversal is the number oriented pairs in
 * the permutation after the reversal have be performed.
 * If more the one reversal have the same maximun score
 * is choosed that have minimum length (2nd criteria).
 *
 * If have not a reversal with score greather 0, is returned
 * a invalid reversal => \rho (-1, -1, INF).
 */
Reversal ExactUnitary::getRevWithMaxOrientedPairs(Permutation &perm)
{

    /* Reversal in the original permutation, no reduced */
    Reversal orig_rev;
    Reversal best_rev(-1, -1, -1, this->param["COST_FUNCTION"]);

    /* Fond reversals in original permuation */
    vector<Reversal> revs_orig;
    /* Find reversals in reduced permutation*/
    vector<Reversal> revs_reduced;
    /* NOP of each reversal */
    vector<ll> nops;

    /* Stores perm in an auxiliar variable*/
    Permutation p2 = perm;

    PRINT("getRevWithMaxOrientedPairs\n");

    /* Get Reversals */
    if (this->param["AllRev"]) {
        this->getAllRev(p2, revs_orig, revs_reduced, nops);
    } else if (this->param["GRIMMRevs"]) {
      // Na versao final, estamos entrando aqui, quando siepel nao eh
      // chamado.
      this->getGRIMMRev(p2, revs_orig, revs_reduced, nops);
      //this->getSiepelRev(p2, revs_orig, revs_reduced, nops);
    } else if (this->param["BreakpointsRevs"]) {
        this->getBreakpointRevs(p2, revs_orig, revs_reduced, nops);
    } else {
        this->getOrientedRev(p2, revs_orig, revs_reduced, nops);
    }

    if (this->param["STATISTICS"] && nops.size() > 0) {
        // Add up the current number of reversals
        this->statistics["sum_revs"] += (double) nops.size();
        ll max_nop = *max_element(nops.begin(), nops.end());
        // Add up the current number of  Opt reversals
        ll count_max_nops = 0;
        for (size_t i = 0; i < nops.size(); i++) {
            if (nops[i] == max_nop) {
                count_max_nops ++;
            }
            //WATCH(nops[i]);
        }
        assert(count_max_nops > 0);
        this->statistics["sum_opt_revs"] += (double) count_max_nops;

        //Increase the number of choose execution
        this->statistics["num_choose_execs"] += 1.0;
        //cout<<"sum_opt_revs:"<<this->statistics["sum_opt_revs"]<<endl;
    }

    //WATCH(revs_orig.size());

    /* If perm have oriented reversal */
    if (revs_orig.size() != 0) {

        /*Gets original permutation, not reduced */
        Permutation orig_perm = perm;
        /* Look for MessType that will be used */
        if (this->param.find("MessType") != this->param.end()) {
            // cout<<"Set MessType: "<<this->param["MessType"]<<endl;
            orig_perm.setParameter("MessType", this->param["MessType"]);
        }
        orig_perm.setNormalMode();
        /* Look for a heuristic that choose a (good) reversal */
        int id_rev = -1;

        if ( this->param["MODE_RAND"]) { // Randomly mode
            int sz = hcr_seq.size();
            int h = rand() % ( sz - 1);
            HeuristicChooseReversal *hcr = hcr_seq[h];
            id_rev = hcr->choose(orig_perm, revs_orig, nops);

            //cout<<"CHOOSED "<<h<<" id_rev "<<id_rev<<endl;

            /* Not Found a good reversal */
            if (id_rev == -1) {
                HeuristicChooseReversal *hcr = hcr_seq[sz - 1];
                id_rev = hcr->choose(orig_perm, revs_orig, nops);

            }

        } else { // Sequential mode
            for (size_t h = 0; h < this->hcr_seq.size(); h++) {
                HeuristicChooseReversal *hcr = hcr_seq[h];
                id_rev = hcr->choose(orig_perm, revs_orig, nops);
                //PRINT("\tchoosed rev id %d, with h = \n\t", id_rev, h);
                /* Found a good reversal */
                if (id_rev != -1) {
                    //WATCH(revs_orig[id_rev]);
                    break;
                }
            }
        }

        /* Verify if was chosen a reversal */
        if (id_rev == -1) {
            cerr<<"getRevWithMaxOrientedPairs: Was no possible find a reversal!"<<endl;
            exit(1);
        }

        best_rev = revs_reduced[id_rev];
    }
    WATCH(best_rev);
    return best_rev;
}



/* This function works just for unsigned permutations */
void ExactUnitary::makeConsecutiveMap(Permutation &perm, vector<int> &begin,
                                      vector<int> &end)
{
    int i;
    int n = perm.size();
    begin = vector<int>(n + 1);
    end = vector<int>(n + 1);

    begin[0] = 0;
    for (i = 1; i <= n; i++) {
        if (perm.at(i) == perm.at(i - 1) + 1)
            begin[i] = begin[i - 1];
        else
            begin[i] = i;
    }

    end[n] = n;
    for (i = n - 1; i >= 0; i--) {
        if (perm.at(i) + 1 == perm.at(i + 1))
            end[i] = end[i + 1];
        else
            end[i] = i;
    }
}

/* : use reduced permutation */
void ExactUnitary::findFramedIntervals(Permutation &perm,
                                       vector<pair<int, int> > &hurdles)
{
    int i, j;
    int n = perm.size();

    vector<int> min_end(n + 1, n + 1);
    vector<int> max_begin(n + 1, 0);

    vector<int> min_end_cycle(n + 1, n + 1);
    vector<int> max_begin_cycle(n + 1, 0);

    /* Current max element */
    vector<ll> curr_max(n + 1);
    vector<ll> less_max(n + 1);
    vector<ll> count_greater(n + 1);
    vector<ll> count_less(n + 1);

    for (i = perm.size() - 1; i >= 0; i--) {
        bool can_continue = true;
        /* Current max element */
        curr_max[i] = perm.at(i);

        /* Current max element that is lower tham perm[i] */
        less_max[i] = 0;

        /* Number elements greather than perm[i] */
        count_greater[i] = 0;

        /* Number elements lower than perm[i] */
        count_less[i] = 0;

        min_end[i] = min_end[i + 1];
        max_begin[i] = max_begin[i + 1];

        /* Find linear framed intervals */
        for (j = i + 1; j <= perm.size(); j++) {
            /*  if(perm.at(j)  != perm.at(j - 1) + 1)
                consecutive = false;*/
            if (perm.at(j) < perm.at(i)) {
                less_max[i] = max(less_max[i], perm.at(j));
                count_less[i] = count_less[i]  + 1;
                /* Don't have more linar framed intervals
                 * that begin in i */
                can_continue = false;
            } else {
                if (can_continue && (perm.at(j) - perm.at(i)) == (j - i)
                        && perm.at(j) > curr_max[i]) {

                    if (j <= min_end[i]) {
                        hurdles.push_back(make_pair(i,j));
                        min_end[i] = j;
                    }
                    max_begin[i] = max(max_begin[i], i);
                }
                count_greater[i] = count_greater[i] + 1;
                curr_max[i] = max (curr_max[i], perm.at(j));
            }
        }
    }

    /* For circularity */
    for (i = perm.size() - 1; i >= 0; i--) {
        for (j = 0; j < i; j++) {
            if (perm.at(j) < perm.at(i)) {
                if (less_max[i] < perm.at(j) && (count_less[i]  == perm.at(j))
                        && (count_greater[i] == (n - perm.at(i)))) {
                    if (j <= min_end[i] && i >= max_begin[i]) {
                        hurdles.push_back(make_pair(i,j));
                    }
                }
                less_max[i] = max(less_max[i], perm.at(j));
                count_less[i] = count_less[i]  + 1;
            } else {
                count_greater[i] = count_greater[i] + 1;
            }
        }
    }
}
