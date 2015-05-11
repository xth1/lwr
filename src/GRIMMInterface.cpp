#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <stddef.h>   /* for offsetof */
#include <unistd.h>
#include <string.h>

#include "mcstructs.h"
#include "uniinvdist.h"
#include "mcrdist.h"
#include "mcread_input.h"
#include "write_data.h"
#include "scenario.h"
#include "testrev.h"
#include "opt_scenario.h"
#include "unsigned.h"
#include "e_malloc.h"
#include "circ_align.h"
#include "texgraph.h"
#include "ext_function.h"
#include "countperms.h"
#include "unsignedhc.h"
#include "matrixmisc.h"

#include "GRIMMInterface.h"

//using namespace std;

#define TRACE(x...)
#define PRINT(x...) TRACE(printf(x));fflush(stdout)
#define WATCH(x) TRACE(cout << #x" = " << x << "\n")

// Reference for the variable definition at mcstructs.h
FILE *outfile;

void initGRIMM()
{
    outfile = stdout;
}


/* Converts Permution to array */
int *to_array(Permutation &perm)
{
    int *p;
    int i;
    int n = (int) perm.size();
    p = (int *)calloc(n, sizeof(int));
    PRINT("PERM\n\t");
    for (i = 0; i < n; i++) {
        p[i] = perm[i + 1];
        PRINT("%d : %d\n", i, p[i]);
    }

    return p;
}

/* Create identity permutation */
int *make_id(int n)
{
    int *p;
    int i;
    p = (int *)calloc(n, sizeof(int));
    PRINT("ID\n\t");
    for (i = 0; i < n; i++) {
        p[i] = i + 1;
        PRINT("%d : %d\n", i, p[i]);
    }

    return p;
}

/*
 * Given a permutation perm, this function returns the first
 * reversal provided by Exact Sorting by signed inversions algorithm
 */
void getAllOptRev(
    Permutation &perm,
    vector<Reversal> &vr,
    int cost_function,
    int scenario_type
)
{
    struct genome_struct *g1, *g2;
    /* Some parameters for sorting by reversals */
    int num_chromosomes = 0;

    int n = (int) perm.size();

    /* Create 'genomes' for each permutation */
    g1 = (struct genome_struct *) malloc(sizeof(struct genome_struct));
    g2 = (struct genome_struct *) malloc(sizeof(struct genome_struct));
    WATCH(g1);
    WATCH(g2);

    /* Current permuation */
    g1->genes = to_array(perm);
    g1->genome_num = 1;
    g1->encoding = NULL;
    g1->gnamePtr = NULL;

    /* Identity permuation */
    g2->genes = make_id(n);
    g2->genome_num = 2;
    g2->encoding = NULL;
    g2->gnamePtr = NULL;
    /* Alloc reversals vector for the maximun possible number of elements */
    int sz = perm.size() * perm.size();
    int **revs = (int **) malloc(sizeof(int *) * sz);
    for (int i = 0; i < sz; i++) {
        revs[i] = (int *) malloc(sizeof(int) * 2);
    }

    /* Cath next reversal */
    int revs_size = 0;
    next_reversal(g1, g2, n, num_chromosomes,
                  scenario_type,
                  revs, &revs_size);

    /* Build the vetor with all reversals*/
    vr = vector<Reversal>(revs_size);
    for (int i = 0; i < revs_size; i++) {
        revs[i][0]++;
        revs[i][1]++;
        Reversal rev = Reversal(revs[i][0], revs[i][1], 0, cost_function);
	// cout << revs[i][0] << "," << revs[i][1] << endl;
	rev.calcCost();
        vr[i] = rev;
    }
    // cout << endl;
// Get a random reversal

//fred the resources
    for (int i = 0; i < sz; i++) {
        free(revs[i]);
    }
    free(revs);
    free(g1->genes);
    free(g1);
    free(g2->genes);
    free(g2);
}

Permutation assignSigns(Permutation p1, Permutation p2, int num_it) {
  struct genome_struct *genome_list, *out_genome;

  int n = (int) p1.size();

  /* Create 'genomes' for each permutation */
  genome_list = (struct genome_struct *) malloc(2 * sizeof(struct genome_struct));
  out_genome = (struct genome_struct *) malloc( sizeof(struct genome_struct));

  /* Identity permuation */
  genome_list[0].genes = to_array(p2);
  genome_list[0].genome_num = 2;
  genome_list[0].encoding = NULL;
  genome_list[0].gnamePtr = NULL;

  /* Current permuation */
  genome_list[1].genes = to_array(p1);
  genome_list[1].genome_num = 1;
  genome_list[1].encoding = NULL;
  genome_list[1].gnamePtr = NULL;


  getUnsignedSignsAssignment(
    n,
    genome_list,
    0,
    1,
    num_it,
    out_genome
  );

  // Build the output permutation
  vector<ll> v(n);

  for (int i = 0; i < n; i++) {
    v[i] = out_genome->genes[i];
  }

  Permutation out_perm(p1.getParameter(), v);

  // Free the resource
  for (int i = 0; i < 2; i++) {
    free(genome_list[i].genes);
  }
  free(genome_list);
  free(out_genome->genes);
  free(out_genome);

  return out_perm;
}


Reversal getNextOptRev(Permutation &p1, Permutation &p2, int scenario_type, int cost_function)
{
    struct genome_struct *g1, *g2;
    /* Some parameters for sorting by reversals */
    int num_chromosomes = 0;

    /* End and begin */
    int end, begin;
    int n = (int) p1.size();
    Reversal rev;

    /* Create 'genomes' for each permutation */
    g1 = (struct genome_struct *) malloc(sizeof(struct genome_struct));
    g2 = (struct genome_struct *) malloc(sizeof(struct genome_struct));
    WATCH(g1);
    WATCH(g2);

    /* First permuation */
    g1->genes = to_array(p1);
    g1->genome_num = 1;
    g1->encoding = NULL;
    g1->gnamePtr = NULL;

    /* Second permuation */
    g2->genes = to_array(p2);
    g2->genome_num = 2;
    g2->encoding = NULL;
    g2->gnamePtr = NULL;
    //printf(" Scenario type %d\n", scenario_type);
    /* Cath next reversal */
    next_opt_reversal(g1, g2, n , num_chromosomes,
                      scenario_type, &begin, &end);

    WATCH(begin);
    WATCH(end);

    free(g1->genes);
    free(g1);
    free(g2->genes);
    free(g2);

    rev = Reversal(begin + 1, end + 1, 0, cost_function);

    return rev;
}

/*
 * Given a permutation perm, this function returns the
 * min_cost reversal provided by Exact Sorting by signed inversions algorithm
 */
Reversal getMinCostOptRev(Permutation &perm, int cost_function, int scenario_type)
{
    vector<Reversal> vr;
    getAllOptRev(perm, vr, cost_function, scenario_type);

    int min_cost_id = 0;

    for (size_t i = 1; i < vr.size(); i++) {
        if (vr[i].cost < vr[min_cost_id].cost)
            min_cost_id = i;
    }

    return vr[min_cost_id];
}

// Alterado por Ulisses
Result sortByGRIMMUnsigned(Permutation &p, int num_it, bool reversed, int cost_function) {
  initGRIMM();
  Result result;
  Permutation p1 = p;
  Permutation p2 = Permutation::makeIdentity(p.getParameter(), p.size());
  p1 = assignSigns(p1, p2, num_it);

  p1.setParameter("signed", true);
  // Alterado por Ulisses
  result = sortByGRIMM(p1, reversed, cost_function);
  p1.setParameter("signed", false);
  result.removeUnitaryReversals();

  return result;
}

// Alterado por Ulisses
Result sortByGRIMM(Permutation &p, bool reversed, int cost_function)
{
    Permutation p1 = p;
    Permutation p2 = Permutation::makeIdentity(p.getParameter(), p.size());
    Result res;
    initGRIMM();

    if (reversed)
        swap(p1, p2);

    while (p1 != p2) {
      Reversal rev = getNextOptRev(p1, p2, 4, cost_function);
        res.reversals.push_back(rev);
        p1.reverse(rev);
        WATCH(rev);
        WATCH(p);
    }

    if (reversed)
        reverse(res.reversals.begin(), res.reversals.end());
    //Alterado por Ulisses
    res = Solver::calcTotalCost(res, cost_function);
    return res;
}
