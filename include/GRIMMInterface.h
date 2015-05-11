#ifndef GRIMM_INTERFACE
#define GRIMM_INTERFACE
#include "HeaderDefault.h"
#include "Permutation.h"
#include "Solver.h"

Reversal getNextOptRev(Permutation &p1, Permutation &p2, int scenario_type, int cost_function);
// Alterado por Ulisses para adicionar cost_function
void getAllOptRev(Permutation &perm, vector<Reversal> &vr, int cost_function, int scenario_type = 4);
void initGRIMM();
Reversal getMinCostOptRev(Permutation &perm, int cost_function, int scenario_type);
// Alterado por Ulisses
Result sortByGRIMM(Permutation &p, bool reversed, int cost_function);
Result sortByGRIMMUnsigned(Permutation &p, int num_it, bool reversed, int cost_function);
#endif
