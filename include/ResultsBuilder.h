#ifndef RESULTS_BUILDER_H
#define RESULTS_BUILDER_H
#include "HeaderDefault.h"
#include "Permutation.h"
#include "SequentialDataAnalyser.h"
#include "Util.h"
#define MAXBUFF 40096

typedef struct {
    string file_name;
    string name;
} Variation;

/* Store variations results of the same size */
typedef struct {
    int size;
    /* File name of variations results */
    vector<Variation> variations;
} SizePerm;

typedef struct {
    ll cost;
    ll rev;
    string seq;
} ResultLine;

class ResultsBuilder
{
public:

    ResultsBuilder (int num_sizes, int num_variations,
                    vector<SizePerm> file_names);
    /* Configuration */


    /* Methods for result files */
    void read_result_line(ifstream &fin, ResultLine &result_line);
    void write_result_line(ofstream &fout, ResultLine &result_line);
    void buildCSV(vector<Table> tables, string out_file_name);

    Table computeTableBySize(vector<Variation> &variations, string out_fname);

private:

    /* Matrix to store results file names:
     * - file_names[n][k]
     * Where n is size of permutations and k is a variation name
     */
    vector<SizePerm> file_names;

    /* Number of permutation sizes there be used */
    int num_sizes;

    /* Number of algorithms variations there be used */
    int num_variations;

};

#endif
