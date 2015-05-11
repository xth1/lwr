#ifndef SEQUENTIAL_DATA_ANALYSER
#define SEQUENTIAL_DATA_ANALYSER 1

#include "HeaderDefault.h"

typedef map<string, double> Field;

typedef struct {
    string name;
    Field field;
} SourceResult;


typedef  map<string,double> msd;

typedef msd OutputLine; // attribute (collum: value
typedef map<string, OutputLine> Table; //solver: line

typedef vector<SourceResult> InstanceResult;

class SequentialDataAnalyser
{

    set<string> sourcesName;
    map<string,Field > sourceMax;
    map<string,Field > sourceMin;
    map<string,Field > sourceSum;

    map<string,Field > sourceMinPercent;
    map<string,Field > sourceMaxPercent;

    ll quantIns;
public:
    SequentialDataAnalyser();

    void addIns(InstanceResult &ir);

    double getMin(string source, string field);
    double getMax(string source, string field);
    double getAvg(string source, string field);
    double getMinPercent(string source, string field);
    double getMaxPercent(string source, string field);

    Table makeTable();
};

#endif
