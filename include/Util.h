
#ifndef UTIL_HEADER
#define UTIL_HEADER
#include "HeaderDefault.h"
#include "Permutation.h"
#include "Solver.h"

struct PSet {
    map<string,PSet> childs;
    string value;
};


typedef PSet NSet;
typedef unsigned long long timestamp_t;

class Util
{

public:
    static string toString(double d);
    static string IntToString(ll d);
    static string removeExtremsSpaces(string s);
    static string readTextFile(string fname);
    static void readLine(istream &stream, string &line);
    static double lg(double d);

    static string extractDirPath(string path);
    static void mkdir(string dir);
    static void cd(string dir);
    static NSet parserSetStream(istream &stream);
    static void cp(string from,string to);
    static void mkAllDir(string path);
    static vector<string> split(string s,char sep);
    static int toInt(string s);
    static Result getResultFromStr(const char *s, int cost_function);
    static vector<ll> getSeqFromStr(const char *s);

    static ll sign(ll v);

    static vector<ll> arrayToVector(int *arr, int sz);

    // Time stuff
    /* Get passed time in microseconds */
    static timestamp_t getTimestamp ();

    /* Get the interval time (in seconds) */
    static double calcIntervalTime(timestamp_t t0, timestamp_t t1);
};

#endif
