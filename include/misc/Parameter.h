#ifndef PARAMETER_CLASS
#define PARAMETER_CLASS 1
#include "HeaderDefault.h"
class Parameter
{
    map<string, double> num_map;
    map<string, string> string_map;

public:
    bool hasString(string key);
    bool hasNum(string key);
    void setString(string key, string val);
    void setNum(string key, double val);
    string getString(string key);
    double getNum(string key);

    static map<string, Parameter> readParameters(
        char *param_list,
        int begin,
        int size
    );
};

#endif