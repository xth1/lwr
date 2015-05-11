#include "Parameter.h"
/*
psd getParameter(string s, string &type) {
  vector<string> vs;
  vs = Util::split(s, '.');
  type = vs[0];

  vs = Util::split(vs[1], ':');
  psd p;
  p.F = vs[0];
  sscanf(vs[1].c_str(),"%lf", &p.S);

  return p;
}

map<string, Parameter> Parameter::readParameters(
    char *param_list,
    int begin,
    int size
) {

    for (int cur_p = begin; curr_p < size; curr_p++) {
        psd p = getParameter(string(argv[p_read]), type);


        if (type == "HCR") {
            hcrP[p.F] = p.S;
        } else if (type == "IMP") {
          improveP[p.F] = p.S;
        } else if (type == "BUILD_SOLUTION") {
          exactP[p.F] = p.S;
        }
    }

}
*/

bool Parameter::hasString(string key)
{
    if (string_map.find(key) != string_map.end()) {
        return true;
    }
    return false;
}

bool Parameter::hasNum(string key)
{
    if (num_map.find(key) != num_map.end()) {
        return true;
    }
    return false;
}

void Parameter::setString(string key, string val)
{
    this->string_map[key] = val;
}

void Parameter::setNum(string key, double val)
{
    this->num_map[key] = val;
}

string Parameter::getString(string key)
{
    if (!hasString(key)) {
        cout<<key + "is not a valid string parameter"<<endl;
        throw key + "is not a valid string parameter";
    }
    return this->string_map[key];
}

double Parameter::getNum(string key)
{
    if (!hasNum(key)) {
        cout<<key + "is not a valid num parameter"<<endl;
        throw key + "is not a valid num parameter";
    }
    return this->num_map[key];
}