#ifndef INSTANCE_GROUP_HEADER
#define INSTANCE_GROUP_HEADER 1
#include "HeaderDefault.h"
#include "Permutation.h"


enum instanceStorageMode {MEMORY_MODE, FILE_MODE};

class InstanceGroup
{

    vector<Permutation> permutations;
    string groupName;
    fstream *fin;
    string fileName;
    instanceStorageMode mode;
    ll N;

    Permutation currentPermutation;
    int currentIndex;
    parameter param;

public:
    InstanceGroup() {};
    InstanceGroup(parameter param, string groupName, string fileName);
    InstanceGroup(string groupName, vector<Permutation> permutations);

    string getName();
    bool next();
    Permutation get();

    ll size();
private:
};
#endif
