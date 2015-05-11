
#include "InstanceGroup.h"

ll InstanceGroup::size()
{
    return N;
}
string InstanceGroup::getName()
{

    return this->groupName;
}

InstanceGroup::InstanceGroup(parameter param,string groupName, string fileName)
{
    this->groupName = groupName;
    this->fileName = fileName;
    //this->
    this->mode = FILE_MODE;
    this->param = param;
    this->currentIndex = 0;

    fin = new fstream();

    fin->open(this->fileName.c_str(),ios::in);

    *(fin)>>(this->N);
}
InstanceGroup::InstanceGroup(string groupName, vector<Permutation> vp)
{
    this->groupName = groupName;

    this->permutations = vp;
    this->mode = MEMORY_MODE;

    this->N = (ll)vp.size();

    //cout<<"N"<<N<<endl;

    this->currentIndex = 0;

    //this->param = param;
}

bool InstanceGroup::next()
{
    if ((ll)currentIndex == N)
        return false;
    if (mode == FILE_MODE) {
        //ifstream fin(this->fileName.c_str());
        Permutation p(this->param);
        *(fin)>>p;
        this->currentPermutation = p;

    } else {
        //cout<<"curr "<<currentIndex<<" N "<<N<<endl;
        this->currentPermutation = this->permutations[currentIndex];
    }

    currentIndex++;
    return true;
}

Permutation InstanceGroup::get()
{

    return this->currentPermutation;
}
