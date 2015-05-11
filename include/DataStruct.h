#ifndef DATA_STRUCT
#define DATA_STRUCT 1
#include "HeaderDefault.h"
template<class T>
class BIT
{
    T *tree;
    int maxVal;
public:
    BIT(int N)
    {
        tree = new T[N+1];
        memset(tree, 0, sizeof(T) * (N+1));
        maxVal = N;
    }
    void update(int idx, T val)
    {
        while (idx <= maxVal) {
            tree[idx] += val;
            idx += (idx & -idx);
        }
    }
    /* Returns the cumulative frequency of index idx */
    T read(int idx)
    {
        T sum=0;
        while (idx>0) {
            sum += tree[idx];
            idx -= (idx & -idx);
        }
        return sum;
    }

    ~BIT()
    {
        delete tree;
    }
};

#endif
