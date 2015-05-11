#include "SequentialDataAnalyser.h"

SequentialDataAnalyser::SequentialDataAnalyser()
{
    quantIns = 0;
}

void SequentialDataAnalyser::addIns(InstanceResult &ir)
{

    Field minField, maxField;

    quantIns++;
    for (size_t i = 0; i < ir.size(); i++) {
        string sn = ir[i].name;

        this->sourcesName.insert(sn);
        Field::iterator it;

        for (it = (ir[i].field).begin(); it != (ir[i].field).end(); it++) {
            //MAX
            if (sourceMax[sn].find(it->F) == sourceMax[sn].end())
                sourceMax[sn][it->F] = it->S;
            else
                sourceMax[sn][it->F] = max(sourceMax[sn][it->F],it->S);

            //MIN
            if (sourceMin[sn].find(it->F) == sourceMin[sn].end())
                sourceMin[sn][it->F] = it->S;
            else
                sourceMin[sn][it->F] = min(sourceMin[sn][it->F],it->S);
            //SUM
            sourceSum[sn][it->F] += it->S;


            if (minField.find(it->F) == minField.end())
                minField[it->F] = it->S;
            else
                minField[it->F] = min(minField[it->F], it->S);

            if (maxField.find(it->F) == maxField.end())
                maxField[it->F] = it->S;
            else
                maxField[it->F] = max(maxField[it->F],it->S);
        }
    }

    for (size_t i = 0; i < ir.size(); i++) {
        string sn = ir[i].name;
        Field::iterator it;
//		cerr<<"Solver "<<sn<<endl;
        for (it = (ir[i].field).begin(); it != (ir[i].field).end(); it++) {


            //Min percent
            if (sourceMinPercent[sn].find(it->F) == sourceMinPercent[sn].end())
                sourceMinPercent[sn][it->F] = 0.00000;

            if (it->S == minField[it->F]) {
                sourceMinPercent[sn][it->F]++;
            }

            //Max percent
            if (sourceMaxPercent[sn].find(it->F) == sourceMaxPercent[sn].end())
                sourceMaxPercent[sn][it->F] = 0.000000;

            if (it->S == maxField[it->F]) {
                sourceMaxPercent[sn][it->F]++;
            }
        }
    }
}

double SequentialDataAnalyser::getMin(string source, string field)
{
    return sourceMin[source][field];
}

double SequentialDataAnalyser::getMax(string source, string field)
{
    return sourceMax[source][field];
}

double SequentialDataAnalyser::getAvg(string source, string field)
{
    return sourceSum[source][field] / (double) quantIns;
}

double SequentialDataAnalyser::getMinPercent(string source, string field)
{
    return sourceMinPercent[source][field] / (double) quantIns;
}

double SequentialDataAnalyser::getMaxPercent(string source, string field)
{
    return sourceMaxPercent[source][field] / (double) quantIns;
}

Table SequentialDataAnalyser::SequentialDataAnalyser::makeTable()
{
    Table table;

    set<string>::iterator it;
    for (it = this->sourcesName.begin(); it != this->sourcesName.end(); it++) {
        OutputLine ol;
        string sn = *it;

        Field::iterator itf;

        for (itf = sourceMax[sn].begin(); itf != sourceMax[sn].end(); itf++) {
            string ss = itf->F + "_max";
            ol[ss] = getMax(sn,itf->F);
        }

        for (itf = sourceMin[sn].begin(); itf != sourceMin[sn].end(); itf++) {
            string ss = itf->F + "_min";
            ol[ss] = getMin(sn,itf->F);
        }

        for (itf = sourceMin[sn].begin(); itf != sourceMin[sn].end(); itf++) {
            string ss = itf->F + "_avg";
            ol[ss] = getAvg(sn,itf->F);
        }

        for (itf = sourceMinPercent[sn].begin(); itf != sourceMinPercent[sn].end(); itf++) {
            string ss = itf->F + "_min_percent";
            ol[ss] = getMinPercent(sn,itf->F);
        }

        for (itf = sourceMaxPercent[sn].begin(); itf != sourceMaxPercent[sn].end(); itf++) {
            string ss = itf->F + "_max_percent";
            ol[ss] = getMaxPercent(sn,itf->F);
        }

        table[sn] = ol;
    }
    return table;
}
