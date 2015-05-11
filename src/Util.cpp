#include "Util.h"

double Util::lg(double d)
{

    return log(d)/log(2.0);
}

string Util::toString(double d)
{
    char buff[100];

    sprintf(buff,"%lf",d);

    return string(buff);
}


string Util::IntToString(ll d)
{
    char buff[100];

    sprintf(buff,"%lld",d);

    return string(buff);
}

int Util::toInt(string s)
{
    int d;
    sscanf(s.c_str(),"%d",&d);

    return d;
}

string Util::removeExtremsSpaces(string s)
{
    int begin,end;

    begin = 0;
    int sz = s.size();
    while (begin < sz && s[begin] ==' ')
        begin++;

    end = sz - 1;
    while (end >0 && s[end] ==' ')
        end--;

    string r = s.substr(begin,end - begin + 1);

    //cout<<"R "<<r<<" "<<begin<<" "<<end<<endl;

    return r;

}

void Util::readLine(istream &stream, string &line)
{
    cout<<"Read line "<<endl;
    getline(stream,line);

    cout<<"Read line "<<line<<endl;
    line = Util::removeExtremsSpaces(line);
}

string Util::extractDirPath(string path)
{

    int sz = path.size();
    int i;
    for (i = sz - 1; i >=0; i--) {
        if (path[i] == '/')
            break;
    }

    string r = path.substr(0,i);

    return r;

}

string Util::readTextFile(string fname)
{
    ifstream fin(fname.c_str());

    string doc = "";
    while (1) {
        string line;
        getline(fin, line);

        if (fin.eof())
            break;
        doc += line + "\n";
    }

    return doc;
}

void Util::mkdir(string dir)
{
    string cmd = "mkdir "+dir;
    system(cmd.c_str());
}
void Util::cd(string dir)
{
    string cmd = "cd "+dir;
    system(cmd.c_str());
}

void Util::cp(string from,string to)
{
    string cmd = "cp " + from +" "+to;
    system(cmd.c_str());
}

void Util::mkAllDir(string path)
{
    for (size_t i = 0 ; i < path.size(); i++) {
        if (path[i] == '/') {
            string sub_path = path.substr(0,i);
            mkdir(sub_path);
        }
    }
    mkdir(path);
}

vector<string> Util::split(string s,char sep)
{
    vector<string> v;
    size_t ini = 0;

    for (size_t i = 0; i < s.size(); i++) {
        if (s[i] == sep) {
            v.push_back(s.substr(ini,i-ini));
            ini = i+1;
        }
    }

    if (ini < s.size())
        v.push_back(s.substr(ini,s.size() - ini));
    return v;
}

ll Util::sign(ll v)
{
    if (v > 0)
        return 1;
    else if ( v < 0)
        return -1;
    return 0;
}

vector<ll> Util::arrayToVector(int *arr, int sz)
{
    vector<ll> v(sz);
    for (int i = 0; i < sz; i++) {
        v[i] = arr[i];
    }
    return v;
}

Result Util::getResultFromStr(const char *s, int cost_function)
{
    string ss = string(s + 2);
    vector<string> rs = Util::split(ss, '[');
    vector<Reversal> revs;
    for (size_t i = 0; i < rs.size(); i++) {
        Reversal rev;
        sscanf(rs[i].c_str(), "%lld,%lld", &rev.begin, &rev.end);
        rev.cost_function = cost_function;
	rev.calcCost();
        // cout<<"rev "<<rev<<endl;
        revs.push_back(rev);
    }
    Result r;
    r.reversals = revs;
    r = Solver::calcTotalCost(r, cost_function);
    return r;
}

vector<ll> Util::getSeqFromStr(const char *s)
{
    vector<string> vs = Util::split(string(s), ',');

    vector<ll> v;
    for (size_t i = 0; i < vs.size(); i++) {
        ll x;
        sscanf(vs[i].c_str(),"%lld",&x);
        v.push_back(x);
    }
    return v;
}

/* Get passed time in microseconds */
timestamp_t Util::getTimestamp ()
{
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t) now.tv_sec * 1000000;
}

double Util::calcIntervalTime(timestamp_t t0, timestamp_t t1)
{
    double time_sec = (t1 - t0) / 1000000.0L;
    return time_sec;
}
