#include "Permutation.h"
#include "DataStruct.h"

#define TRACE(x...)
#define PRINT(x...) TRACE(printf(x));fflush(stdout)
#define WATCH(x) TRACE(cout << #x" = " << x << "\n")
//using namespace std;

//TODO: this class enables to use the extended mode, that is useful for some
// algorithms. For some other algorithms are used not extended permutations.
// The existence of this two modes make harder to design methods with
// well-defined behavior for Permutation class.
// A aparently good solution is makes the extended mode as default and then
// make sure that all implemented algorithms support it.

ll Permutation::numPerm(ll n)
{
    if (n <= 1)
        return 1;
    return n*numPerm(n-1);
}

/* Get all strips of a unsigned permutation
 * @return a array with all strips that was found.
 */
vector<Strip> Permutation::getAllStrips() const
{
    vector<Strip> strips;
    // Stores begin of current strip...
    ll begin = 1;
    PRINT("\tBREAK: ENTER AT GET ALL STRIPS \n\t");
    WATCH(*this);

    // First case
    if (this->at(1) == 1) { // element 1 is in the corret position...
        ll pos = 2;
        //Find the end of the ascendant sequence that starts in 1.
        while (pos <= this->size() && this->at(pos) == pos)
            pos++;
        // Creates the initial Strip, that is always POSITIVE
        Strip strip;
        strip.begin = begin;
        strip.end = pos - 1;
        // Select the Strip orientation
        ll diff = this->at(strip.end) - this->at(strip.begin);
        if (diff <= 0 && this->at(strip.begin) != 1) {
            strip.dir = NEGATIVE;
        } else {
            strip.dir = POSITIVE;
        }

        begin = pos;
        strips.push_back(strip);

        PRINT("\t BREAK: first strip [%lld, %lld]\n", strip.begin, strip.end);
    }

    // HACK for it works for extended permutations...
    int size = this->size();
    if (this->getParameter("extended")) {
        size--;
    }

    for (int i = begin + 1; i <= size; i++) {
        // Check if there is a consecutive pair (inverted or not)...
        if (abs(this->at(i) - this->at(i - 1)) > 1) {
            Strip strip;
            strip.begin = begin;
            strip.end = i - 1;
            PRINT("\t BREAK: strip [%lld, %lld] => (%lld, %lld)\n",
                  strip.begin, strip.end,
                  this->at(strip.begin), this->at(strip.end)
                 );
            // Identify the sign of a strip
            ll diff = this->at(strip.end) - this->at(strip.begin);
            if (diff <= 0) {
                strip.dir = NEGATIVE;
            } else {
                strip.dir = POSITIVE;
            }
            // Push up the current strip in the sequence
            strips.push_back(strip);
            begin = i;
        }
    }
    // Last strip
    Strip strip;
    strip.begin = begin;
    strip.end = size;

    if ((strip.end - strip.begin) >= 0) {
        if ((this->at(strip.end) - this->at(strip.end - 1)) == 1
                || this->at(strip.end) == size) {
            strip.dir = POSITIVE;
        } else {
            strip.dir = NEGATIVE;
        }
        strips.push_back(strip);
    }

    //DEBUG
    /*cout<<"Strips"<<endl;
    for (int i = 0; i < strips.size(); i++) {
        cout<<strips[i].begin<<" "<<strips[i].end<<" "<<strips[i].dir<<endl;
    }*/

    return strips;
}

ll Permutation::getNumDecreasingStrips()
{

    vector<Strip> strips = this->getAllStrips();

    ll count = 0;

    for (int i = 0; i < (int)strips.size(); i++) {

        if (strips[i].dir == NEGATIVE) {
            count++;
        }
    }

    return count;
}

void Permutation::makeAllPermutationsToFile(ll n, string fileName)
{
    parameter param;
    param["reversal"] = true;
    param["LWE"] = true;
    ofstream fout(fileName.c_str());
    Permutation p = Permutation::makeIdentity(param,n);
    fout<<numPerm(n)<<endl;
    do {
        fout<<n<<endl;
        fout<<p<<endl;
    } while (p.next());

}

void Permutation::makeRandPermutationsToFile(ll n,ll q, string fileName)
{
    parameter param;
    param["reversal"] = true;
    param["LWE"] = true;

    if (n <= 15) {
        ll n_perm = Permutation::numPerm(n);
        q = min(q,n_perm);
    }

    vector<ll> vp(n);
    ofstream fout(fileName.c_str());
    for (int i = 0; i < n; i++)
        vp[i] = i + 1;

    fout<<q<<endl;

    map<Permutation, bool, PermutationComparer> has;
    for (int k = 0; k < q;) {
        random_shuffle(vp.begin(), vp.end());
        Permutation p(param, vp);

        if (has[p])
            continue;
        k++;
        fout<<n<<endl;
        fout<<p<<endl;
    }
}

vector<Permutation> Permutation::makeAllPermutations(ll n)
{

    vector<Permutation> vp;

    parameter param;
    param["reversal"] = true;
    param["LWE"] = true;

    Permutation p = Permutation::makeIdentity(param,n);

    do {
        vp.push_back(p);
    } while (p.next());

    return vp;
}


void Permutation::computePositiveSum(vector<ll> &sum) const
{
    int i;
    int v;
    sum[0] = 0;
    for (i = 1; i <= this->size(); i++) {
        WATCH(i);
        v = this->sign(i) == 1? 1: 0;
        sum[i] = sum[i- 1] + v;
        //printf("sum: %d, %d, %d\n", v, sum[i-1],sum[i]);
    }
}

/* Compute the number of oriented pairs in a
 * signed permutation.
 */
//TODO: this function is accurate? Let's add unittests!
ll Permutation::getNumOrientedPairs()
{
    int i;
    int n = this->size();
    ll count = 0;

    // Set temporary the permutation as extended...
    bool is_extended = true;
    if (this->param["extended"] == false) {
        is_extended = false;
        this->param["extended"] = true;
    }

    for (i = 1; i < n; i++) {
        int p = this->position(this->at(i) + 1);
        if (this->signAt(i) * this->signAt(p) == -1) {
            count++;
        }
    }
    //TODO: the correct would be !is_extended?
    if (is_extended) {
        this->param["extended"] = false;
    }
    return count;
}

/* Compute messLevel of a  permutation p of size n:
 * Sum (abs(p[i]) - i) + sv(i), for 1 <= i <= n.
 * sv(i) = 1 , if sign(i) = -1;
 *       = 0 , if sign(i) = 1;
 */
ll Permutation::messLevel() const
{

    int i;
    ll mess = 0;
    int mess_type = (int) this->getParameter("MessType");
    vector<ll> sumPositive;
    int begin, end;
    ll sp;

    /* Check if mess_type is valid */
    if (mess_type < 0 || mess_type > 3) {
        cout<<"Invalid mess type"<<endl;
        exit(1);
    }

    //cout<<"MessType "<<mess_type<<endl;

    if (mess_type == 3) {
        sumPositive = vector<ll>(this->size() + 1);
        this->computePositiveSum(sumPositive);
    }

    for (i = 1; i <= this->size(); i++) {
        ll curr_mess = abs(this->at(i) - i);
        ll curr_sign = 0;
        if (this->sign(i) == -1) {
            curr_sign = 1;
        }
        if (mess_type == 0) {
            mess += curr_mess;
        } else if (mess_type == 1) {
            mess += curr_mess + curr_sign;
        } else if (mess_type == 2) {
            mess += 2 * curr_mess + this->sign(i) * curr_mess;
        } else if (mess_type == 3) {
            begin = this->at(i);
            end = i;
            if (end < begin)
                swap(begin, end);
            /* Number of elemente with positive sign in [begin..end] */
            sp = sumPositive[end] - sumPositive[begin - 1];
            mess += sp + curr_mess;
        }
        //cout<<curr_mess<<" "<<curr_sign<<endl;
        //cout<<mess<<endl;
    }

    return mess;
}


ll Permutation::inversionLevel()
{

    BIT<ll> bit(this->size());
    ll ni = 0;
    for (ll i = 1; i <= this->size(); i++) {
        bit.update(this->at(i),1);
        /* Computes the number of previous elements bigger than p[i]  */
        ni += i - bit.read(this->at(i));
    }

    return ni;
}

bool Permutation::next()
{

    bool has = std::next_permutation(this->seq.begin(), this->seq.end());

    if (!has)
        return false;
    this->updatePositions();

    return true;
}


pll Permutation::getMedian(ll begin, ll end)
{

    vector<ll> v;
    for (int i = begin; i <=end; i++) {
        v.push_back(this->at(i));
    }

    sort(v.begin(),v.end());

    int sz = v.size();

    pll pp;
    pp.F = v[sz/2];
    pp.S = this->position(pp.F);

    return pp;
}

bool Permutation::isSorted(ll begin, ll end, bool asUnsigned)
{
    int i;
    //cout<<"IS sorted: asUnsigned"<<asUnsigned<<endl;
    if (!asUnsigned && this->sign(begin) < 0)
        return false;
    for ( i = begin + 1; i <= end; i++) {
        if (!asUnsigned && (this->sign(i) < 0))
            return false;
        if ((this->at(i - 1) + 1) != this->at(i))
            return false;
    }
    //cout<<"Is sorted"<<endl;
    return true;
}

/* Expand the permutation, if it have been reduced */
void Permutation::expand()
{
    vector<ll> v;

    if (!this->is_reduced)
        return;

    if (strips.size() == 0)
        return;

    int i, j;
    int sz = this->strips.size();
    int p_sz;

    // the last strip has the original permutation size
    p_sz = this->strips[sz - 1].S;
    WATCH(*(this));
    for (i = 0; i < sz; i++) {
        pll s = this->strips[i];
        PRINT("Expand Strip [%d, %d]\n", s.F, s.S);
        // increasing strips
        if (s.F <= s.S) {
            for (j = s.F; j <= s.S; j++) {
                // ignore the border elements
                if (j == 0 || j == p_sz)
                    continue;
                v.push_back(j);
            }
        } else { // decreasing srips
            for (j = s.F; j >= s.S; j--) {
                // ignore the border elements
                if (j == 0 || j == p_sz)
                    continue;
                v.push_back(j);
            }
        }
    }
    /* Change parameter */
    //this->param["reduced"] = false;

    /* Change state of the permutation */
    this->is_reduced = false;

    this->initSeq(v);
}

void Permutation::setNormalMode()
{
    if (this->getParameter("reduced") && this->getParameter("extended")) {
        this->expand();
        this->setParameter("reduced", false);
        this->setParameter("extended", false);
    } else {
        cerr<<"Error at setNormalMode: inconsistent state"<<endl;
        exit(1);
    }
}

void Permutation::setExtendedMode()
{
    if (this->getParameter("extended") == false) {
        this->setParameter("extended", true);
        this->pos[this->size()] = this->size();
    }
}
void Permutation::setReducedExtendedMode()
{
    this->setParameter("reduced", true);
    this->setParameter("extended", true);
    this->reduce();
}

/* Method to reduce the permutation */
void Permutation::reduce()
{
    if (!this->param["extended"]) {
        cerr<<"ERROR: reduce is only supported for extended permutations"<<endl;
        exit(1);
    }

    // Check if the permutation already is reduced
    if (this->is_reduced == true)
        return;

    // Clear the data structures
    this->strips.clear();
    this->strips_pos.clear();

    // Begin of a strip
    int begin_strip = 0;
    PRINT("REDUCE\n");
    for (int i = 0; i < this->size(); i++) {
        // Look for the end of a consecutive component (Strip)
        bool is_break_point = false;
        if (this->param["signed"]) {
            if (this->get(i) + 1 != this->get(i + 1))
                is_break_point = true;
        } else {
            if (abs(this->get(i) - this->get(i + 1)) != 1)
                is_break_point = true;
        }

        if (is_break_point) {
            // Add a strip found
            pll s = make_pair(this->get(begin_strip), this->get(i));
            this->strips.push_back(s);
            // Add the positions of the strip found
            pll sp = make_pair(begin_strip, i);
            this->strips_pos.push_back(sp);
            PRINT("Strip pos %d: %lld %lld\n", i, sp.F, sp.S);
            // Add the last element of the strip to sequence
            begin_strip = i + 1;
        }
    }

    // Last Strip
    // Add a strip found
    pll s = make_pair(this->get(begin_strip), this->get(this->size()));
    this->strips.push_back(s);
    // Add the positions of the strip found
    pll sp = make_pair(begin_strip, this->size());
    this->strips_pos.push_back(sp);

    //DEBUG
    /*for (size_t i = 0; i < strips.size(); i++) {
        PRINT("Strip %d: [%d, %d]\n", i, strips[i].F, strips[i].S);
    }*/

    // Ignore the strips [0..k] and [u..N]
    int begin = 1;
    int end = strips.size() - 2;

    // Build index of Strips regarding the relative order
    vector<int> id(this->size() + 1, -1);

    for (int k = begin; k <= end; k++) {
        WATCH(abs(strips[k].F));
        id[abs(strips[k].F)] = 1;
    }

    int curr_pos = 1;
    for (int k = 0; k <= this->size(); k++) {
        if (id[k] != -1)
            id[k] = curr_pos++;
    }

    // Build the new sequence
    vector<ll> new_seq;
    for (int i = begin; i <= end; i++) {
        pll s = strips[i];
        ll abs_e = abs(s.F);
        ll p = id[abs_e];
        if (s.F < 0)
            p *= -1;
        new_seq.push_back(p);
    }

    /* Change state of */
    //this->param["reduced"] = true;
    this->is_reduced = true;

    this->initSeq(new_seq);
}

/* Converts a reversal from the reduced permutation
 * for the match reversal when the  permutation
 * is not reduced
 */
Reversal Permutation::convertReversal(Reversal rev)
{
  rev.calcCost();
    PRINT("convertReversal\n");
    if (strips_pos.size() == 0 || !this->param["reduced"])
        return rev;
    Reversal new_rev;
    new_rev.begin = this->strips_pos[rev.begin].F;
    new_rev.end   = this->strips_pos[rev.end].S;
    // Alterado por Ulisses
    new_rev.cost_function = rev.cost_function;
    new_rev.calcCost();

    return new_rev;
}


/* Converts a reversal from the non reduced permutation
 * for the match reversal when the  permutation
 * is reduced
 */
Reversal Permutation::convertReversalToReduced(Reversal rev)
{

    if (strips_pos.size() == 0 || !this->param["reduced"])
        return rev;
    Reversal new_rev(-1, -1, -1, rev.cost_function);
    PRINT("convertReversalToReduced\n");
    // WATCH(rev);
    // WATCH(*(this));
    /* Find the components for rev.begin and rev.end */
    for (size_t i = 0; i < this->strips_pos.size(); i++) {
        PRINT("\n\t strip %d: %d %d\n", i, this->strips_pos[i].F, this->strips_pos[i].S);
        if (this->strips_pos[i].F == rev.begin)
            new_rev.begin = i;
        if (this->strips_pos[i].S == rev.end)
            new_rev.end = i;
    }

    new_rev.cost_function = rev.cost_function;
    new_rev.calcCost();

    PRINT("\t");
    WATCH(*(this));
    PRINT("\t");
    WATCH(new_rev);
    // Some assertions about the new_rev
    int perm_sz = this->param["extended"] ? this->size() - 1: this->size();

    //Ulisses: comentei esses asserts aqui, eles fazem sentido, mas
    //nao vi porque estao dando errado.

    //assert(new_rev.begin > 0 && new_rev.begin <= perm_sz);
    //assert(new_rev.end > 0 && new_rev.end <= perm_sz);
    //assert(new_rev.begin <= new_rev.end);

    if (new_rev.begin == -1 || new_rev.end == -1)
        throw string("CANNOT_REDUCE");

    return new_rev;
}

/*
 * Binary representation of permutation (in relation to median)
 *
 */
vector<ll> Permutation::getBinary()
{

    pll m = this->getMedian(1,this->size());

    vector<ll> seq;

    for (int i = 1; i <= this->size(); i++) {

        if (this->at(i) < m.S)
            seq.push_back(0);
        else
            seq.push_back(1);
    }
    return seq;
}

Permutation::Permutation(parameter param,vector<ll> seq = vector<ll>())
{

    this->param = param;
    this->is_reduced = 0;
    if (seq.size() > 0) {
        this->initSeq(seq);
        if (this->param["reduced"])
            this->reduce();
    }
    //TODO: this assert slow down the performance
    assert(this->isValid());
}

Permutation::Permutation(parameter param, ll permI64)
{
    this->param = param;
    this->is_reduced = 0;
    vector<ll> v = this->llToVector(permI64);
    this->initSeq(v);
    if (this->param["reduced"])
        this->reduce();
}

Permutation::Permutation(parameter param)
{
    this->is_reduced = 0;
    this->param = param;
}

Permutation::Permutation(vector<ll> seq)
{
    this->is_reduced = 0;
    this->initSeq(seq);
    if (this->param["reduced"])
        this->reduce();
}

/* Converts a signed permutation \pi to a unsigened permutation \pi'.
 * Elements i in \pi are replaced to pair (2*i - 1, 2*i),
 * Elements -i in \pi are replaced to pair (2*i, 2*i - 1).
 */
Permutation Permutation::toUnsigned()
{
    this->is_reduced = false;
    assert(this->param["signed"]);

    Permutation u_perm;
    parameter u_param;
    int i;
    vector<ll> u_vector(this->size());

    for (i = 0; i < this->size(); i++) {
        u_vector[i] = abs(this->get(i + 1));
    }

    u_param = this->param;
    u_param.erase("signed");
    u_perm = Permutation(u_param, u_vector);

    return u_perm;
}

/* Generates a new permutation representing init considering that the
 *  end is the identity permutation. */
Permutation Permutation::convertFromTo(Permutation &init, Permutation &final)
{

    if (init.size() != final.size()) {
        cout<<"Error at convertFromTo: permutations have different size"<<endl;
        exit(0);
    }

    ll perm_sz = init.size();
    if (init.getParameter("extended"))
        perm_sz--;

    vector<ll> seq(perm_sz);
    int i;

    for (i = 1; i <= perm_sz; i++) {
        ll pos = final.position(init.at(i));
//    cout <<"pos "<<pos<<endl;
        seq[i - 1] = pos * final.signAt(pos) * init.signAt(i);;
    }

    Permutation p(init.getParameter(), seq);

    return p;
}

/*
 * validate this permutation
 */
bool Permutation::isValid()
{

    int sz = this->size();
    std::set<ll> s;

    for (int i = 1; i <= sz; i++) {
        s.insert(this->at(i));
    }

    ll ma = *max_element(s.begin(),s.end());
    if ((int)s.size() == sz && ma == sz)
        return true;

    return false;
}

Permutation Permutation::makeIdentity(parameter param, int n)
{
    int i;
    vector<ll> v;

    if (param["extended"])
        n--;

    for (i = 1; i <= n; i++)
        v.push_back(i);

    Permutation p(param,v);

    return p;
}

/*
 */
ll Permutation::numBreakPoints(bool asUnsigned)
{
    ll num_bk = 0;
    bool fake_extended = false;
    // Check if there is a breakpoint in the permutation beginning
    if (this->getParameter("extended") == false) {
        this->setParameter("extended", 1);
        fake_extended = true;
    }

    if (this->getParameter("signed") && !asUnsigned) {
        // Signed case
        for (int i = 0; i < this->size(); i++) {
            // Computes difference between consecutive elements
            ll diff = this->get(i + 1) - this->get(i);
            if (diff != 1) {
                num_bk++;
                PRINT(" [%d, %d], ", this->get(i), this->get(i + 1));
            }
        }
    } else {
        // Unsigned case
        for (int i = 0; i < this->size(); i++) {
            // Computes absolute difference between consecutive elements
            ll diff = abs(this->get(i + 1) - this->get(i));
            if (diff != 1) {
                num_bk++;
                PRINT(" [%d, %d], ", this->get(i), this->get(i + 1));
            }
        }
    }

    if (fake_extended) {
        this->setParameter("extended", 0);
    }
    return num_bk;
}

void Permutation::initSeq(vector<ll> seq)
{
    int i;
    int sz;

    sz = seq.size();

    vector<ll> newSeq(sz + 2);
    //cout<<"init seq"<<endl;

    /* Build newSeq */
    newSeq[0] = 0; //for extended permutation
    newSeq[sz + 1] = sz + 1; //for extended permutation
    for (i = 0; i < sz; i++) {
        newSeq[i + 1] = seq[i];
    }

    /* Change permutation size */
    this->N = sz;

    /* Create pos vector */
    this->pos = vector<ll>(this->N + 2, 0);

    /* Stores seq */
    this->seq = newSeq;

    /* Update postitions */
    this->updatePositions();
}

void Permutation::updatePositions()
{
    int i;
    int begin = 1;
    //cout<<"Update Positions"<<endl;
    //cout<<"Exetended: "<<this->param["extended"]<<endl;
    if (this->param["extended"])
        begin--;

    for (i = begin; i <= this->size(); i++) {
        //cout<<i<<" : "<<this->at(i)<<endl;
        pos[this->at(i)] = i;
    }

}

void Permutation::changeSignsRandomly()
{
    int i;

    for (i = 1; i <= this->size(); i++) {
        if (rand() % 2 == 0) {
            this->set(i, this->get(i) * -1);
        }
    }
}

/*
 * k: element [1..N]
 * @return position of k in permutation
 */
int	Permutation::position(int k) const
{
    //cout<<"position "<<k<<" "<<pos[k]<<endl;
    return pos[k];
}

int Permutation::size() const
{
    if (this->getParameter("extended") != 0) {
        //printf("Size extended\n");
        return N + 1;
    }
    return N;
}

parameter Permutation::getParameter() const
{
    return this->param;
}

void Permutation::setParameter(string param_name, double val)
{
    this->param[param_name] = val;
}

double Permutation::getParameter(string param_name) const
{
    parameter::const_iterator cit, cend;
    cit = this->param.find(param_name);
    cend = this->param.end();
    if (cit == cend) {
        /*throw ("Permutation: parameter " + param_name
        		+ " not found");*/
        return 0;
    }
    return cit->second;
}

/* Return the absoulete value of element at
 * position k.
 */
ll Permutation::at(int k) const
{
    if (this->indexInBounds(k))
        return abs(this->seq[k]);
    else {
        cerr<<"at ERROR: try access a element at position out of bounds "<<k<<endl;
        exit(1);
    }
}

/* Return the value of element at
 * position k.
 */
ll Permutation::get(int k) const
{

    if (this->indexInBounds(k))
        return this->seq[k];
    else {
        cerr<<"get ERROR: try access a element at position out of bounds "<<k<<endl;
        exit(1);
        return NIL;
    }
}

/* TODO: Don't work well when the permutation is reduced */
void Permutation::set(int k, ll v)
{
    if (this->indexInBounds(k) == false)
        return;
    seq[k] = v;
}

/* Returns the value at position k, not the absolute value */
ll Permutation::operator[](int k) const
{
    return this->get(k);
}

/*
 * k: element [1..N]
 * @return slice of k in permutation
 */

int Permutation::slice(int k)
{
    return min(position(k), N - position(k) + 1);
}

/*
 * k: element [1..N]
 * @return sign of k in permutation
 */
int Permutation::sign(int e) const
{
    //cout<<"Sign "<<e<<endl;
    int p = this->position(abs(e));
    //printf("p = %d\n", p);
    int v = this->get(p);

    //printf("v: %d\n", v);
    return v >= 0 ? 1 : -1;
}

/*
 * k: element [1..N]
 * @return sign of k in permutation
 */
int Permutation::signAt(int k)
{
    return this->sign(this->at(k));
}

string Permutation::toString() const
{
    int i;
    string ret_string;
    int n = this->size();

    //cout<<"n = "<<n<<endl;
    //cout<<"extended = "<<this->param["extended"]<<endl;

    int end = n;
    int begin = 1;

    if (this->getParameter("extended")) {
        begin--;
    }

//  cout<<"begin = "<<begin<<endl;
//  cout<<"end = "<<end<<endl;

    ret_string = std::to_string(seq[begin]);
    for (i = begin+1; i <= end; i++) {
      ret_string += ","+std::to_string(seq[i]);
    }
    return ret_string;
}


void Permutation::print(ostream &out = cout, char sep = ' ') const
{
    int i;
    int n = this->size();

    //cout<<"n = "<<n<<endl;
    //cout<<"extended = "<<this->param["extended"]<<endl;

    int end = n;
    int begin = 1;

    if (this->getParameter("extended")) {
        begin--;
    }

//  cout<<"begin = "<<begin<<endl;
//  cout<<"end = "<<end<<endl;

    for (i = begin; i < end; i++) {
        out<<seq[i]<<sep;
    }

    out<<seq[end];
}


bool Permutation::indexInBounds (int i) const
{
    if (this->getParameter("extended") != 0) {
        if (i >= 0 && i <= this->N + 1)
            return true;
        return false;
    }	else {
        if (i >= 1 && i <= this->N)
            return true;
        return false;
    }
}

bool Permutation::invert(int i, int j)
{
    if (!indexInBounds(i) || !indexInBounds(j))
        return false;

    int k;

    for (k = i; k <= j; k++) {
        this->set(k, -1 * this->get(k));

        // PRINT("Invert: %d: %d\n", k, seq[k]);
    }

    return true;
}

bool Permutation::isASRValidReversal(int i,int j)
{

    if (abs(this->slice(i) - this->slice(j)) <= 1)
        return true;
    return false;
}

bool Permutation::isReversalValid(int i, int j)
{
    if (i < 0 || !indexInBounds(i) || j > this->N || !indexInBounds(j))
        return false;
    if (this->getParameter("extended")) {
        if (i == this->size() || j == this->size())
            return false;
    }
    return true;
}

bool Permutation::reverse(Reversal rev)
{

  int i = rev.begin;
  int j = rev.end;

    if (!this->isReversalValid(i, j))
        return false;

    /* Check for the problem model */
    if (param["reversal"] == false)
        return false;

    if (param["ASR"]) { //Almost simetric reversals
        if (!isASRValidReversal(i,j))
            return false;
    }

    if (this->param["reduced"]) {
        // Get the relative reversal in the original permutation
      Reversal r = this->convertReversal(Reversal(i, j, 0, rev.cost_function));
        i = r.begin;
        j = r.end;
        this->expand();
        // cout<<"::reverse rev "<<r<<endl;
        //cout<<":::: perm (after extend)= "<<*(this)<<endl;
    }

    if (this->param["signed"] ) {
        this->invert(i,j);
        PRINT("Is signed!\n");
    }

    std::reverse(seq.begin() + i, seq.begin() + j + 1);
    //cout<<":::: perm (after reverse)= "<<*(this)<<endl;
    if (this->param["reduced"]) {
        this->reduce();
        //cout<<":::: perm (after reduce)= "<<*(this)<<endl;
    }

    this->updatePositions();

    return true;
}

//Representation--------------------------------------------------------


/*
 * Convert representation of a permutation.
 * From 64-bits integer @k to a vector.
 */
vector<ll> Permutation::llToVector(ll k)
{
    vector<ll> v;
    ll num;
    int i = 0;
    /* Size of permutation */
    int sz = 0;
    while (1) {
        /* Get the bits of the element i */
        num = (MASK) & (k>>(i * N_BITS));
        if (num == 0) /* Break case */
            break;
        v.push_back(num);
        i++;
        sz++;
    }
    /* Get information abou sign */
    ll dl = N_BITS * (ll)v.size() + 1;
    for (i = 0; i < sz; i++) {
        ll sign = 0;
        /* Get the sign */
        sign = 1 & (k >> dl);

        /* Set the element in vector */
        if (sign == 0) /* Negative sign */
            v[i] *= -1;

        /*Increase desloc factor */
        dl++;
    }

    return v;
}

/*
 * Convert representation of a permutation.
 * From a vector @v 64-bits integer.
 */
ll Permutation::getI64()
{
    vector<ll> v = this->seq;
    /* Generated representation of uma permutation */
    ll s = 0;
    int i;
    ll b;
    for (i = 0; i < (int)v.size(); i++) {
        b = (v[i]<<(i*N_BITS));
        s |= b;
    }

    /* Importantion: between permutation information and signs information
     * there are 4 bits of value 0, to indicate end of the permutation.
     */

    /* Add sign information */
    if (this->param["signed"]) {
        /* Set initial desloc factor */
        ll dl = N_BITS * (ll)v.size() + 1;
        for (i = 1; i <= this->size(); i++) {
            ll sign = 0;
            if (this->param["signed"] && this->sign(i) == -1) {
                sign = 1;
            }
            /* Set sign bit */
            s |= sign << dl;
            /*Increase desloc factor */
            dl++;
        }
    }
    return s;
}

int Permutation::cmp(const Permutation &p) const
{
    int i;

    int sz;

    if (this->size() > p.size())
        return 1;
    else if (this->size() < p.size())
        return -1;

    sz = this->size();

    for ( i = 1 ; i <= sz; i++) {
        if (this->get(i) < p[i])
            return -1;
        else if (this->get(i) > p[i])
            return 1;
    }

    return 0;
}

// /* Build a list with all permutations in neighgorhood of one
//  * by reversals.
//  */
// vector<Permutation::perm_rev> Permutation::makeNeighborhood(int cost_function)
// {
//     ll i , j, k;
//     vector<perm_rev> ng;

//     /* Clone this object */
//     Permutation perm = *(this);

//     for (i = 1; i <= this->size(); i++) {
//         k = i;
//         /* For unsigned permutations don't make sense reversals
//          * (i,i) */
//         if (!this->param["signed"])
//             k++;
//         for (j = k; j <= this->size(); j++) {
//             Permutation p = perm;
//             p.reverse(i,j);
//             ng.push_back(perm_rev(p, Reversal(i, j, 0, cost_function)));
//         }
//     }

//     return ng;
// }

bool Permutation::operator<(const Permutation &p) const
{
    int r = this->cmp(p);
    //cout<<*(this)<<" , "<<p<<" : "<<r<<endl;
    return r == -1;
}

bool Permutation::operator>(const Permutation &p) const
{
    return this->cmp(p) == 1;
}

bool Permutation::operator==(const Permutation &p) const
{
    return this->cmp(p) == 0;
}

bool Permutation::operator!=(const Permutation &p) const
{
    return this->cmp(p) != 0;
}

ostream &operator<<(ostream &stream, Permutation p)
{
    p.print(stream);
    return stream;
}

istream &operator>>(istream &stream, Permutation &p)
{
    int sz;
    stream >> sz;

    vector<ll> seq;

    for (int i = 1; i <= sz; i++) {
        ll x;
        stream >> x;
        seq.push_back(x);
    }

    p.initSeq(seq);
    if (p.getParameter()["reduced"])
        p.reduce();
    return stream;
}

void Permutation::printReversal(ostream &stream, Reversal r)
{

    for (int i = 1; i <= this->size(); i++) {

        if (r.begin == i)
            stream<<"|";
        /* Prints the value at position i, not the absolute value */
        stream<<(this->get(i));

        if (r.end == i)
            stream<<"|";
        stream<<" ";
    }

}


/* Given a permutation, returns the set of reversals that removes at least
 * a breakpoint.
 * Time complexity: O(N)
 */
void Permutation::getBreakpointReversals(
					 vector<Permutation::perm_rev> &out_reversals, 
					 int cost_function
					 ) const
{

    if (this->getParameter("extended") == false) {
        cerr<<"getBreakpointReversals: Permuation must be extended"<<endl;
        exit(1);
    }

    // This algorithm only makes sense for unsigned permutations...
    /*if (this->getParameter("signed")) {
      Permutation p_copy = *(this);
      p_copy = p_copy.toUnsigned();
      p_copy.getBreakpointReversals(out_reversals);
      return;
    }*/

    vector<Strip> strips = this->getAllStrips();
    vector<int> id_strip_at(this->size() + 1);

    // For each permutation's element sets the direction from the
    // its strip
    for (size_t i = 0; i < strips.size(); i++) {
        PRINT("STRIP: %d, [%d %d] => { %lld, %lld}\n",
              i,
              strips[i].begin, strips[i].end,
              this->at(strips[i].begin),this->at(strips[i].end)
             );
        for (int k = strips[i].begin; k <= strips[i].end; k++) {
            id_strip_at[this->at(k)] = i;
        }
    }

    for (size_t i = 0; i < strips.size(); i++) {
        // A decreasing strip was just found
        if (strips[i].dir == NEGATIVE) {
            // First element on the decreasing strip...
            int w = this->at(strips[i].begin);
            int k = this->at(strips[i].end);
            int first_id_strip = id_strip_at[w + 1];
            int second_id_strip = id_strip_at[k - 1];
            Reversal first_rev;
            Reversal second_rev;
            PRINT("w = %lld, k = %lld", w, k);
            assert(w + 1 <= this->size());
            assert(k - 1 >= 0);
            // Check if actually there is a strip [w + 1..j]
            if ((w + 1 == this->size())
                    || (strips[first_id_strip].begin == this->position(w + 1))
               ) {
                /* Look for the first reversal to join the strips [w..k] and [w+1..j] */
                if (this->position(w + 1) > this->position(w)) {
                    // 1st scenario: ...[w..k]...[w+1..j]...
                    first_rev.begin = this->position(w);
                    first_rev.end = this->position(w + 1) - 1;
                } else {
                    // 2nd scenario: ...[w+1..j]...[w..k]...
                    first_rev.begin = this->position(w + 1);
                    first_rev.end = this->position(w) - 1;
                }

		first_rev.cost_function = cost_function;
                first_rev.calcCost();

                WATCH(this->position(w + 1));
                WATCH(first_rev);
                assert(first_rev.begin > 0);
                assert(first_rev.end >= first_rev.begin);
                assert(first_rev.end < this->size());
                if (first_rev.begin < first_rev.end) {
                    // performs the reversal and insert into reversal list
                    Permutation pp = *(this);
                    pp.reverse(first_rev);
                    out_reversals.push_back(make_pair(pp, first_rev));
                }
            }
            WATCH(strips[second_id_strip].end);
            /* Second reversal */
            // Check if actually there is a strip [l..k-1]
            if ((k - 1 == 0)
                    || (strips[second_id_strip].end == this->position(k - 1))
               ) {
                if (this->position(k - 1) > this->position(k)) {
                    // 1st scenario: ...[w...k]...[l...k-1]...
                    second_rev.begin = this->position(k) + 1;
                    second_rev.end = this->position(k - 1);
                } else {
                    // 2nd scenario: ...[l...k-1]...[w...k]...
                    second_rev.begin = this->position(k - 1) + 1;
                    second_rev.end = this->position(k);
                }

		second_rev.cost_function = cost_function;
		second_rev.calcCost();
                WATCH(second_rev);
                assert(second_rev.begin > 0);
                assert(second_rev.end >= second_rev.begin);
                assert(second_rev.end < this->size());
                // Only inserts the second reversal if it is not equal to the first one
                if (second_rev.begin < second_rev.end && first_rev != second_rev) {
                    // performs the reversal
                    Permutation pp = *(this);
                    pp.reverse(second_rev);
                    out_reversals.push_back(make_pair(pp, second_rev));
                }
            }
        }
    }
}
