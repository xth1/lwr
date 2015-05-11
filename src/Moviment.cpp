#include "Moviment.h"

// Alterado por Ulisses toda vez que aparecer cost_function

double Moviment::getCost()
{
    return NIL;
}


// Eu mantive o cost porque algumas funcoes antigas queriam chamar a
// opcao 2. Jah Jah eu preciso mudar isso.
Reversal::Reversal(ll begin,ll end,ll cost, int cost_function)
{
  assert(cost_function >= 1);
  assert(cost_function <= 3);  

    this->begin = begin;
    this->end = end;
    this->cost_function = cost_function;
    this->calcCost();
}

// Reversal::Reversal(ll begin,ll end, int cost_function)
// {
//   assert(cost_function >= 1);
//   assert(cost_function <= 3);  
//     this->begin = begin;
//     this->end = end;
//     this->cost_function = cost_function;
//     this->calcCost();
// }

double Reversal::getCost()
{
    return (double)this->cost;
}

ostream &operator<<(ostream &stream, Reversal r)
{
    stream<<"Reversal: "<<"["<<r.begin<<" "<<r.end<<" "<<r.cost<<"]";
    return stream;
}

istream &operator>>(istream &stream, Reversal &r)
{
    stream>>r.begin>>r.end>>r.cost;
    return stream;
}

bool Reversal::operator==(Reversal &rev) const
{
    if (rev.begin == this->begin && rev.end == this->end) {
        assert(rev.cost == this->cost);
        return true;
    }
    return false;
}

bool Reversal::operator!=(Reversal &rev) const
{
    return !(*(this) == rev);
}

bool Reversal::operator<(Reversal &rev) const
{
    if (this->begin < rev.begin)
        return true;
    else if (this->begin == rev.begin)
        return this->end < rev.end;
    return false;
}

void Reversal::calcCost()
{
  assert(this->cost_function >= 1);
  assert(this->cost_function <= 3);
  if (this->cost_function == 1) {
    this->cost = this->end - this->begin + 1;
  } else if (this->cost_function == 2) {
    this->cost = pow(this->end - this->begin + 1, 2);
  } else if (this->cost_function == 3) {
    this->cost = log2(this->end - this->begin + 1) + 1;
  }
}
