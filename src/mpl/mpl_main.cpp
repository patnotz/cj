#include <iostream>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>

using namespace std;
using namespace boost;

struct A {
  A() :
    value(0) {
  }
  A(const int val) :
    value(val) {
  }
  int value;
  void evaluate() {
    cout << "Hello from A with value = " << value << endl;
  }
  template <class S> void setup(S & s) {
    cout << "setup A" << endl;
  }
};

struct B {
  B() :
    value(0), a(0) {
  }
  B(const int val) :
    value(val), a(0) {
  }
  int value;
  int a;
  void evaluate() {
    cout << "Hello from B with value = " << value*a << endl;
  }
  template <class S> void setup(S & s) {
    a = (*fusion::find<A>(s)).value;
    cout << "setup B, a = " << a << endl;
  }
};

struct C {
  C() :
    value(0) {
  }
  C(const int val) :
    value(val) {
  }
  int value;
  void evaluate() {
    cout << "Hello from C with value = " << value << endl;
  }
  template <class S> void setup(S & s) {
    cout << "setup C" << endl;
  }
};

struct D {
  D() :
    value(0) {
  }
  D(const int val) :
    value(val) {
  }
  int value;
  void evaluate() {
    cout << "Hello from D with value = " << value << endl;
  }
  template <class S> void setup(S & s) {
    cout << "setup D" << endl;
  }
};

template <class S>
struct manager {
  template <class F>
  void apply(F & f) {
    fusion::for_each(seq,f);
  }
  S seq;
};

template <class S>
struct evaluator {
  evaluator(S & s) : seq(s) {}
  template<class T>
  void operator()(T & t) const {
    t.setup(seq);
    t.evaluate();
  }
  S & seq;
};

int main(int argc, char * argv[]) {
  typedef fusion::list<A, B, C> Seq;
  manager<Seq> m;
  if(fusion::has_key<A>(m.seq))
    *fusion::find<A>(m.seq) = A(9);
  if(fusion::has_key<B>(m.seq))
    *fusion::find<B>(m.seq) = B(2);
  if(fusion::has_key<C>(m.seq))
    *fusion::find<C>(m.seq) = C(7);
  evaluator<Seq> e(m.seq);
  m.apply(e);

  //typedef typename fusion::result_of::find<Seq, B>::type Pos;
  //typedef typename fusion::result_of::insert<Seq,Pos,D>::type Seq2;
  //Seq2 s2;
  //manager<Seq2> m2;
}
