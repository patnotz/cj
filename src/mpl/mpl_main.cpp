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
struct DoSetup {
  DoSetup(S & s) : seq(s) {}
  template<class T>
  void operator()(T & t) const {
    t.setup(seq);
  }
  S & seq;
};

struct DoEvaluate {
  template <class T>
  void operator()(T & t) const {
    t.evaluate();
  }
};

template <class S>
struct manager {
  manager() : seq(S()), do_setup(seq) {}
  void run() {
    fusion::for_each(seq,do_setup);
    fusion::for_each(seq,do_evaluate);
  }
  template <class T>
  void set(T t) {
    *fusion::find<T>(seq) = t;
  }
  template <class T>
  T & get() {
    return *fusion::find<T>(seq);
  }

  S seq;
  DoSetup<S> do_setup;
  DoEvaluate do_evaluate;
};

int main(int argc, char * argv[]) {
  typedef fusion::cons<A, fusion::cons<B, fusion::cons<C> > > Seq;
  manager<Seq> m;
  m.set(A(9));
  m.set(B(2));
  m.set(C(7));
  m.run();

  typedef fusion::cons<D, Seq> Seq2;
  manager<Seq2> m2;
  m2.set( m.get<A>() );
  m2.set( m.get<B>() );
  m2.set( m.get<C>() );
  m2.set( D(1) );
  m2.run();
}
