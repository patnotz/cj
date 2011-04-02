#include <iostream>
#include <string>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>

using namespace std;
using namespace boost;

struct Base {
  Base(const int value_arg, const string name_arg) :
    value(value_arg), name(name_arg) {
  }
  int value;
  string name;
  void evaluate() {
    cout << "Hello from " << name << " with value = " << value << endl;
  }
};

struct A : public Base{
  A() : Base(0, "A") {}
  A(const int val) : Base(val, "A") {}
  template <class S>
  void setup(S & s) {
    cout << "setup " << name << endl;
  }
};

struct B : public Base {
  B() : Base(0, "B"), b(0) {}
  B(const int val) : Base(val, "B"), b(val) {}
  int b;
  template <class S>
  void setup(S & s) {
    int a = (*fusion::find<A>(s)).value;
    value = a * b;
    cout << "setup B, a = " << a << endl;
  }
};

struct C : public Base {
  C() : Base(0,"C") {}
  C(const int val) : Base(val, "C") {}
  template <class S>
  void setup(S & s) {
    cout << "setup " << name << endl;
  }
};

struct D : public Base {
  D() : Base(0,"D") {}
  D(const int val) : Base(val, "D") {}
  template <class S>
  void setup(S & s) {
    cout << "setup " << name << endl;
  }
};

struct X;
struct Y;
typedef fusion::map<
    fusion::pair<X, fusion::list<C> >,
    fusion::pair<Y, fusion::list<D> > > opt;

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

  typedef fusion::result_of::value_at_key<opt,Y>::type ys;
  typedef fusion::result_of::begin<ys>::type yi;
  typedef fusion::result_of::value_of<yi>::type X_type;
  typedef fusion::cons<X_type, Seq> Seq2;
  manager<Seq2> m2;
  m2.set( m.get<A>() );
  m2.set( m.get<B>() );
  m2.set( m.get<C>() );
  m2.set( X_type(1) );
  m2.run();
}
