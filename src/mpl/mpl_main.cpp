#include <iostream>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/container.hpp>

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
};

struct B {
  B() :
    value(0) {
  }
  B(const int val) :
    value(val) {
  }
  int value;
  void evaluate() {
    cout << "Hello from B with value = " << value << endl;
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
};

template<class T>
struct active_trait {
  static const bool value = false;
};

template<>
struct active_trait<A> {
  static const bool value = true;
};

template<>
struct active_trait<C> {
  static const bool value = true;
};

struct evaluator {
  template<class T>
  void operator()(T & t) const {
    t.evaluate();
  }
};

int main(int argc, char * argv[]) {
  typedef fusion::list<A, B, C> Seq;
  Seq seq;
  fusion::at_c<0>(seq) = A(9);
  fusion::at_c<1>(seq) = B(2);
  fusion::at_c<2>(seq) = C(7);
  fusion::for_each(seq, evaluator());
}
