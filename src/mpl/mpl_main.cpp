#include <iostream>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/transform.hpp>

using namespace std;
using namespace boost;

struct A {
	static void evaluate() { cout << "Hello from A" << endl; }
};

struct B {
	static void evaluate() { cout << "Hello from B" << endl; }
};

struct C {
	static void evaluate() { cout << "Hello from C" << endl; }
};

template <class T>
struct active_trait {
	static const bool value = false;
};

template <>
struct active_trait<A> {
	static const bool value = true;
}

template <>
struct active_trait<C> {
	static const bool value = true;
}

struct evaluator {
	template <class T>
	void operator()(T) const {
		T::evaluate();
	}
};

int main( int argc, char * argv[] )
{
	typedef mpl::vector<A, B, C> all_items;
	typedef mpl::transform<all_items, active_traits> items;
	mpl::for_each<items>(evaluator());

}
