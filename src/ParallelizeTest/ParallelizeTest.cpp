
#include <list>
#include <iostream>
#include <vector>
#include <chrono>

#include "TimeDecorator.hpp"

#if COMPILER_MSVC
#  define DISABLE_OPTIMISATIONS() __pragma( optimize( "", off ) )
#  define ENABLE_OPTIMISATIONS() __pragma( optimize( "", on ) )
#elif COMPILER_GCC
#  define DISABLE_OPTIMISATIONS() \
        _Pragma( "GCC push_options" ) \
        _Pragma( "GCC optimize (\"O0\")" )
#  define ENABLE_OPTIMISATIONS() _Pragma( "GCC pop_options" )
#elif COMPILER_CLANG
#  define DISABLE_OPTIMISATIONS() _Pragma( "clang optimize off" )
#  define ENABLE_OPTIMISATIONS() _Pragma( "clang optimize on" )
#endif

using namespace std;

/* Global */

int NOMBRE_ELEM = 10000000;

/*------------- Class -------------*/

class Element;

class VisitorTest {
public:
    virtual void visit(Element &) = 0;
};

class Element {
protected:
    int val;
    bool even;

public:

    Element() : val(0), even(true) {};

    Element(int val) : val(val), even(val % 2 == 0) {};

    virtual void accept(VisitorTest &visitor) {
        visitor.visit(*this);
    }

    void divide(int n) {
        val /= n;
    }

    void multiply(int n) {
        val *= n;
    }

    const int value() const {
        return val;
    }

    const bool isEven() const {
        return even;
    }
};

class ConcreteVisitorTest : public VisitorTest {
private:
    int coef;

public:
    ConcreteVisitorTest(int coef) : coef(coef) {};

    virtual void visit(Element &element) override {
        if (element.isEven())
            element.divide(coef);
        else
            element.multiply(coef);
    };
};

/*------------- Init -------------*/

void init_vector(vector<Element> &vector) {
    vector.clear();

    cout << "Init vector.. " << endl;

    for (int i = 0; i < NOMBRE_ELEM; i++) {
        vector.push_back(Element(i));
        // cout << vector[i].value() << " ";
    }
    // cout << endl;
}

void init_list(list<Element> &list) {
    list.clear();

    cout << "Init list.. " << endl;

    for (int i = 0; i < NOMBRE_ELEM; i++) {
        list.push_back(Element(i));
        // cout << list[i].value() << " ";
    }
    // cout << endl;
}


/*------------- Fonctions de test -------------*/

// VECTOR

void simple_vector(vector<Element> &v, VisitorTest& visitor) {
    for (int i = 0; i < v.size(); i++) {
        visitor.visit(v[i]);
    }
}

void openmp_vector(vector<Element> &v, VisitorTest& visitor) {
    #pragma omp parallel for
    for (int i = 0; i < v.size(); i++) {
        visitor.visit(v[i]);
        // cout << v[i].value() << " ";
    }
}

// LIST

void simple_list(list<Element> &l, VisitorTest& visitor) {
    auto iter = l.begin();
    for (; iter != l.end(); iter++) {
        visitor.visit(*iter);
    }
}

void openmp_list_copy(list<Element> &l, VisitorTest& visitor) {
    std::vector<Element*> elements;
    auto iter = l.begin();
    for (iter = l.begin(); iter != l.end(); ++iter)
        elements.push_back(&(*iter));

    #pragma omp parallel for
    for (int i = 0; i < elements.size(); i++) {
        visitor.visit(*(elements[i]));
        // cout << l[i].value() << " ";
    }
}

/*------------- Tests finaux -------------*/

void vector_test() {
    vector<Element> v;
    ConcreteVisitorTest visitor(2);
    float time;

    // Simple thread test
    init_vector(v);
    cout<<"Simple vector..   ";
    time = count_time(simple_vector, v, visitor);
    cout<<" Done in "<<time<<" ms."<<endl;

    // OpenMP test
    init_vector(v);
    cout<<"OpenMP vector..   ";
    time = count_time(openmp_vector, v, visitor);
    cout<<" Done in "<<time<<" ms."<<endl;
}

void list_test() {
    list<Element> l;
    ConcreteVisitorTest visitor(2);
    float time;

    // Simple thread test
    init_list(l);
    cout<<"Simple list..   ";
    time = count_time(simple_list, l, visitor);
    cout<<" Done in "<<time<<" ms."<<endl;

    //OpenMP list copy
    init_list(l);
    cout<<"OpenMP list copy..   ";
    time = count_time(openmp_list_copy, l, visitor);
    cout<<" Done in "<<time<<" ms."<<endl;
}

/*------------- Main -------------*/

int main(int argc, char const *argv[]) {

    if (argc > 1) {
        NOMBRE_ELEM = atoi(argv[1]);
    }

    cout << "Launching tests with " << NOMBRE_ELEM << " elements.\n";
    
    vector_test();
    
    list_test();
    
    return 0;
}