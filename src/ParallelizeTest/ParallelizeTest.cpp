
#include <list>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

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
        //for(int i = 0; i < 10000; i++)
        val /= n;
    }

    void multiply(int n) {
        for (int i = 0; i < 10000; i++)
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

/* Global */

int NOMBRE_ELEM = 10000000;

vector<Element> elements;

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

void simple_vector(vector<Element> &v, VisitorTest &visitor) {
    for (int i = 0; i < v.size(); i++) {
        visitor.visit(v[i]);
    }
}

void openmp_vector(vector<Element> &v, VisitorTest &visitor) {
#pragma omp parallel for
    for (int i = 0; i < v.size(); i++) {
        visitor.visit(v[i]);
        // cout << v[i].value() << " ";
    }
}

// LIST

void simple_list(list<Element> &l, VisitorTest &visitor) {
    auto iter = l.begin();
    for (; iter != l.end(); iter++) {
        visitor.visit(*iter);
    }
}

void openmp_list_copy(list<Element> &l, VisitorTest &visitor) {
    std::vector<Element *> elements;
    for (auto iter = l.begin(); iter != l.end(); ++iter)
        elements.push_back(&(*iter));

    #pragma omp parallel for
    for (int i = 0; i < elements.size(); i++) {
        visitor.visit(*(elements[i]));
        // cout << l[i].value() << " ";
    }
}

void openmp_list_copy2(list<Element> &l, VisitorTest &visitor) {
    std::vector<Element *> items;

    items.reserve(l.size());
    //put the pointers in the vector
    transform(l.begin(), l.end(), back_inserter(items),
              [](Element &n) { return &n; }
    );

    #pragma omp parallel for
    for (int i = 0; i < items.size(); i++) {
        visitor.visit(*(items[i]));
    }
}


void openmp_list_nowait(list<Element> &l, VisitorTest &visitor) {
    list<Element>::iterator iter;
    #pragma omp parallel private(iter)
    {
        //Every thread has its own copy of the for loop
        //And each one has a copy of iter
        for (iter = l.begin(); iter != l.end(); iter++) {
        #pragma omp single nowait
            {
                //Only one per iteration (single), and the other don't wait (nowait)
                iter->accept(visitor);
            }
        }
    }
    //Very long solution, why ?
    // Maybe because simulated tasks are too short
}

void openmp_list_task(list<Element> &l, VisitorTest &visitor) {
    #pragma omp parallel
    {
    #pragma omp single
        {
            for (auto iter = l.begin(); iter != l.end(); ++iter) {
    #pragma omp task firstprivate(iter)
                iter->accept(visitor);
            }
        }
    }
    //Only with openmp 3.0, but very VERY long, why?
    //Maybe too much tasks...
}

/*------------- Tests finaux -------------*/

void vector_test() {
    vector<Element> v;
    ConcreteVisitorTest visitor(2);
    float time;

    // Simple thread test
    v.assign(elements.begin(), elements.end());
    cout << "Simple vector..   ";
    time = count_time(simple_vector, v, visitor);
    cout << " Done in " << time << " ms." << endl;

    // OpenMP test
    v.clear();
    v.assign(elements.begin(), elements.end());
    cout << "OpenMP vector..   ";
    time = count_time(openmp_vector, v, visitor);
    cout << " Done in " << time << " ms." << endl;
}

void list_test() {
    list<Element> l;
    ConcreteVisitorTest visitor(2);
    float time;

    // Simple thread test
    l.assign(elements.begin(), elements.end());
    cout << "Simple list..   ";
    time = count_time(simple_list, l, visitor);
    cout << " Done in " << time << " ms." << endl;

    //OpenMP list copy
    l.clear();
    l.assign(elements.begin(), elements.end());
    cout << "OpenMP list copy..   ";
    time = count_time(openmp_list_copy, l, visitor);
    cout << " Done in " << time << " ms." << endl;

    //OpenMP list copy2
    l.clear();
    l.assign(elements.begin(), elements.end());
    cout << "OpenMP list copy2..   ";
    time = count_time(openmp_list_copy2, l, visitor);
    cout << " Done in " << time << " ms." << endl;


    //VEry long solutions......
    /*
    //OpenMP list nowait
    init_list(l);
    cout<<"OpenMP list nowait..   ";
    time = count_time(openmp_list_nowait, l, visitor);
    cout<<" Done in "<<time<<" ms."<<endl;

    //OpenMP list task
    init_list(l);
    cout<<"OpenMP list task..   ";
    time = count_time(openmp_list_task, l, visitor);
    cout<<" Done in "<<time<<" ms."<<endl;

    */
}

/*------------- Main -------------*/

int main(int argc, char const *argv[]) {

    if (argc > 1) {
        NOMBRE_ELEM = atoi(argv[1]);
    }

    cout << "Launching tests with " << NOMBRE_ELEM << " elements.\n";

    init_vector(elements);

    vector_test();

    list_test();

    return 0;
}