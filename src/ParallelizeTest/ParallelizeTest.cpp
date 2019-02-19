
#include <list>
#include <iostream>
#include <vector>
#include <chrono>

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

void init_vector(vector<Element> &vector) {
    vector.clear();

    cout << "Init vector.. " << endl;

    for (int i = 0; i < 10000000; i++) {
        vector.push_back(Element(i));
        // cout << vector[i].value() << " ";
    }

    // cout << endl;
}

void init_list(list<Element> &list) {
    list.clear();

    cout << "Init list.. " << endl;

    for (int i = 0; i < 10000000; i++) {
        list.push_back(Element(i));
        // cout << list[i].value() << " ";
    }

    // cout << endl;
}


void vector_test() {
    vector<Element> vector;
    ConcreteVisitorTest visitor(2);


    // Simple thread test

    init_vector(vector);

    auto start_time = chrono::high_resolution_clock::now();

    cout << "Simple thread.. " << endl;

    for (int i = 0; i < vector.size(); i++) {
        visitor.visit(vector[i]);
        // cout << vector[i].value() << " ";
    }

    // cout << endl;

    auto end_simple_time = chrono::high_resolution_clock::now();

    cout << "Simple thread done in "
         << chrono::duration_cast<chrono::milliseconds>(end_simple_time - start_time).count()
         << " ms" << endl;





    // OpenMP test

    init_vector(vector);

    auto start_openmp_time = chrono::high_resolution_clock::now();

    cout << "OpenMP.. " << endl;

    #pragma omp parallel for
    for (int i = 0; i < vector.size(); i++) {
        visitor.visit(vector[i]);
        // cout << vector[i].value() << " ";
    }

    // cout << endl;

    auto end_openmp_time = chrono::high_resolution_clock::now();

    cout << "OpenMP done in "
         << chrono::duration_cast<chrono::milliseconds>(end_openmp_time - start_openmp_time).count()
         << " ms" << endl;
}


void list_test() {
    list<Element> list;
    ConcreteVisitorTest visitor(2);


    // Simple thread test

    init_list(list);

    auto start_time = chrono::high_resolution_clock::now();

    cout << "Simple thread.. " << endl;

    auto iter = list.begin();

    for (; iter != list.end(); iter++) {
        visitor.visit(*iter);
        // cout << list[i].value() << " ";
    }

    // cout << endl;

    auto end_simple_time = chrono::high_resolution_clock::now();

    cout << "Simple thread done in "
         << chrono::duration_cast<chrono::milliseconds>(end_simple_time - start_time).count()
         << " ms" << endl;





    // OpenMP test

    init_list(list);

    auto start_openmp_time = chrono::high_resolution_clock::now();

    cout << "OpenMP.. " << endl;

    std::vector<Element*> elements;
    for (iter = list.begin(); iter != list.end(); ++iter)
        elements.push_back(&(*iter));

    #pragma omp parallel for
    for (int i = 0; i < elements.size(); i++) {
        visitor.visit(*(elements[i]));
        // cout << list[i].value() << " ";
    }

    // cout << endl;

    auto end_openmp_time = chrono::high_resolution_clock::now();

    cout << "OpenMP done in "
         << chrono::duration_cast<chrono::milliseconds>(end_openmp_time - start_openmp_time).count()
         << " ms" << endl;
}

int main(int argc, char const *argv[]) {

    vector_test();

    list_test();

    return 0;
}