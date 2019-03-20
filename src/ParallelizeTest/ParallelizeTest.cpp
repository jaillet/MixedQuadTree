#include <thread>
#include <list>
#include <set>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <deque>
#include <map>
#include <omp.h>
#include <string.h>

#include "tbb/parallel_for_each.h"
#include "tbb/task_scheduler_init.h"

#include "TimeDecorator.hpp"
#include "CPUInfo.hpp"

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
class ElementBis;

class VisitorTest {
public:
    virtual void visit(Element &) = 0;
    virtual void visit(ElementBis &) = 0;
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

class ElementBis {
protected:
    bool even;
    set<int> values;

public:

    ElementBis() : even(true) { values.insert(0); }
    ElementBis(int val) : even(val % 2 == 0) { values.insert(val); };

    ElementBis(ElementBis & el) : even(el.isEven()) {
        for(auto ite = el.values.begin(); ite != el.values.end(); ite++) {
            values.insert(*ite);
        }
    };

    ElementBis(const ElementBis & el) : even(el.isEven()) {
        for(auto ite = el.values.begin(); ite != el.values.end(); ite++) {
            values.insert(*ite);
        }
    };

    virtual void accept(VisitorTest &visitor) {
        visitor.visit(*this);
    }

    void divide2() {
        int res = (*values.begin()) / 2;
        values.insert(res);
        even = res % 2 == 0;
    }

    const int value() const {
        return (*values.begin());
    }

    const bool isEven() const {
        return even;
    }

    bool operator <(ElementBis const& b)const {
        if (*values.begin() < *b.values.begin()) return true;
        else if (*values.begin() > *b.values.begin()) return false;
        else {
            if (values.size() > b.values.size()) return false;
            else if (values.size() < b.values.size()) return true;
            return false;
        }
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

    void visit(ElementBis & element) {};
};

class ConcreteSetVisitorTest : public VisitorTest {
private:
    set<ElementBis> & data;

public:
    ConcreteSetVisitorTest(set<ElementBis> & data) : data(data) {};

    void visit(Element &) {};

    virtual void visit(ElementBis & element) override {
        if (element.isEven()) {
            element.divide2();
            data.insert(element);
        }
    };
};

/* Global */

int NOMBRE_ELEM = 10000000;
int NOMBRE_ITER = 1;
int NOMBRE_THREAD = omp_get_max_threads();

vector<Element> elements;

map<string, double> time_average;

/*------------- Init -------------*/

void init_vector(vector<Element> &vector) {
    vector.clear();
    vector.reserve(NOMBRE_ELEM);

    //cout << "Init vector.. " << endl;

    //#pragma omp parallel for
    for (int i = 0; i < NOMBRE_ELEM; i++) {
        vector[i] = Element(i);
        //vector.emplace_back(i);
        //cout << vector[i].value() << " ";
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

void init_set(set<ElementBis> &set) {
    set.clear();

    for (int i = 1; i <= NOMBRE_ELEM; i++) {
        set.insert(ElementBis(i));
    }
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

void inteltbb_vector(vector<Element> &v, VisitorTest &visitor) {
    tbb::parallel_for_each(v.begin(), v.end(), [&] (Element &element) { visitor.visit(element);});
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
    //Very long solution, because each thread iterates over the list
    //Overhead strong, maybe OK with few elements and long process
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
    //Only with openmp 3.0 (tasks)
    //Only one thread who iterates, and creates one task for each element
}

void inteltbb_list(list<Element> &v, VisitorTest &visitor) {
    tbb::parallel_for_each(v.begin(), v.end(), [&] (Element &element) { visitor.visit(element);});
}

// DEQUE

void simple_deque(deque<Element> &l, VisitorTest &visitor) {
    auto iter = l.begin();
    for (; iter != l.end(); iter++) {
        visitor.visit(*iter);
    }
}

void openmp_deque(deque<Element> &l, VisitorTest &visitor) {
#pragma omp parallel for
    for (int i = 0; i < l.size(); i++) {
        visitor.visit(l[i]);
        // cout << l[i].value() << " ";
    }
}

void inteltbb_deque(deque<Element> &v, VisitorTest &visitor) {
    tbb::parallel_for_each(v.begin(), v.end(), [&] (Element &element) { visitor.visit(element);});
}

/*------------- Tests finaux -------------*/

void simple_vector_test(int nb_iteration) {
    vector<Element> v;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "Simple_vector ";
    for (int i = 0; i < nb_iteration; i++) {
        v.clear();
        v.assign(elements.begin(), elements.end());
        time = count_time(simple_vector, v, visitor);
        time_average["Simple vector"] += time;
    }
    time_average["Simple vector"] /= nb_iteration;
    cout << time_average["Simple vector"] << " ms." << endl;
}

void openmp_vector_test(int nb_iteration) {
    vector<Element> v;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "OpenMP_vector ";
    for (int i = 0; i < nb_iteration; i++) {
        v.clear();
        v.assign(elements.begin(), elements.end());
        time = count_time(openmp_vector, v, visitor);
        time_average["OpenMP vector"] += time;
    }
    time_average["OpenMP vector"] /= nb_iteration;
    cout << time_average["OpenMP vector"] << " ms." << endl;
}

void inteltbb_vector_test(int nb_iteration) {
    vector<Element> v;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "IntelTBB_vector ";
    for (int i = 0; i < nb_iteration; i++) {
        v.clear();
        v.assign(elements.begin(), elements.end());
        time = count_time(inteltbb_vector, v, visitor);
        time_average["IntelTBB vector"] += time;
    }
    time_average["IntelTBB vector"] /= nb_iteration;
    cout << time_average["IntelTBB vector"] << " ms." << endl;
}

void vector_test(int nb_iteration) {
    // Simple thread test
    //simple_vector_test(nb_iteration);

    // OpenMP test
    openmp_vector_test(nb_iteration);

    // IntelTBB test
    inteltbb_vector_test(nb_iteration);
}

void simple_list_test(int nb_iteration) {
    list<Element> l;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "Simple_list ";
    for (int i = 0; i < nb_iteration; i++) {
        l.clear();
        l.assign(elements.begin(), elements.end());
        time = count_time(simple_list, l, visitor);
        time_average["Simple list"] += time;
    }
    time_average["Simple list"] /= nb_iteration;
    cout << time_average["Simple list"] << " ms." << endl;
}

void openmp_list_copy_test(int nb_iteration) {
    list<Element> l;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "OpenMP_list_copy ";
    for (int i = 0; i < nb_iteration; i++) {
        l.clear();
        l.assign(elements.begin(), elements.end());
        time = count_time(openmp_list_copy, l, visitor);
        time_average["OpenMP list copy"] += time;
    }
    time_average["OpenMP list copy"] /= nb_iteration;
    cout << time_average["OpenMP list copy"] << " ms." << endl;
}

void openmp_list_copy2_test(int nb_iteration) {
    list<Element> l;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "OpenMP_list_copy2 ";
    for (int i = 0; i < nb_iteration; i++) {
        l.clear();
        l.assign(elements.begin(), elements.end());
        time = count_time(openmp_list_copy2, l, visitor);
        time_average["OpenMP list copy2"] += time;
    }
    time_average["OpenMP list copy2"] /= nb_iteration;
    cout << time_average["OpenMP list copy2"] << " ms." << endl;
}

void inteltbb_list_test(int nb_iteration) {
    list<Element> l;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "IntelTBB_list ";
    for (int i = 0; i < nb_iteration; i++) {
        l.clear();
        l.assign(elements.begin(), elements.end());
        time = count_time(inteltbb_list, l, visitor);
        time_average["IntelTBB list"] += time;
    }
    time_average["IntelTBB list"] /= nb_iteration;
    cout << time_average["IntelTBB list"] << " ms." << endl;
}

void list_test(int nb_iteration) {
    // Simple thread test
    //simple_list_test(nb_iteration);

    // OpenMP list copy
    openmp_list_copy_test(nb_iteration);

    // OpenMP list copy2
    openmp_list_copy2_test(nb_iteration);

    // IntelTBB test
    inteltbb_list_test(nb_iteration);
}

void simple_deque_test(int nb_iteration) {
    deque<Element> d;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "Simple_deque ";
    for (int i = 0; i < nb_iteration; i++) {
        d.clear();
        d.assign(elements.begin(), elements.end());
        time = count_time(simple_deque, d, visitor);
        time_average["Simple deque"] += time;
    }
    time_average["Simple deque"] /= nb_iteration;
    cout << time_average["Simple deque"] << " ms." << endl;
}

void openmp_deque_test(int nb_iteration) {
    deque<Element> d;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "OpenMP_deque ";
    for (int i = 0; i < nb_iteration; i++) {
        d.clear();
        d.assign(elements.begin(), elements.end());
        time = count_time(openmp_deque, d, visitor);
        time_average["OpenMP deque"] += time;
    }
    time_average["OpenMP deque"] /= nb_iteration;
    cout << time_average["OpenMP deque"] << " ms." << endl;
}

void inteltbb_deque_test(int nb_iteration) {
    deque<Element> d;
    ConcreteVisitorTest visitor(2);
    float time;

    cout << "IntelTBB_deque ";
    for (int i = 0; i < nb_iteration; i++) {
        d.clear();
        d.assign(elements.begin(), elements.end());
        time = count_time(inteltbb_deque, d, visitor);
        time_average["IntelTBB deque"] += time;
    }
    time_average["IntelTBB deque"] /= nb_iteration;
    cout << time_average["IntelTBB deque"] << " ms." << endl;
}

void deque_test(int nb_iteration) {
    // Simple thread test
    //simple_deque_test(nb_iteration);

    // OpenMP test
    openmp_deque_test(nb_iteration);

    // IntelTBB test
    inteltbb_deque_test(nb_iteration);
}

/* ############## TEST set ################### */

void simple_set(set<ElementBis> & elements) {
    set<ElementBis> result;
    set<ElementBis> new_els;


    do {
        new_els.clear();
        ConcreteSetVisitorTest visitor(new_els);

        auto ite = elements.begin();
        while (ite != elements.end()) {
            ElementBis el (*ite);

            el.accept(visitor);

            if (!(el).isEven()) {
                result.insert(el);
            }
            //cout << new_els.size() << endl;
            ite++;
        }

        std::swap(new_els, elements);

    } while (!new_els.empty());

    std::cout << result.size() << endl;

    int res = 0;
    for (auto ite = result.begin(); ite != result.end(); ite++) {
        cout << (*ite).value() << " ";
        res += (*ite).value();
    }
    cout << endl;

    cout << "Res : " << res << endl;
}

void simple_set_test(int nb_iteration) {
    set<ElementBis> elements;
    float time;

    cout << "Simple_set ";
    for (int i = 0; i < nb_iteration; i++) {
        init_set(elements);
        time = count_time(simple_set, elements);
        time_average["Simple set"] += time;
    }
    time_average["Simple set"] /= NOMBRE_ITER;
    cout << time_average["Simple set"] << " ms." << endl;

}

/*------------- Main -------------*/

int main(int argc, char const *argv[]) {

    if (argc > 1) {
        NOMBRE_ELEM = atoi(argv[1]);
    }

    if (argc > 2) {
        if (atoi(argv[2]) <= NOMBRE_THREAD && atoi(argv[2]) > 0) {
            NOMBRE_THREAD = atoi(argv[2]);
        } else {
            cerr << "Invalid number of threads or not supported by computer" << endl;
            exit(1);
        }
    }

    omp_set_num_threads(NOMBRE_THREAD);
    tbb::task_scheduler_init test(NOMBRE_THREAD);

    cout << "Launching tests with " << NOMBRE_ELEM << " elements and " << NOMBRE_THREAD << " threads" << endl;
    
    CPUInfo cinfo;
    //6 lines of info
    std::cout << "Processor : " << std::endl;
    std::cout << std::boolalpha; //Cout true for 1, false for 0
    std::cout << "CPU vendor = " << cinfo.vendor() << std::endl;
    std::cout << "CPU Brand String = " << cinfo.model() << std::endl;
    std::cout << "# of cores = " << cinfo.cores() << " ";
    std::cout << "# of logical cores = " << cinfo.logicalCpus() << std::endl;
    std::cout << "# of thread (std::thread) = " << std::thread::hardware_concurrency() << std::endl;
    std::cout << "Is CPU Hyper threaded = " << cinfo.isHyperThreaded() << std::endl;
    
    init_vector(elements);

    if (argc > 3) {
        if (strcmp(argv[3], "vector") == 0) vector_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "simple_vector") == 0) simple_vector_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "openmp_vector") == 0) openmp_vector_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "inteltbb_vector") == 0) inteltbb_vector_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "list") == 0) list_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "simple_list") == 0) simple_list_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "openmp_list_copy") == 0) openmp_list_copy_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "openmp_list_copy2") == 0) openmp_list_copy2_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "inteltbb_list") == 0) inteltbb_list_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "deque") == 0) deque_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "simple_deque") == 0) simple_deque_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "openmp_deque") == 0) openmp_deque_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "inteltbb_deque") == 0) inteltbb_deque_test(NOMBRE_ITER);
        else if (strcmp(argv[3], "simple_set") == 0) simple_set_test(NOMBRE_ITER);
        else {
            cerr << "Invalid argument" << endl;
            exit(1);
        }

    } else {
        vector_test(NOMBRE_ITER);

        list_test(NOMBRE_ITER);

        deque_test(NOMBRE_ITER);
    }

    return 0;
}