#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <vector>

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

void endMsg(){
    std::cout << std::flush;
    std::cerr << "use: ./meanTime number\n";
    std::cerr << "where:\n";
    std::cerr << "number is the number of times we execute mesher_roi program\n";
}

void fillVectorToEvaluate(char arg, std::vector<std::string>& vec) {
    switch(arg) {
        case 't': {
            vec.push_back("Transition Patterns");
            break;
        }
        default: break;
    }
}

int main(int argc, char const *argv[])
{
    
    const char* cmd = "./mesher_roi -p ../data/a.poly -a 5 -u ../output/a";

    if (argc < 3) {
        endMsg();
        return 1;
    }

    int number = std::atoi(argv[1]);

    std::string arguments = argv[2];
    char arg;
    std::vector<std::string> processToEvaluate;
    for(unsigned i = 0; i < arguments.size(); ++i) {
        arg = arguments[i];
        fillVectorToEvaluate(arg, processToEvaluate);
    }

    /*Debugging
    for (std::vector<std::string>::iterator i = processToEvaluate.begin(); i != processToEvaluate.end(); ++i)
    {
        std::cout << *i << std::endl;
    }
    std::cout << number << std::endl;
    */


    return 0;
}