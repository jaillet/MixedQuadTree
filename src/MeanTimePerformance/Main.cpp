#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <map>
#include <regex>

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    //std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);

    FILE* pipe = popen(cmd, "r");
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        result += buffer.data();
    }
    auto returnCode = pclose(pipe);
    std::cout << "EXIT CODE : " << returnCode << std::endl;
    return result;
}

void endMsg(){
    std::cout << std::flush;
    std::cerr << "use: ./meanTime number [t]\n";
    std::cerr << "where:\n";
    std::cerr << "number is the number of times we execute mesher_roi program\n";
    std::cerr << "Then at least one of the following list:\n";
    std::cerr << "t : to compute the mean of the process Transition Patterns\n";
    //TODO add doc
}

void fillMapToEvaluate(char arg, std::map<std::string, float>& map) {

    //Always look total :
    map["Generation done"] = 0.f;

    //TODO add options for every process...
    //Or check all every time?
    switch(arg) {
        case 't': {
            map["Transition Patterns"] = 0.f;
            break;
        }
        default: break;
    }
}

void addCorrespondingTime(const std::string& ouputOfMesher, std::map<std::string, float>& map) {

    //Parse output and find lines that contain a word that is a key of map
    //Then store corresponding time (after ' in ')

    std::string processPattern; //Example "Transition Patterns"
    std::smatch match;
    for(auto& kv : map) {
        processPattern = kv.first;
        std::regex rgx(processPattern + " in ([0-9]+)");
        if (std::regex_search(ouputOfMesher.begin(), ouputOfMesher.end(), match, rgx)) {
            map[processPattern] += std::stoi(match[1].str().c_str());
        }
        else {
            //No match, weird?
            std::cerr << "No matching found for regex : " << processPattern << " in ([0-9]+)" << "in the string :\n" << ouputOfMesher << std::endl;
        }
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
    std::map<std::string, float> mapProcessToTime;
    for(unsigned i = 0; i < arguments.size(); ++i) {
        arg = arguments[i];
        fillMapToEvaluate(arg, mapProcessToTime);
    }

    /*Debugging
    for (std::vector<std::string>::iterator i = mapProcessToTime.begin(); i != mapProcessToTime.end(); ++i)
    {
        std::cout << *i << std::endl;
    }
    std::cout << number << std::endl;
    */

    for (int i = 0; i < number; ++i)
    {
        std::string output = exec(cmd);
        addCorrespondingTime(output, mapProcessToTime);
    }

    //Divide to get the mean
    float nb = (float) number;
    for(auto& kv : mapProcessToTime) {
        mapProcessToTime[kv.first] /= nb;
    }

    std::cout << "Test ended. After " << number << " executions, here are the results :\n";
    for(auto& kv : mapProcessToTime) {
        std::cout << " * " << kv.first << " in " << kv.second << " ms" << std::endl;
    }




    return 0;
}