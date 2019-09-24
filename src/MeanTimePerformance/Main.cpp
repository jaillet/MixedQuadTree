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

    std::string cmds(cmd);
    cmds = cmds + " 2> /dev/null";

    FILE* pipe = popen(cmds.c_str(), "r");
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        result += buffer.data();
    }
    auto returnCode = pclose(pipe);
    //std::cout << "EXIT CODE : " << returnCode << std::endl;
    return result;
}

void endMsg(){
    std::cout << std::flush;
    std::cerr << "use: ./meanTime number [options]\n";
    std::cerr << "where:\n";
    std::cerr << "number is the number of times we execute mesher_roi program\n";
    std::cerr << "options are the options of mesher_roi (except -p and -u), default is : -a 7\n";
}

void fillMapToEvaluate(std::map<std::string, float>& map) {
    //Initialize the map to store meantime for each process :
    map["generateGridMesh"] = 0.f;
    map["Refine Quad"] = 0.f;
    map["Balanced mesh"] = 0.f;
    map["Transition Patterns"] = 0.f;
    map["generateQuadtreeMesh"] = 0.f;
    map["detectFeatureQuadrants"] = 0.f;
    map["linkElementsToNodes"] = 0.f;
    map["detectInsideNodes"] = 0.f;
    map["WriteQuadtreeMesh"] = 0.f;
    map["ProjectCloseToBoundary"] = 0.f;
    map["RemoveOnSurfaceSafe"] = 0.f;
    map["linkElementsToNodes"] = 0.f;
    map["ShrinkToBoundary"] = 0.f;
    map["ApplySurfacePatterns"] = 0.f;
    map["Generation done"] = 0.f;
}

bool addCorrespondingTime(const std::string& ouputOfMesher, std::map<std::string, float>& map) {

    //Parse output and find lines that contain a word that is a key of map
    //Then store corresponding time (after ' in ')

    std::string processPattern; //Example "Transition Patterns"
    std::smatch match;

    bool added = true;

    for(auto& kv : map) {
        processPattern = kv.first;
        std::regex rgx(processPattern + " in ([0-9]+)");
        if (std::regex_search(ouputOfMesher.begin(), ouputOfMesher.end(), match, rgx)) {
            map[processPattern] += std::stoi(match[1].str().c_str());
            if (!added) {
                std::cerr << "WARNING: added one for process " << processPattern << " but not for the processPatern before it......\n";
                std::cerr << "Mean will be wrong !!!!" << std::endl; 
            }
        }
        else {
            //No match, weird?
            std::cerr << "No matching found for regex : " << processPattern << " in ([0-9]+)" << "in the string :\n" << ouputOfMesher << std::endl;
            added = false;
        }
    }

    return added;

}

int main(int argc, char const *argv[])
{

    if (argc < 2) {
        endMsg();
        return 1;
    }

    std::string options = "-a 9";
    const int numberOfExecutions = std::atoi(argv[1]);
    
    if (argc > 2) {
        //Options for mesher_roi are provided
        //We want to concatenate all the remaining arg
        options = "";
        for (int i = 2; i < argc; i++) {
            options += argv[i];
            options += " ";
        }
    }

    std::string cmd = "./mesher_roi -p ../data/a.poly " + options + " -u ../output/a";

    std::cout << "Running next command " << numberOfExecutions << " times : \n\t" << cmd << std::endl;

    //Initialize map
    std::map<std::string, float> mapProcessToTime;
    fillMapToEvaluate(mapProcessToTime);

    //Executions
    int NUMBER_GOOD = 0;
    for (int i = 0; i < numberOfExecutions; ++i)
    {
        std::string output = exec(cmd.c_str());
        if (addCorrespondingTime(output, mapProcessToTime)) NUMBER_GOOD++;
        //std::cout << "\r" << "Number of executions : " << i+1; Don't know why this prints AFTER for loop
    }
    std::cout << std::endl;

    //Divide to get the mean
    float nb = (float) NUMBER_GOOD;
    for(auto& kv : mapProcessToTime) {
        mapProcessToTime[kv.first] /= nb;
    }

    std::cout << "Test ended. After " << numberOfExecutions << " executions, with " << NUMBER_GOOD << " good executions, here are the results :\n\n";
    for(auto& kv : mapProcessToTime) {
        std::cout << " * " << kv.first << " in " << kv.second << " ms" << std::endl;
    }

    return 0;
}