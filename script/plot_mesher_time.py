import glob
import matplotlib.pyplot as plt
import numpy as np

folderName = "analyse_mesher_1557901/"
maxRL = 10
nbThreads = 12
nbTries = 10

#dataA = {} # list indexed by nbThread of dictionnary
# exemple : dataA = [ { "customReductionTBB" : [ list indexed by level, contains meanTime] }  ]
#dataS = {}


def getData(option):
    data = {}
    for nbThread in range(1,nbThreads+1):
    
        data[nbThread] = {}
        files = glob.glob(folderName + "analyse_mesher_" + option + "/thread_" + str(nbThread) + "/*")

        for fileStr in files:
            file = open(fileStr, 'r')

            methodName = fileStr[:-4].split("/")[-1]

            #Sum of all time for each level
            data[nbThread][methodName] = [0 for i in range(maxRL)]

            for number, line in enumerate(file):
                line = line.split(' ')
                if line[0] == "Level":
                    data[nbThread][methodName][int(line[1])] += int(line[3])

        #Calculate mean
        for method in data[nbThread]:
            data[nbThread][method] = [x / nbTries for x in data[nbThread][method]]

    return data

dataA = getData('a')
dataS = getData('s')

"""
def plot_threads_times():
    counter = 1

    fig = plt.figure()
    fig.suptitle("Graph of average thread execution time for different number of elements", fontsize=16)

    #Fullscreen:
    mng = plt.get_current_fig_manager()
    #mng.resize(*mng.window.maxsize())

    #List of (pyplot lines, label) -> for legend
    #do it once (same legend for everyone??)
    legendLines = []
    legendLabel = []
    done = False

    for nb_elements, nb_elements_values in sorted(data.items()):

        #100 elem = 0..
        if nb_elements == 100:
            continue

        #plt.figure(counter)
        ax = plt.subplot(2, 3, counter)

        plt.title(str(nb_elements) + " elements.")
        plt.ylabel("Execution time (ms)")
        plt.xlabel("Number of threads")

        for type, type_values in nb_elements_values.items():
            #type_value : {nbThread : timeValue}
            #type : OpenMP_vector, IntelTBB_deque, etc..

            x = list(type_values.keys())
            x.sort()

            print(x)
            print(type_values)

            if len(x) == len(type_values.items()):
                lst = list(type_values.items())
                lst.sort(key=takeFirst)
                values = list()

                for val in lst:
                    values.append(val[1])
                l = plt.plot(x, values, '-x', label=type, markersize=15)

                #For figure legend:
                if not done:
                    legendLines.append(l)
                    legendLabel.append(type)

            for nb_threads in x:
                print(nb_elements, type, nb_threads, type_values[nb_threads])

        #plt.legend()
        counter += 1
        done = True

    ax.legend(bbox_to_anchor=(1.15, 0.5), loc='lower left', borderaxespad=0.)
    plt.show()

"""