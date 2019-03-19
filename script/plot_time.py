import glob
import matplotlib.pyplot as plt
import numpy as np


def takeFirst(elem):
    return elem[0]

processor = ""
data = {}
data2 = {}

onlyfiles = glob.glob("analyse_time_2019-03-07_14:39:03/time_size*.txt")

for file in onlyfiles:
    file_object = open(file, 'r')

    nb_elements = 0
    nb_threads = 0

    for number, line in enumerate(file_object):

        if number >= 1 and number <= 6:
            #Lines for processor spec
            processor += line
        else :
            words = line.split(' ')

            if len(words) > 1:

                if words[0] == 'Launching':
                    nb_elements = int(words[3])
                    nb_threads = int(words[6])

                    if nb_elements not in data:
                        data[nb_elements] = {}

                    if nb_threads not in data2:
                        data2[nb_threads] = {}

                else:
                    type = words[0]
                    if type not in data[nb_elements]:
                        data[nb_elements][type] = {}

                    if type not in data2[nb_threads]:
                        data2[nb_threads][type] = {}

                    time = float(words[1])
                    
                    data[nb_elements][type][nb_threads] = time
                    data2[nb_threads][type][nb_elements] = time

print(processor)

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


def plot_elements_times():
    counter = 1

    fig = plt.figure()
    fig.suptitle("Graph of the average execution time of a structure according to the number of elements for a different number of threads", fontsize=16)

    #Fullscreen:
    mng = plt.get_current_fig_manager()
    #mng.resize(*mng.window.maxsize())

    #List of (pyplot lines, label) -> for legend
    #do it once (same legend for everyone??)
    legendLines = []
    legendLabel = []
    done = False

    for nb_threads, nb_threads_values in sorted(data2.items()):

        #plt.figure(counter)
        ax = plt.subplot(2, 3, counter)

        plt.title(str(nb_threads) + " threads.")
        plt.ylabel("Execution time (ms)")
        plt.xlabel("Number of elements")
        #plt.yscale('log')
        plt.xscale('log')

        for type, type_values in nb_threads_values.items():
            #type_values : {nbElem : timeValue}
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

                #x.pop(0)
                #values.pop(0)

                l = plt.plot(x, values, '-x', label=type, markersize=15)

                #For figure legend:
                if not done:
                    legendLines.append(l)
                    legendLabel.append(type)

            for nb_elements in x:
                print(nb_threads, type, nb_elements, type_values[nb_elements])

        #plt.legend()
        counter += 1
        done = True

        ax.legend()
    plt.show()

def plot_strong_scaling_speedup():
    #With the same value of nbElem, change the number of thread for one method
    #https://www.kth.se/blogs/pdc/2018/11/scalability-strong-and-weak-scaling/

    #speedup = t1 / tN, t1 = time with one thread, tN with N threads

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
        plt.ylabel("Speedup")
        plt.xlabel("Number of threads")

        for type, type_values in nb_elements_values.items():
            #type_value : {nbThread : timeValue}
            #type : OpenMP_vector, IntelTBB_deque, etc..

            #x = nbThreads
            x = list(type_values.keys())
            x.sort()

            print(type)
            print(x)
            print(type_values)

            #Time with one thread (change 2 to 1 after)
            t1 = type_values[2]

            speedup = list()
            for val in type_values.values():
                #Hack because val sometimes equals 0...
                if val == 0:
                    val = 0.01
                speedup.append(t1 / val)

            l = plt.plot(x, speedup, '-x', label=type, markersize=15)

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


#plot_elements_times()
#plot_threads_times()
plot_strong_scaling_speedup()
