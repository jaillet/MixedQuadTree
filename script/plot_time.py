import glob
import matplotlib.pyplot as plt
import numpy as np


def takeFirst(elem):
    return elem[0]


data = {}

onlyfiles = glob.glob("analyse_time_2019-03-07_14:39:03/time_size*.txt")

for file in onlyfiles:
    file_object = open(file, 'r')

    nb_elements = 0
    nb_threads = 0

    for line in file_object:
        words = line.split(' ')

        if len(words) > 1:

            if words[0] == 'Launching':
                nb_elements = int(words[3])
                nb_threads = int(words[6])

                if nb_elements not in data:
                    data[nb_elements] = {}

            else:
                #Now like : "simple_list"
                #And not like : "simple list"
                type = words[0]
                if type not in data[nb_elements]:
                    data[nb_elements][type] = {}

                time = float(words[1])
                #time = float(words.pop(0))
                
                data[nb_elements][type][nb_threads] = time

counter = 1

#print(data[100])

for nb_elements, nb_elements_values in data.items():

    plt.figure(counter)

    plt.title(str(nb_elements) + " elements.")
    plt.ylabel("Execution time")
    plt.xlabel("Number of threads")

    for type, type_values in nb_elements_values.items():
        #x = nombre de threads
        x = list(type_values.keys())
        x.sort()

        print(x)

        if len(x) == len(type_values.items()):
            lst = list(type_values.items())
            lst.sort(key=takeFirst)
            values = list()
            for val in lst:
                values.append(val[1])
            plt.plot(x, values, label=type)

        for nb_threads, time in type_values.items():
            print(nb_elements, type, nb_threads, time)

    plt.legend()
    plt.show(counter)

    counter += 1
