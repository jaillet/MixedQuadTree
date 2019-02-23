# Comparaison lib de parallélisation

## Open MP

### Avantages
* Simple d'utilisation
 
### Inconvénients

## Intel TBB

### Avantages
* Container thread-safe

### Inconvénients

# Script pour la consommation mémoire

Utilise valgrind.  

Lance N fois le programme mesher_roi sur a.poly avec l'option -s i (i allant de 1 à N).
Sauvegarde le résutat dans build/memory_usage_mesher_roi_N_DATE.  
Les fichiers massif.out.*.i correspondent aux résultats détaillés de l'analyse mémoire avec i niveaux de raffinement.  
Le fichier memory_usage contient la pic de mémoire utilisée pour chaque niveaux de raffinement.

# Passage de list à deque

Utilisation mémoire légèrement plus importante.

## Utilisation mémoire avec list

Command : mesher_roi -p ../data/a.poly -s N

|N  | Peak memory usage|
|---|------------------|
|1  | 115.816 Ko|
|2  | 129.616 Ko|  
|3  | 166.096 Ko|  
|4  | 265.664 Ko|  
|5  | 534.928 Ko|  
|6  | 1.18233 Mo|  
|7  | 2.57091 Mo|  
|8  | 5.66521 Mo|  
|9  | 12.0442 Mo|  
|10 | 24.0443 Mo|  


Command : mesher_roi -p ../data/a.poly -a N

|N  | Peak memory usage|
|---|------------------|
1  | 115.816 Ko
2  | 129.616 Ko
3  | 166.28 Ko
4  | 265.56 Ko
5  | 564.976 Ko
6  | 1.55675 Mo
7  | 5.06608 Mo
8  | 17.9289 Mo
9  | 68.0225 Mo
10 | 263.253 Mo


## Utilisation mémoire avec deque


Command : mesher_roi -p ../data/a.poly -s N

|N  | Peak memory usage|
|---|------------------|
|1  | 119.712 Ko|
|2  | 134.368 Ko|
|3  | 169.856 Ko|
|4  | 264.68 Ko|
|5  | 543.824 Ko|
|6  | 1.23563 Mo|
|7  | 2.75078 Mo|
|8  | 6.0956 Mo|
|9  | 12.9715 Mo|
|10 | 25.9963 Mo|

Command : mesher_roi -p ../data/a.poly -a N

|N  | Peak memory usage|
|---|------------------|
1  | 119.792 Ko
2  | 134.368 Ko
3  | 169.968 Ko
4  | 264.68 Ko
5  | 570.192 Ko
6  | 1.57748 Mo
7  | 5.12051 Mo
8  | 18.2691 Mo
9  | 68.9012 Mo
10 | 267.182 Mo




