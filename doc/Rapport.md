# Productions

* Scripts d'analyse : `script/analyse_time.sh̀̀̀`, `script/memory_usage.sh̀̀̀` et `script/massif_analyser.awk`
* Tracés des graphes : `script/plot_time.py` 
* Programme de test : `build/parallelize_test`

# Comparaison de librairies de parallélisation

## Open MP

### Avantages
* Simple d'utilisation
 
### Inconvénients
* Difficile d'utiliser des itérateurs

## Intel TBB

### Avantages
* Container thread-safe
* Tasks

### Inconvénients

# Comparaisons de consommation mémoire

## Script pour l'analyse mémoire

Programme `script/memory_usage.sh [N]`

Utilise valgrind.  

Lance N fois le programme `mesher_roi` sur a.poly avec l'option -s i (i allant de 1 à N).
Sauvegarde le résutat dans `build/memory_usage_mesher_roi_N_DATE`.  
Les fichiers massif.out.*.i correspondent aux résultats détaillés de l'analyse mémoire avec i niveaux de raffinement.  
Le fichier memory_usage contient le pic de mémoire utilisée pour chaque niveaux de raffinement.

## Passage de list à deque

Comparaison entre vector / list / deque : https://baptiste-wicht.com/posts/2012/12/cpp-benchmark-vector-list-deque.html  

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


# Comparaisons de temps d'exécutions

## Script pour l'analyse

Programme `script/analyse_time.sh`

Lance le programme `build/parallelize_test` pour différentes valeurs d'éléments (100 à 10^9), et différents nombre de threads (8 à 1).
Sauvegarde le résultat dans `script/analyse_time_[DATE]/time_size_[NBELEM]_thread_[NBTRHEAD]`

## Affichage des résultats :

Fichier script/plot_time.py à executer avec python3 (pyplot)

3 types de graphes :

* elements_times : Pour chaque méthode, un graphique présentant le temps d'execution en fonction du nombre d'élements, pour chaque nombre de threads.

* thread_times : Pour chaque nombre d'élement, un graphique présentanat les temps d'éxécution en fonction des threads, pour chaque méthode

* strong_scaling : Idem que précedemment, mais avec l'accélération (TempsAvec1Thread / TempsAvecNThread)

# Analyse du programme Mesher_roi

## Classe Mesher

### Variables de classe

```cpp
vector<MeshPoint> points;
vector<Quadrant> Quadrants;
set<QuadEdge> QuadEdges;
list<RefinementRegion *> regions;
```

### Fonction generateMesh (ou refineMesh)

* generateMesh : Si on part d'un niveau de raffinement égal à 0
* refineMesh : Si un maillage à un certain niveau de raffinement est déjà existant (commence donc à ce niveau de raffinement)

#### Paramètres

	- input : Ref sur Polyline
	- rl : niveau de raffinement
	- name : le nom de l'output (si spécifié)
	- all_reg : liste des ptr des régions à raffiner (Surface et/ou allRegion et/ou Boundary(pour refineMesh))
	- decoration : bool pour savoir si on ajoute des décorations au fichier VTK
	
Paramètres supplémentaires pour refineMesh :

	- roctli : liste des quadrants de départ
	- gt : transformation géométrique
	- minrl : le raffinement minimum des quadrants du maillage
	- omaxrl : le raffinement maximum des quadrants du maillage


#### Algo

* rotateGridMesh(input, all_reg, gt) : applique si nécessaire la rotation dans gt sur input
* generateGridMesh(input) : Crée et stocke dans ```Quadrants``` les Quadrant de départ (à partir de la polyline)
* generateQuadTreeMesh(rl, input, all_reg, name) : Algo de base pour la génération du maillage (Pour ```refineMesh```, appelle splitQuadrants)

Suite identique pour refineMesh et generateMesh??? à re-checker.

### Fonction generateQuatdreeMesh

Fonction appelé lorsque l'on souhaite raffiner un poly.

#### Paramètres
    - rl : niveau de raffinement souhaité
    - input : la polyline
    - all_reg : liste de RafinementRegion
    - name : nom de l'output

#### Algo
##### Init
        - temp_Quadrants : list des Quadrant actuel
        - new_Quadrants : list Quadrant vide
        - new_pts : liste Point3D vide
        - sv : SplitVisitor
        - i : niveau actuel de raffinement
##### Corps
Faire 

* vide `new_pts`
* tant que `tmp_Quadrants` n'est pas vide
    * `iter` = `tmpQuadrants[0]`
    * `to_refine` <- false
    * computeMaxDistance, ie met à jour max_dis dans `iter`
    * on regarde si `iter` a besoin d'être raffiné (on maj `to_refine` en conséquence). Cela dépend de la statégie donc de `all_reg` (raffinement pour tous et/ou ceux en intersection avec la polyline)
    * si le quadrant ne doit pas être raffiné
        * on ajoute le quadrant à la fin de `new_Quadrants`
    * sinon si il doit être raffiné
        * creation `vector<vector<Point3D>> clipping_coords` et ajout au SplitVisitor
        * creation `vector<vector<unsigned>> split_elements` et ajout au SplitVisitor : les nouveaux elements créés par le SV
        * iter->accept(sv) : Le SV regarde si le quadrant doit être splitté, et sauvegarde les nouveaux elements dans split_elements
        * si `inter_edges` est vide (ie inner quad == la polyline n'intersecte pas et l'élément est à l'intérieur)
            * on raffine le quadrant, donc pour chaque élément splittés  
                * on créée un nouveau Quadrant à partir des éléments splittés de raffinement `i+1`
                * on l'ajoute à la fin de la liste `new_Quadrants`
        * sinon (la polyline intersecte ou est dehors de la polyline)
            * pour chaque element splitté
                * on créée un nouveau Quadrant à partir des éléments splittés de raffinement `i+1`
                * Création de l'IntersectionVisitor : A faire avant ?
                * Ajout des infos du quadrant à l'IntersectionVisitor
                * si le quadrant est en intersection avec une input face
                    * on l'ajoute à la fin de la liste `new_Quadrants`
                * sinon on doit regarder s'il est à l'exterieur ou à l'interieur
                    * si il est à l'intérieur, on l'ajoute à la liste `new_Quadrants`
    * on enleve `iter` de `tmp_Quadrants`
* swap `tmp_Quadrants`, `new_Quadrants`
* insert dans `points` le contenu de `new_pts`
* pré-incrémente `i` 
  
tant que `new_pts` n'est pas vide