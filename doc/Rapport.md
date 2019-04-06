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

##### Premier niveau

<pre>
Faire  
    Vider <b>new_pts</b>  
    Boucle tant que <b>tmp_Quadrants</b> n'est pas vide sur niveau 2  
    Swap <b>tmp_Quadrants</b> et <b>new_Quadrants</b>
    Si <b>new_pts</b> vide
        break
    Ajout dans <b>points (var classe Mesher)</b> le contenu de <b>new_pts</b>
    Incrémente <b>i</b>
Tant que <b>new_pts</b> n'est pas vide 
</pre>

##### Deuxième niveau

<pre>
Tant que <b>tmp_Quadrants</b> n'est pas vide
    <b>iter</b> <- début d'itérateur de <b>tmp_Quadrants</b>
    
    computeMaxDistance sur <b>iter</b> (Quadrant)
        Lecture <b>points</b> (var classe Mesher)
        Lecture <b>pointindex</b> (var classe Quadrant)
        Lecture <b>sub_elements</b> (var classe Quadrant)
        Ecriture <span style="color:orange"><b>max_dis</b></span> (var classe Quadrant)
        
    <b>to_refine</b> <- (*all_reg.begin())->intersectsQuadrant(points, *iter)
        TODO
    
    S'il n'est pas <b>to_refine</b>
        Ajout de <b>iter</b> dans <b>new_quadrants</b>
    Sinon
        Init ref <b>inter_edges</b> avec <b>iter->intersected_edges</b> (var classe Quadrant)
        Init <b>clipping_coords</b> et set à <b>sv.clipping</b> (SplitVisitor)
        Init <b>split_elements</b> et set à <b>sv.new_eles</b> (SplitVisitor)
        
        Applique <b>sv</b> (SplitVisitor) sur <b>iter</b> (Quadrant)
            Lecture <b>iter.pointindex</b>
            Insertion <b>sv.new_eles</b> (ref vers <b>split_elements</b>)
            Lecture <b>sv.points</b> (ref vers <b>points</b> var classe Mesher)          
            Insertion / Lecture <b>sv.new_pts</b> (ref vers <b>new_pts</b>) 
            Insertion / Suppression / Lecture <b>sv.edges</b> (ref vers <b>QuadEdges</b> var classe Mesher)         
            Insertion <b>sv.clipping</b> (ref vers <b>clipping_coords</b>)
        
        Si <b>inter_edges</b> (ref vers <b>iter->inteersected_edges</b>) vide
            Pour tout <b>split_elements</b>
                Init Quadrant <b>o</b> à partir de <b>split_elements</b>
                Insertion <b>new_Quadrants</b> du quadrant <b>o</b> 
        Sinon
            Pour tout <b>split_elements</b>
                Init Quadrant <b>o</b> à partir de <b>split_elements</b>
                Init iv (IntersectionVisitor)
                Set <b>iv.ply</b> (ref vers <b>input</b>)
                Set <b>iv.edges</b> (ref vers <b>inter_edges</b>)
                Set <b>iv.coords</b> (ref vers <b>clipping_coords</b> du split element)
                
                Applique <b>iv</b> (IntersectionVisitor) sur <b>iter</b> (Quadrant)
                    Insertion / Lecture <b>o.intersected_edges</b>
                    Lecture <b>iv.ply.mVertices</b> (ref vers input.mVertices)
                    Lecture <b>iv.edges</b> (ref vers iter->intersected_edges)
                    Lecture <b>iv.ply.mEdges</b> (ref vers input.mEdges)
                    
                    call intersectsEdge
                        Lecture <b>iv.coords</b> (ref vers <b>clipping_coords</b> du split element)
                        Lecture <b>iv.ply.mVertices</b> (ref vers input.mVertices)
                        Lecture <b>iv.ply.mEdges</b> (ref vers input.mEdges)
                        
                        call computePosition
                            Lecture <b>iv.ply.mEdges</b> (ref vers input.mEdges)
                        
                
                        
                        TODO cpliGeneralCase
                    
                Si retour fonction vrai
                    Insertion <b>new_Quadrants</b> du quadrant <b>o</b>  
                Sinon
                    Apelle fonction isItIn
                        TODO
                                            
                    Si retour fonction vrai
                        Insertion <b>new_Quadrants</b> du quadrant <b>o</b>  
                
    Retire <b>iter</b> de <b>tmp_quadrants</b>  
    
</pre>

#### Algo 2

##### First level

<pre>
For each level until rl  
    Clear <b>new_pts</b>  
    While until <b>tmp_Quadrants</b> is empty (go to Second level)
    Swap <b>tmp_Quadrants</b> and <b>new_Quadrants</b>
    If <b>new_pts</b> empty
        break
    Insert in <b>points (var class Mesher)</b> the content of <b>new_pts</b>
</pre>

##### Second level

<pre>
While until <b>tmp_Quadrants</b> is empty
    <b>iter</b> <- begin iterator of <b>tmp_Quadrants</b>
    
    For each RafinementRegion in <b>all_reg</b>
        <b>region_rl</b> <- getRafinementLevel of RafinementRegion
        If <b>region_rl</b> is lower than <b>i</b> (raffinement level) then continue
        If <b>region_rl</b> is equal or lower than <b>iter</b>  getRafinementLevel then continue
        If <b>reg_iter</b> intersectsQuantrant on points en <b>iter</b> then <b>to_refine</b> <- true
        
    If not <b>to_refine</b>
        Insert <b>iter</b> in <b>new_quadrants</b>
    Else (to refine)
        Init ref <b>inter_edges</b> with <b>iter->intersected_edges</b> (var class Quadrant)
        Init <b>qrl</b> with <b>iter</b> rafinement level
        Init <b>clipping_coords</b> and set <b>sv.clipping</b> (SplitVisitor)
        Init <b>split_elements</b> and set <b>sv.new_eles</b> (SplitVisitor)
        
        Apply <b>sv</b> (SplitVisitor) on <b>iter</b> (Quadrant)
            Read <b>iter.pointindex</b>
            Insert <b>sv.new_eles</b> (ref <b>split_elements</b>)
            Read <b>sv.points</b> (ref <b>points</b> var class Mesher)          
            Insert / Read <b>sv.new_pts</b> (ref <b>new_pts</b>) 
            Insert / Remove / Read <b>sv.edges</b> (ref <b>QuadEdges</b> var class Mesher)         
            Insert <b>sv.clipping</b> (ref <b>clipping_coords</b>)
        
        If <b>inter_edges</b> (ref <b>iter->inteersected_edges</b>) empty
            For each <b>split_elements</b>
                Init Quadrant <b>o</b> based on <b>split_elements</b>
                Insert <b>o</b> quadrant in <b>new_Quadrants</b>
        Else
            For each <b>split_elements</b>
                Init Quadrant <b>o</b> based on <b>split_elements</b>
                Init iv (IntersectionVisitor)
                Set <b>iv.ply</b> (ref <b>input</b>)
                Set <b>iv.edges</b> (ref <b>inter_edges</b>)
                Set <b>iv.coords</b> (ref <b>clipping_coords</b> of split element)
                
                Apply <b>iv</b> (IntersectionVisitor) on <b>iter</b> (Quadrant)
                    Insert / Read <b>o.intersected_edges</b>
                    Read <b>iv.ply.mVertices</b> (ref input.mVertices)
                    Read <b>iv.edges</b> (ref iter->intersected_edges)
                    Read <b>iv.ply.mEdges</b> (ref input.mEdges)
                    
                    Call intersectsEdge
                        Read <b>iv.coords</b> (ref <b>clipping_coords</b> of split element)
                        Read <b>iv.ply.mVertices</b> (ref input.mVertices)
                        Read <b>iv.ply.mEdges</b> (ref input.mEdges)
                        
                        Call computePosition
                            Read <b>iv.ply.mEdges</b> (ref input.mEdges)
                        
                        Call cpliGeneralCase
                            TODO 
                    
                If return true
                    Insert <b>o</b> quadrant in <b>new_Quadrants</b>
                Else
                    Call fonction isItIn
                        TODO
                                            
                    If return true
                        Insert <b>o</b> quadrant in <b>new_Quadrants</b>  
                
    Remove <b>iter</b> from <b>tmp_quadrants</b>  
    
</pre>

##### Concurrent access for second level

Insert in <b>new_quadrants</b>

SplitVisitor (set class var)

Insert / Read <b>sv.new_pts</b> (ref <b>new_pts</b>) 
Insert / Remove / Read <b>sv.edges</b> (ref <b>QuadEdges</b> var class Mesher)  

Remove <b>iter</b> from <b>tmp_quadrants</b>  

In SplitVisitor, insert new points in <b>new_pts</b> and use the indice of this points to update mid point of edges.  
So, this should be fixed for multi-threading because multiple thread can insert points at the same time.


