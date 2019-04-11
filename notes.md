#PROCHAINE REUNION

13h LUNDI 1 AVRIL

# SplitVisitor
Lorsque visit() accès concurrent :

* new_pts (list<Point3D>) : var de la classe pour stocker les nouveaux pts créés 
* edges (set<QuadEdge>) : suppression et insertion d'éléments 
* points (vector<MeshPoint>) : accès en lecture 
* new_eles (vector<vector<unsigned int> >) : var de la classe pour stocker les nouveaux pts créés
* clipping (vector<vector<Point3D> >) : var de la classe pour stocker les nouveaux pts créés

Problème potentiel, ce sert de points (ref vers Mesher.points) et new_eles (ref vers data de boucle) pour calculer le nombre de point existant et attribué le numéro suivant lors de la création d'un nouveau.

# IntersectionsVisitor
Lorsque visit() accès concurrent :

* select_edges (bool) :  
true, check avec des points spécifiques (-> coords & edges)  
false, check avec tous les points (-> points)
* ply (Polyline) : accès en lecture / copies
* quandrant->intersected_edges (list<unsigned int>) : accès lecture et écriture
* edges (list<unsigned int>) : accès en lecture
* coords (vector<Point3D>) : accès en lecture
* points (vector<MeshPoint>) : accès en lecture

# Dans la première boucle

Variable **i** correspond au niveau de raffinement actuel.  
Variable **new_pts** permet de sortir de la première boucle :
- A chaque début d'itération elle est clear
- Est utilisé par le **SplitVisitor** qui lui ajoute des éléments.
- A chaque fin d'itération son contenu est ajouté à la variable **points (var de classe)**

**SplitVisitor** accès concurrent : 
* new_pts (list<Point3D>) : réinit après traitement en entier du quadrant (donc plusieurs itérations donc ajout de pts)
* edges (set<QuadEdge>) : container global
* points (vector<MeshPoint>) : container global
* new_eles (vector<vector<unsigned int> >) : réinit avant chaque appel
* clipping (vector<vector<Point3D> >) : réinit avant chaque appel

## Premier niveau

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

## Deuxième niveau

<pre>
Tant que <b>tmp_Quadrants</b> n'est pas vide
    <b>iter</b> <- début d'itérateur de <b>tmp_Quadrants</b>
    
    computeMaxDistance sur <b>iter</b> (Quadrant)
        Lecture <b>points</b> (var classe Mesher)
        Lecture <b>pointindex</b> (var classe Quadrant)
        Lecture <b>sub_elements</b> (var classe Quadrant)
        Ecriture <b>max_dis</b> (var classe Quadrant)
        
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
                    Lecture <b>o.intersected_edges</b>
                    Lecture <b>input.mVertices</b>
                    TODO
                    
                Si retour fonction vrai
                    Insertion <b>new_Quadrants</b> du quadrant <b>o</b>  
                Sinon
                    Apelle fonction isItIn
                        TODO
                                            
                    Si retour fonction vrai
                        Insertion <b>new_Quadrants</b> du quadrant <b>o</b>  
                
    Retire <b>iter</b> de <b>tmp_quadrants</b>  
    
</pre>



# Compte rendu réunion

Regarder si on peut écrire avec plusieurs threads en même temps.
  
Effectuer des tests inteltbb et openmp avec nb thread fixe.  
Regarder usage mémoire et temps.  
Faire plusieurs execution + moyenne + tracé courbe.


Faire un container set avec plusieurs thread qui le modifie et s'en sert pour le modifier dans les visiteurs.  
Et tester. Comment vérifier les valeurs ? 
Regarder si c'est possible de tracer avec inteltbb / openmp.
Vérifier si option pour executer en sequentiel puis parallele.




Analyser le code existant voir là ou ça coince.  


# Analyse temps / structures / threads

# Analyse performance et parallélisation de set avec accès concurrent

## Méthode de test

Reproduction d'une boucle similaire à la première boucle de generateQuadtreeMesh dans Mesher.cpp.  


Algo + Reduction.  
Chaque thread a edge et new pts local (new pts indices commence au precedent).  
Reduction avec map pour chaque thread pour inserer new pts et edges dans global.

ECRITE ALGO ET LEUR ENVOYER
IMPLEMENTER EN OPENMP ET INTELTBB


# TBB FOR ALGO
For params : a -s 15

Parallel version with 8 threads :
         * level 0 in 1 ms
         * level 1 in 0 ms
         * level 2 in 0 ms
         * level 3 in 0 ms
         * level 4 in 1 ms
         * level 5 in 3 ms
         * level 6 in 5 ms
         * level 7 in 10 ms
         * level 8 in 31 ms
         * level 9 in 49 ms
         * level 10 in 87 ms
         * level 11 in 188 ms
         * level 12 in 345 ms
         * level 13 in 676 ms
         * level 14 in 1287 ms
- number of Points : 856771
- number of QuadEdge : 444155 

Sequentail version :
         * level 0 in 35 ms
         * level 1 in 0 ms
         * level 2 in 0 ms
         * level 3 in 1 ms
         * level 4 in 2 ms
         * level 5 in 7 ms
         * level 6 in 14 ms
         * level 7 in 30 ms
         * level 8 in 62 ms
         * level 9 in 128 ms
         * level 10 in 265 ms
         * level 11 in 543 ms
         * level 12 in 1140 ms
         * level 13 in 2318 ms
         * level 14 in 4727 ms
- number of Points : 856758
- number of QuadEdge : 444155