# SplitVisitor
Lorsque visit() accès concurrent : 
* new_pts (list<Point3D>) : var de la classe pour stocker les nouveaux pts créés 
* edges (set<QuadEdge>) : suppression et insertion d'éléments 
* points (vector<MeshPoint>) : accès en lecture 
* new_eles (vector<vector<unsigned int> >) : var de la classe pour stocker les nouveaux pts créés
* clipping (vector<vector<Point3D> >) : var de la classe pour stocker les nouveaux pts créés

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
