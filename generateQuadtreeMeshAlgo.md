# Explication fonction Mesher::generateQuatdreeMesh

Fonction appelé lorsque l'on souhaite raffiner un poly.

## Paramètres
    - rl : niveau de raffinement souhaité
    - input : la polyline
    - all_reg : liste de RafinementRegion
    - name : nom de l'output
    
## Variable de classe
- points : ??
    
## Algo
### Init
        - temp_Quadrants : list des Quadrant actuel
        - new_Quadrants : list Quadrant vide
        - new_pts : liste Point3D vide
        - sv : SplitVisitor ???
        - i : niveau actuel de raffinement
### Corps
Faire 

* vide **new_pts**
* tant que **tmp_Quadrants** n'est pas vide
    * **iter** <- tmp_Quadrants[0]
    * **to_refine** <- false
    * computeMaxDistance, ie met à jour max_dis dans **iter**
    * on regarde si **iter** a besoin d'être raffiné (on maj **to_refine** en conséquence). Cela dépend de la statégie donc de **all_reg** (raffinement pour tous et/ou ceux en intersection avec la polyline)
    * si le quadrant ne doit pas être raffiné
        * on ajoute le quadrant à la fin de **new_Quadrants**
    * sinon si il doit être raffiné
        * ??? SplitVisito
        * si **inter_edges** est vide (ie inner quad)
            * on raffine le quadrant, donc pour chaque élément splittés  
                * on créée un nouveau Quadrant à partir des éléments splittés de raffinement **i+1**
                * on l'ajoute à la fin de la liste **new_Quadrants**
        * sinon
            * pour chaque element splitté
                * on créée un nouveau Quadrant à partir des éléments splittés de raffinement **i+1**
                * ??? IntersectionsVisitor
                * si le quadrant est en intersection avec une input face
                    * on l'ajoute à la fin de la liste **new_Quadrants**
                * sinon on doit regarder s'il est à l'exterieur ou à l'interieur
                    * si il est à l'intérieur, on l'ajoute à la liste **new_Quadrants**
    * on enleve **iter** de **tmp_Quadrants**
* swap **tmp_Quadrants**, **new_Quadrants**
* insert dans **points** le contenu de **new_pts**
* pré-incrémente **i** 
  
tant que **new_pts** n'est pas vide 
    
      
## Questions
Dans **Quadrant**, pointindex est correspond à ?  
Dans **Quadrant**, que fais computeMaxDistance ?
Quand est maj **new_pts** ? Dans les parties visiteur apperemment

Dans **generateQuadtreeMeshAlgo**, intéret de faire continue et pop dans if ?
```
while
    if 
        ---
        tmp_Quadrants.pop()
        continue
    else
        ---
    tmp_Quandrants.pop()
end while
```