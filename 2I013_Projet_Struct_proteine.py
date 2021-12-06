import math
import os
import glob


print("\nPour executer la fonction principale, veuillez entrer la fonction prog () dans la console.")
print("\nPour supprimer tous les motifs, veuillez entrer clear_motif_files () dans la console.")

### Partie 1 : Initialisation des alphabets et des paves
### ----------------------------------------------------


### Alphabet : Angles.


def lire_fichier_pos (fichier) :
    
    """ char -> tuple (list [float] ** 3, int)
    
        Retourne les positions en 3 dimensions des atomes c alpha du fichier .pdb et le nombre de TER - 1. """
    
    f = open (fichier)
    line = f.readline ()
    lx = []
    ly = []
    lz = []
    mark = 0
    nbter = 0
    
    while line :
        if line[:6] == "ATOM  " and line [13:16] == "CA " :
            if mark != int (line [22:26]) :
                mark = int (line [22:26])
            else :
                line = f.readline ()
                continue
            
            lx.append (float (line[31:38]))
            ly.append (float (line[39:46]))
            lz.append (float (line[47:54]))
        elif line[:3] == "TER" :
            lx.append (None)
            ly.append (None)
            lz.append (None)
            nbter += 1
        
        line = f.readline ()
    
    f.close ()
    
    return (lx, ly, lz, nbter - 1)


def calcul_angles_plans (lx, ly, lz) :

    """ list [float] * list [float] * list [float] -> (list [float], list [int])
    
        Retourne l'angle oriente des c alpha a partir de son precedent, de lui-meme et de ses 2 suivants et retourne les delimitations des TER. """
        
    lnx = []
    lny = []
    lnz = []
    lpm = []
    angles = []
    
    for i in range(1, len(lx) - 1) :
        if lx[i - 1] == None or lx[i] == None or lx[i + 1] == None :
            lnx.append (None)
            lny.append (None)
            lnz.append (None)
            lpm.append (None)
        else :
            lnx.append((ly[i] - ly[i - 1]) * (lz[i + 1] - lz[i - 1]) - (lz[i] - lz[i - 1]) * (ly[i + 1] - ly[i - 1]))
            lny.append((lz[i] - lz[i - 1]) * (lx[i + 1] - lx[i - 1]) - (lx[i] - lx[i - 1]) * (lz[i + 1] - lz[i - 1]))
            lnz.append((lx[i] - lx[i - 1]) * (ly[i + 1] - ly[i - 1]) - (ly[i] - ly[i - 1]) * (lx[i + 1] - lx[i - 1]))
            lpm.append(- lnx[i - 1] * lx[i] - lny[i - 1] * ly[i] - lnz[i - 1] * lz[i])
    
    nb = 0
    maxisi = []
    
    for i in range(1, len(lnx)):
        if lnx[i - 1] == None or lnx[i] == None or lx[i + 2] == None :
            if nb != 0 :
                maxisi.append (nb)
                nb = 0
            
            continue
            
        alpha = math.acos ((lnx[i - 1] * lnx[i] + lny[i - 1] * lny[i] + lnz[i - 1] * lnz[i]) / (math.sqrt (lnx[i - 1] * lnx[i - 1] + lny[i - 1] * lny[i - 1] + lnz[i - 1] * lnz[i - 1]) * math.sqrt (lnx[i] * lnx[i] + lny[i] * lny[i] + lnz[i] * lnz[i]))) * 180 / math.pi
        
        if lnx[i - 1] * lx[i + 2] + lny[i - 1] * ly[i + 2] + lnz[i - 1] * lz[i + 2] + lpm[i - 1] < 0 :
            angles.append (- alpha)
        else:
            angles.append (alpha)
        
        nb += 1
	
    return (angles, maxisi)


### Alphabet : Distances.
    

def distance3d (x1, x2, x3, y1, y2, y3) :
    
    """ float ** 6 -> float 
    
        Retourne la distance euclidienne de x (x1, x2, x3) a y (y1, y2, y3). """
        
    return math.sqrt (((x1 - y1) ** 2) + ((x2 - y2) ** 2) + ((x3 - y3) ** 2))


def list_distances (lx, ly, lz, pas) :
    
    """ list [float] * list [float] * list [float] * int -> list [float]
    
        Retourne la liste des distances entre un residu situe a i et un a i + pas. """
    
    distances = []
    lxf = []
    lyf = []
    lzf = []
    i = 1
    
    while i < len (lx) - 3 :
        if lx[i + 2] == None :
            i = i + 2 + 1
        else :
            lxf.append (lx[i])
            lyf.append (ly[i])
            lzf.append (lz[i])
            
        i += 1
    
    for i in range (len (lxf) - pas) :
        distances.append (distance3d (lxf [i], lyf [i], lzf [i], lxf [i + pas], lyf [i + pas], lzf [i + pas]))
    
    return distances


### Paves pour les angles et les distances.
    

def creer_paves (texte, discretisation) :
    
    """ list [float] * float -> list [set {float}]
    
        Retourne les paves des donnes des atomes c alpha avec une discretisation. """
    
    mi = texte[0]
    ma = texte[0]
    
    for k in texte[1:] :
        if k < mi :
            mi = k
        elif k > ma :
            ma = k
	
    d = discretisation
    nb = 1
    
    while d < 1:
        d = d * 10
        nb *= 10
    
    length = int ((int (ma * nb) // d) - (int (mi * nb) // d) + 1)
    s = []
    
    for i in range (length) :
        s.append (set ())
        
    for l in texte :
        s[int (int (l * nb) // d - (int (mi * nb) // d))].add (l)
        
    return s


def paves_agrandis (paves, degenerescence) :
    
    """ list [set {float}] * float -> list [set {float}]
        PS : degenerescence de 1 revient a ne pas appliquer cette fonction
    
        Retourne les paves en tenant compte de la degenerescence. """
     
    np = []        
        
    for i in range (len (paves)) :
        ens = set ()
        
        for k in range(degenerescence) :
            ens = ens | paves[(i + k) % len(paves)]
            
        if ens != set ():
            np.append (ens)    
        
    return np


### Partie 2 : Algorithme de Triade non relationnel
### -----------------------------------------------


def Init (s, paves) :
    
    """ list [char] * list [set {char}] -> list [set {int}]
        s : chaine a etudiee
        
        Retourne le tableau v qui contient les indices des elements de s dans paves. """
    
    v = []
    
    for lettre in s :
        ve = set ()
        
        for indice, ensemble in enumerate (paves) :
            if lettre in ensemble :
                ve.add (indice)
            
        v.append (ve)
        
    return v


def Dual (v) :
    
    """ list [set {int}] -> list [list [int]]
        v : tableau issu de Init
        
        Retourne le dual p de v sous forme de pile. """
        
    p = [[]]

    for indice, ensemble in enumerate (v) :      
        for k in ensemble :
            if k + 1 != len (p) :
                for i in range (k + 1 - len (p)) :
                    p.append ([])
            
            p[k].append (indice)
    
    return p


def BuildQa (v, p, pas) :
    
    """ list [set {int}] * list [list [int]] * int -> list [list [int]]
        v : tableau issu de Init, p : pile dual de v
        
        Retourne le tableau q issu du tri de sorte que p[i] est dans q[v[p[i] + m]]. """

    copy_p = p.copy ()
    q = []
    card = []
    
    for i in range (len (p)) :
        q.append ([])
        card.append (0)
        
    for sous_tab in copy_p :
        sous_tab.reverse ()
        
        for k in sous_tab :
            if k + pas < len (v) :
                for l in v[k + pas] :
                    q[l].append (k)
                    card[l] += 1
        
        for ci, c in enumerate (card) :
            if c < 2 and c != 0 :
                card[ci] = 0
                
                for j in range (c):
                    q[ci].pop ()
            if card[ci] != 0 :
                q[ci].append (-1)
                
            card[ci] = 0        
                    
    for qi in range (len (q)) :
        if q[qi] != [] :
            q[qi].pop ()
            
    nq = []
    
    for tq in q:
        if tq != []:
            nq.append (tq)
        
    return nq    


def FiltreQa (qa):
    
    """ list [list [int]] -> list [list [int]]
    
        Retourne qa filtre. """
    
    ll = []
    
    for sl in qa:
        nsl = [[]]
        
        for k in sl:
            if k == -1:
                nsl.append ([])
            else:
                nsl[len (nsl) - 1].append (k)
        
        aj = []
        
        for i1 in range (len (nsl)):
            b = False
            
            if len (aj) == 0:
                aj.append (nsl [i1])
                continue
            
            for taj in aj:
                if set (nsl [i1]) <= set (taj):
                    b = True
                    break
            
            if b:
                continue
            else:
                raj = aj.copy ()
                
                for traj in raj:
                    if set (traj) < set (nsl[i1]):
                        aj.remove (traj)
                        
                aj.append (nsl [i1])
        
        ll.append (aj)
    
    if len (ll) != 0:
        nl = [ll[0]]
        nv = [kl for kl in ll[0]]
    else:
        nl = []
        nv = []
    
    for i in range (1, len (ll)) : 
        nl.append ([])
        
        for tl in ll [i] :
            pnv = nv.copy ()
            b = True
            
            for tnv in pnv :
                if set (tnv) <= set (tl) :
                    nv.remove (tnv)
                    
                    for j in range (len (nl)):
                        if tnv in nl[j]:
                            nl[j].remove (tnv)
                  
                elif set (tnv) > set (tl):
                    b = False 
                    break
            
            if b:
                nl[len (nl) - 1].append (tl)
                nv.append (tl)
                    
    while [] in nl:
        nl.remove ([])
    
    rl = []
    
    for tll in nl :
        rl.append ([])
        
        for ktll in tll :
            for k in ktll :
                rl[len (rl) - 1].append (k)
            
            rl[len (rl) - 1].append (-1)
        
        rl[len (rl) - 1].pop ()
    
    return rl     


def BuildP (qb) :
    
    """ list [list [int]] -> list [set {int}]
        q : tableau trie issu de Build
        
        Retourne le tableau p issu de Qb. """
        
    v = [set ()]
    i = 0
    rq = [sq[::-1] for sq in qb]
    
    for sous_tab in rq :
        for k in sous_tab :
            if k == -1 :
                i += 1
                v.append (set ())
                
                continue
            
            v[i].add (k)
        
        i += 1    
        v.append (set ())
    
    v.pop ()    
        
    return v


def FilterP (p) :
    
    """ list [list [int]] -> list [list [int]]
        p : tableau a filtrer
        
        Retourne le tableau p filtre. """
    
    s = set ();    
    
    for i in range (len (p)) :
        for j in range (i + 1, len (p)) :
            if p[i] == p[j] :
                s.add (j)
            if p[i] < p[j] :
                s.add (i)
            if p[i] > p[j] :
                s.add (j)
    
    copy_p = []
                
    for l in range (len (p)) :
        if l in s :
            continue
        
        copy_p.append (p[l])
        
    return copy_p


def VectPf (pf) :
    
    """ list [set {char}] -> list [set {int}]
    
        Retourne le tableau v issu de pf avec le meme principe que Init. """
 
    v = []
    

    for indice, ensemble in enumerate (pf) :
        for k in ensemble :
            if k + 1 > len(v) :
                for i in range (k + 1 - len(v)) :
                    v.append (set ())
            
            v[k].add (indice) 

    return v   
   
    
def Exe_non_rela (s, maxisc, lnbter, paves, k, quorum) :
    
    """ list [char] * list [int] * list [int] * list [set {char}] * int * int -> tuple (list [set {int}], int)
            
        Retourne les positions des (max)k-motifs repetes et de leurs extensions et retourne la taille maximum possible d'un motif si k = -1. (non relationnel) """

    v = Init (s, paves)
    p = Dual (v)
    q = BuildQa (v, p, 1)
    q = FiltreQa (q)
    p = BuildP (q)
    p = FilterP (p)
    
    if k == -1 :
        p = cut (p, maxisc, lnbter, 2, quorum)
    
    if len (p) == 0 :
        return (p, 0)
    
    if k != -1 :
        lpuiss2 = [2 ** l for l in range (1, int (math.log2 (k)))]
        
        for i in lpuiss2 :
            v = VectPf (p)
            p = Dual (v)
            q = BuildQa (v, p, i)
            q = FiltreQa (q)
            p = BuildP (q)
            p = FilterP(p)
            
        ppre = p
        kp = k
    else :
        cptpuiss2 = 1
        
        while len (p) != 0 :
            ppre = p
            cptpuiss2 *= 2
            v = VectPf (p)
            p = Dual (v)
            q = BuildQa (v, p, cptpuiss2)
            q = FiltreQa (q)
            p = BuildP (q)
            p = FilterP(p)
            p = cut (p, maxisc, lnbter, cptpuiss2 * 2, quorum)
            
        kp = cptpuiss2
            
    return (ppre, kp)


### Partie 2 bis : algorithme de Triade relationnel
### -----------------------------------------------
   
    
### Fonction outil pour BuildQb.
    

def relations_pos_conti (ld, pavesd):
    
    """ tuple (float) * list [set {float}] -> list [set {float}]
    
        Retourne les relations entre deux positions contigues. """
        
    l = []
    
    for k in ld:
        l.append (set ())
        
        for i in range (len (pavesd)):
            if k in pavesd[i]:
                l[len (l) - 1].add (i)
    
    return l    


### Fonction a executer apres FiltreQa et avant BuildP.
    

def BuildQb (pos_conti, qa) :
     
    """ list [set {int}] * list [list [int]] * int -> list [list [int]]
    
        Retourne qb a partir de qa filtree et des relations entre deux positions contigues. """

    copy_qa = qa.copy ()
    qb = []
    card = []
    
    for i in range (len (qa)) :
        qb.append ([])
        card.append (0)
        
    for sous_tab in copy_qa :
        
        for k in sous_tab[::-1] :
            if k != -1:
                for j in pos_conti[k] :
                    if j < len (qb) :
                        qb[j].append (k)
                        card[j] += 1
            if k == -1 or k == sous_tab[0] :
                for ci, c in enumerate (card) :
                    if c < 2 and c != 0 :
                        for b in range (c) :
                            qb[ci].pop ()
                    if qb[ci] != [] and qb[ci][len (qb[ci]) - 1] != -1 :
                        qb[ci].append (-1)
                
                    card[ci] = 0
    
    nqb = []
    
    for l in range (len (qb)) :
        if qb[l] == [] :
            continue
        else:
            qb[l].pop ()
            nqb.append (qb[l])
        
    return nqb


def Exe_rela (s, maxisc, lnbter, lx, ly, lz, pavesa, discretisationd, degenerescenced, k, quorum) :
    
    """ list [char] * list [int] * list [int] * list [char] * list [char] * list [char] * list [set {char}] * float * int * int * int -> tuple (list [set {int}], int)
            
        Retourne les positions des (max)k-motifs repetes et de leurs extensions et retourne la taille maximum possible d'un motif si k = -1. (relationnel) """

    ld =  list_distances(lx, ly, lz, 1)
    pavesd = creer_paves (ld, discretisationd)
    pavesd = paves_agrandis (pavesd, degenerescenced)
    v = Init (s, pavesa)
    p = Dual (v)
    qa = BuildQa (v, p, 1)
    qa = FiltreQa (qa)
    pos_conti = relations_pos_conti (ld, pavesd)
    qb = BuildQb (pos_conti, qa)
    p = BuildP (qb)
    p = FilterP (p)
    
    if k == -1 :
        p = cut (p, maxisc, lnbter, 2, quorum)
    
    if len (p) == 0 :
        return (p, 0)
    
    if k != - 1 :
        for i in range (k - 2):
            ld =  list_distances(lx, ly, lz, i + 2)
            pavesd = creer_paves (ld, discretisationd)
            pavesd = paves_agrandis (pavesd, degenerescenced)
            v = VectPf (p)
            p = Dual (v)
            qa = BuildQa (v, p, i + 2)
            qa = FiltreQa (qa)
            pos_conti = relations_pos_conti (ld, pavesd)
            qb = BuildQb (pos_conti, qa)
            p = BuildP (qb)
            p = FilterP (p)
            
        ppre = p
        kp = k
    else :
        cptnnr = 1
        
        while len (p) != 0 :
            ppre = p
            cptnnr += 1
            ld =  list_distances(lx, ly, lz, cptnnr)
            pavesd = creer_paves (ld, discretisationd)
            pavesd = paves_agrandis (pavesd, degenerescenced)
            v = VectPf (p)
            p = Dual (v)
            qa = BuildQa (v, p, cptnnr)
            qa = FiltreQa (qa)
            pos_conti = relations_pos_conti (ld, pavesd)
            qb = BuildQb (pos_conti, qa)
            p = BuildP (qb)
            p = FilterP (p)
            p = cut (p, maxisc, lnbter, cptnnr + 1, quorum)
        
        kp = cptnnr
    
    return (ppre, kp)


### Partie 3 : Lecture des resultats de l'algorithme de Triade par creation de fichiers .pdb
### ----------------------------------------------------------------------------------------


### Lecture du proteine a etudier.
    

def lire_fichier (fichier) :
    
    """ char -> tuple (char, set {int : list [char]})
    
        Retourne le titre HEADER et le dictionnaire du fichier en associant les numeros des residus et leur lignes. """
        
    d = {}
    f = open (fichier, "r")
    line = f.readline ()
    i = -1
    
    if line :
        title = line
        line = f.readline ()
        
    pll = False
    pird = 0
    ecart_hetatm = 0     ### Compteur de decalage dans le cas ou le fichier .pdb possede les chaines d'atomes ayant des coupures de HETATM. (c'est le cas pour 6j66.pdb)
    
    while line :
        if line [:4] == "ATOM" :
            ir = int (line [22:26])
            
            if ir - i > 1 and pll:
                ecart_hetatm += ir - (i + 1)
            
            if not pll :
                ecart_hetatm = 0
                pir = ir + pird
                pll = True
                
            if i != ir :
                d [ir - pir - ecart_hetatm] = [line]
                i = ir
            else :
                d [ir - pir - ecart_hetatm].append (line)
        elif line [:3] == "TER" :
            pll = False
            pird = - (ir - pir - ecart_hetatm) - 1
            line = f.readline ()
            continue
        
        line = f.readline ()
    
    f.close ()
    
    return (title, d)


### Creation des fichiers .pdb. (motifs)


def creer_fichier (fichier, new_fichier, ep, maxisc, k):
    
    """ char * char * set {int} * list [int] * int -> void
    
        Retourne le fichier .pdb representant le k-motif repete d'un ensemble ep (indices des extentions). """
        
    d = lire_fichier (fichier)
    g = open (new_fichier, "w")
    g.write (d[0])
    nb = ord (chr (33))
    inds = list (ep)
    inds.sort ()
    deci = 0
    cpti = 0
    
    for val in inds :
        while val >= maxisc[cpti] :
            cpti += 1
            deci += 3
        
        for i in range (k) :
            for line in d[1][val + 1 + deci + i] :
                line_modif = "".join ((line[:21], chr(nb), line[22:]))
                g.write (line_modif)
        
        nb += 1
        
        if val != inds [len (inds) - 1] :
            g.write ("TER\n")
            
        if nb == ord (chr (127)) :
            break
    
    g.write ("END\n")
    g.close ()
    
    return

 
### Fonction principale de creation de pdb (non relationnelle) pour 1 proteine.   
    

def main (fichier, discretisationa, degenerescencea, k) :
    
    """ char * float * int * int -> tuple (int, int)
    
        Retourne le nombre de fichiers crees, la taille maximum possible d'un motif si k = -1 et cree tous les fichiers representant les (max)k-motifs repetes en tenant compte des parametres. (non relationnel) """
    
    (lx, ly, lz, nbter) = lire_fichier_pos (fichier)
    lnbter = [nbter]
    tup = calcul_angles_plans (lx, ly, lz)
    s = tup[0]
    maxisi = tup[1]
    maxisc = [maxisi[0]]
    
    for i in range (1, len (maxisi)) :
        maxisc.append (maxisc[i - 1] + maxisi[i])
        
    pavesa = creer_paves (s, discretisationa)
    pavesa = paves_agrandis (pavesa, degenerescencea)
    (result, maxk) = Exe_non_rela (s, maxisc, lnbter, pavesa, k, 1)
    
    if k != - 1 :
        result = cut (result, maxisc, lnbter, maxk, 1)     ### Dans le cas ou une proteine ou les chaines de residus ne sont pas reliees en continues. (separees par des TER comme par exemple la proteine 3zwu.pdb)
        
    nb = 0
    
    for ens in result :
        if len (ens) != 0 :
            nb += 1
            new_fichier = "{}-motif{}.pdb".format (maxk, nb)
            creer_fichier (fichier, new_fichier, ens, maxisc, maxk)
    
    return (len (result), maxk)


### Fonction principale de creation de pdb (relationnelle) pour 1 proteine.
    

def main2 (fichier, discretisationa, degenerescencea, discretisationd, degenerescenced, k) :
    
    """ char * float * int * float * int * int -> tuple (int, int)
    
        Retourne le nombre de fichiers crees, la taille maximum possible d'un motif si k = -1 et cree tous les fichiers representant les (max)k-motifs repetes en tenant compte des parametres. (relationnel) """
    
    (lx, ly, lz, nbter) = lire_fichier_pos (fichier)
    lnbter = [nbter]
    tup = calcul_angles_plans (lx, ly, lz)
    s = tup[0]
    maxisi = tup[1]
    maxisc = [maxisi[0]]
    
    for i in range (1, len (maxisi)) :
        maxisc.append (maxisc[i - 1] + maxisi[i])
        
    pavesa = creer_paves (s, discretisationa)
    pavesa = paves_agrandis (pavesa, degenerescencea)
    (result, maxk) = Exe_rela (s, maxisc, lnbter, lx, ly, lz, pavesa, discretisationd, degenerescenced, k, 1)
    
    if k != - 1 :
        result = cut (result, maxisc, lnbter, maxk, 1)     ### Dans le cas ou une proteine ou les residus ne sont pas reliees en continue. (separees par des TER comme par exemple la proteine 3zwu.pdb)

    nb = 0
    
    for ens in result :
        if len (ens) != 0 :
            nb += 1
            new_fichier = "{}-motif{}.pdb".format (maxk, nb)
            creer_fichier (fichier, new_fichier, ens, maxisc, maxk)
    
    return (len (result), maxk)


### Partie 3 : Travail sur une liste de proteine
### --------------------------------------------


def creer_fichier2 (list_fichiers, maxisc, lnbter, new_fichier, ep, k):
    
    """ list [char] * list[int] * list [int] * char * set {int} * int -> void
    
        Retourne le fichier .pdb representant k-motif repetes d'un ensemble ep (indices des extensions) en prenant compte plusieurs fichiers .pdb. """
        
    ld = []
    
    for c in list_fichiers :
        tup = lire_fichier (c)
        ld.append (tup[1])
    
    g = open (new_fichier, "w")
    nb = ord (chr (33))
    inds = list (ep)
    inds.sort ()
    cpt = 0
    dec = 0
    cpti = 0
    deci = 0
    lnbterc = lnbter.copy ()
    ival = 0
    nbe = 0
    endw = False
    
    while ival < len (inds) :
        while inds[ival] >= maxisc[cpt + cpti] :
            if lnbterc[cpt] == 0 : 
                g.write ("END FILE {}\n".format (list_fichiers[cpt].upper ()))
                endw = True
                dec = maxisc[cpt + cpti]
                cpt += 1
                deci = 0
                nbe = 0
            else :
                cpti += 1
                deci += 3
                lnbterc[cpt] -= 1
        
        if inds[ival] != inds[0] and not endw :
            g.write ("TER\n")
            
        endw = False
        
        for i in range (k) :
            for line in ld [cpt][inds[ival] + 1 + deci - dec + i] :
                line_modif = "".join ((line[:21], chr(nb), line[22:]))
                g.write (line_modif)
                
        nb += 1
        nbe += 1
        
        if nb == ord (chr (127)):
            break
        
        if nbe % ((ord (chr (127)) - ord (chr (33))) // len (list_fichiers)) == 0 and len (inds) >= ord (chr (127)) - ord (chr (33)) :
            while lnbterc[cpt] != 0 :
                cpti += 1
                lnbterc[cpt] -= 1
            
            while ival < len (inds) and inds[ival] < maxisc[cpt + cpti] :
                ival += 1
                
            ival -= 1
            g.write ("END OF {}\n".format (list_fichiers[cpt].upper ()))
            endw = True
            dec = maxisc[cpt + cpti]
            cpt += 1
            deci = 0
            nbe = 0
            
        ival += 1
                
    if cpt < len (list_fichiers) :        
        g.write ("END FILE {}\n".format (list_fichiers[cpt].upper ()))
    
    g.close ()
    
    return


### Fonction filtre pour eviter qu'un motif original ne deborde sur 2 proteines ou pour separer deux chaines de residus delimitees par des TER.


def cut (p, maxisc, lnbter, k, quorum):
    
    """ list[set {int}] * list [int] * list [int] * int * int -> list [set {int}]
    
        Retourne la liste des indices des motifs et de leurs extensions sans interception entre les fichiers .pdb et en respectant le quorum."""
    
    res = []
    
    for ens in p :
        res.append (set())
        tq = [0] * len (lnbter)
        
        for val in ens :
            debut = 0
            ind = 0
            cpti = 0
            
            while ind < len (lnbter) :
                cpti += lnbter[ind]
                
                if debut <= val and val < maxisc[ind + cpti] :
                    if val + k - 1 < maxisc[ind + cpti] :
                        res[len(res) - 1].add (val)
                        tq[ind] = 1
                    
                    break
                else:
                    debut = maxisc[ind + cpti]
                    
                ind += 1
                    
        nbq = 0
        
        
        for x in tq :
            nbq += x
                    
        if nbq < quorum or len (res[len (res) - 1]) <= 1:
            res.pop ()
                          
    return res 

### Fonction principale de creation de pdb (non relationnelle) pour une liste de proteines.


def main3 (list_fichiers, discretisationa, degenerescencea, k, quorum) :
    
    """ list[char] * float * int * int * int -> tuple (int, int)
    
        Retourne le nombre de fichiers crees, la taille maximum possible d'un motif si k = -1 et cree tous les fichiers representant les (max)k-motifs repetes en tenant compte des parametres. (non relationnel) """
    
    maxisi = []
    s = []
    lnbter = []
    
    for c in list_fichiers:
        (lx, ly, lz, nbter) = lire_fichier_pos (c)
        lnbter.append (nbter)
        tup = calcul_angles_plans (lx, ly, lz)
        s += tup[0]
        maxisi += tup[1]
        
    maxisc = [maxisi[0]] 
    
    for i in range (1, len (maxisi)) :
        maxisc.append (maxisc[i - 1] + maxisi[i])
    
    pavesa = creer_paves (s, discretisationa)
    pavesa = paves_agrandis (pavesa, degenerescencea)
    (result, maxk) = Exe_non_rela (s, maxisc, lnbter, pavesa, k, quorum)
    
    if k != -1 :
        result = cut (result, maxisc, lnbter, maxk, quorum)
    
    nb = 0
    
    for ens in result :
        if len (ens) != 0 :
            nb += 1
            new_fichier = "{}-motif{}.pdb".format (maxk, nb)
            creer_fichier2 (list_fichiers, maxisc, lnbter, new_fichier, ens, maxk)
    
    return (len (result), maxk)
    

### Fonction principale de creation de pdb (relationnelle) pour une liste de proteines.


def main4 (list_fichiers, discretisationa, degenerescencea, discretisationd, degenerescenced, k, quorum) :
    
    """ list[char] * int * int * int * int -> tuple (int, int)
    
        Retourne le nombre de fichiers crees, la taille maximum possible d'un motif si k = -1 et cree tous les fichiers representant les (max)k-motifs repetes en tenant compte des parametres. (relationnel) """
    
    maxisi = []
    s = []
    lx = []
    ly = []
    lz = []
    lnbter = []
    
    for c in list_fichiers:
        (px, py, pz, nbter) = lire_fichier_pos (c)
        lnbter.append (nbter)
        lx += px
        ly += py
        lz += pz
        tup = calcul_angles_plans (px, py, pz)
        s += tup[0]
        maxisi += tup[1]
     
    maxisc = [maxisi[0]] 
    
    for i in range (1, len (maxisi)) :
        maxisc.append (maxisc[i - 1] + maxisi[i])
    
    pavesa = creer_paves (s, discretisationa)
    pavesa = paves_agrandis (pavesa, degenerescencea)
    (result, maxk) = Exe_rela (s, maxisc, lnbter, lx, ly, lz, pavesa, discretisationd, degenerescenced, k, quorum)
    
    if k != -1 :
        result = cut (result, maxisc, lnbter, maxk, quorum)
    
    nb = 0
    
    for ens in result :
        if len (ens) != 0 :
            nb += 1
            new_fichier = "{}-motif{}.pdb".format(maxk, nb)
            creer_fichier2 (list_fichiers, maxisc, lnbter, new_fichier, ens, maxk)
    
    return (len (result), maxk)


### Partie 4 : Similitudes et fusions entre motifs 
### ----------------------------------------------
    

def fusion_motifs_similaires (result, seuil):
    
    """ list[set(int)] * float -> list[set(int)]
        seuil : entre 0 et 1 inclus  (distance de Jaccard)
    
        Retourne les positions des k-motifs repetes et de leurs extensions apres fusion de ces derniers. """
        
    new_result = []
    ajoute = False
    
    for ind1 in range (len (result)) :
        ens_app = result[ind1]
        test = False
        indt = 0
        
        while indt < len (new_result) and not test :
            if ens_app <= new_result[indt] :
                test = True
            indt += 1   
            
        if test:
            continue
        
        for ind2 in range (len (result)) :
            if ind1 == ind2 :
                continue
            
            ajoute = False
            
            if (1 - (len (result[ind1] & result[ind2]) / len (result[ind1] | result[ind2]))) < seuil :
                ind3 = 0
                
                while ind3 < len (new_result) and not ajoute :
                    if result[ind2] <= new_result[ind3] :
                        new_result[ind3] |= ens_app
                        ajoute = True
                    ind3 += 1
                    
                
                if not ajoute :
                    ens_app = ens_app | result[ind2]
                else:
                    break
                    
        if not ajoute:
            new_result.append (ens_app)
        
    return new_result


### Fonction similaire a main4 mais avec fusion des motifs repetes entre eux.


def main_final (list_fichiers, discretisationa, degenerescencea, discretisationd, degenerescenced, k, quorum, seuil) :
    
    """ list[char] * float * int * float * int * int * float -> tuple (int, int)
        seuil : entre 0 et 1 inclus  (distance de Jaccard)
    
        Retourne le nombre de fichiers crees, la taille maximum possible d'un motif si k = -1 et cree tous les fichiers representant les (max)k-motifs repetes avec fusion des motifs selon la distance de Jaccard (seuil). (relationnel) """
    
    maxisi = []
    s = []
    lx = []
    ly = []
    lz = []
    lnbter = []
    
    for c in list_fichiers:
        (px, py, pz, nbter) = lire_fichier_pos (c)
        lnbter.append (nbter)
        lx += px
        ly += py
        lz += pz
        tup = calcul_angles_plans (px, py, pz)
        s += tup[0]
        maxisi += tup[1]
     
    maxisc = [maxisi[0]] 
    
    for i in range (1, len (maxisi)) :
        maxisc.append (maxisc[i - 1] + maxisi[i])
        
    pavesa = creer_paves (s, discretisationa)
    pavesa = paves_agrandis (pavesa, degenerescencea)
    (result, maxk) = Exe_rela (s, maxisc, lnbter, lx, ly, lz, pavesa, discretisationd, degenerescenced, k, quorum)
    
    if k != -1 :
        result = cut (result, maxisc, lnbter, maxk, quorum)
    
    result = fusion_motifs_similaires (result, seuil)
    nb = 0
    
    for ens in result :
        if len (ens) != 0 :
            nb += 1
            new_fichier = "{}-motif{}.pdb".format (maxk, nb)
            creer_fichier2 (list_fichiers, maxisc, lnbter, new_fichier, ens, maxk)
    
    return (len (result), maxk)


### Annexe
### ------


def prog () :
    
    """ Fonction d'execution de toutes les fonctions principales en evitant les erreurs de frappes """
    
    c = int (input ("Votre choix :\n1 : algorithme a un pdb non relationnel\n2 : algorithme a un pdb relationnel\n3 : algorithme a plusieurs pdb non relationnel\n4 : algorithme a plusieurs pdb relationnel\n5 : algorithme a plusieurs pdb relationnel avec fusion des motifs\n\n"))                 
    
    while c not in {1, 2, 3, 4, 5} :
            c = int (input ("Votre choix :\n1 : algorithme a un pdb non relationnel\n2 : algorithme a un pdb relatioennel\n3 : algorithme a plusieurs pdb non relationnel\n4 : algorithme a plusieurs pdb relationnel\n5 : algorithme a plusieurs pdb relationnel avec fusion des motifs\n\n"))                 
    
    nb = 1
    
    if c in {3, 4, 5} :
        nbplu =  len ([k for k in os.listdir () if ".pdb" in k])
        print ("\nVotre nombre de fichiers .pdb ? Sachant qu'il y a :", nbplu, "fichiers .pdb dans le repertoire.")
        nb = int (input ())
        
        while nb < 1 :
            print ("\nErreur de nombre de fichiers. Recommencez.")
            print ("Votre nombre de fichiers .pdb ? Sachant qu'il y a :", nbplu, "fichiers .pdb dans le repertoire.")
            nb = int (input ())
            
    l = []
    
    for i in range (nb) :
        lp = [k for k in os.listdir () if ".pdb" in k]
        print ("\nLes fichiers .pdb presents dans le repertoire sont :")
        
        for nk in lp :
            print (nk)
            
        if nb == 1 :
            print ("\nLe nom du fichier .pdb ?")
        else :
            print ("\nLe nom du fichier .pdb numero", i + 1, "?")
        
        el = input ()
        
        while not os.path.exists (el) :
            print ("\nLe fichier .pdb n'existe pas. Recommencez.")
            
            if nb == 1 :
                print ("Le nom du fichier .pdb ?")
            else :
                print ("Le nom du fichier .pdb numero", i + 1, "?")            
            
            el = input ()
                
        l.append (el)
    
    discretisationa = float (input ("Discretisation des angles ?\n\n"))

    while discretisationa <= 0.0 :
        print ("\nDiscretisation invalide. Recommencez.")
        print ("Discretisation des angles ?")
        discretisationa = float (input ())
        
    degenerescencea = int ( float (input ("Degenerescence des angles ?\n\n")))

    while degenerescencea < 1 :
        print ("\nDegenerescence invalide. Recommencez.")
        print ("Degenerescence des angles ?")
        degenerescencea = int (float (input ()))
        
    if c in {2, 4, 5} :
        discretisationd = float (input ("Discretisation des distances ?\n\n"))

        while discretisationd <= 0.0 :
            print ("\nDiscretisation invalide. Recommencez.")
            print ("Discretisation des distances ?")
            discretisationd = float (input ())
        
        degenerescenced = int ( float (input ("Degenerescence des distances ?\n\n")))

        while degenerescenced < 1 :
            print ("\nDegenerescence invalide. Recommencez.")
            print ("Degenerescence des distances ?")
            degenerescenced = int (float (input ()))
            
    k = int (float (input ("Taille des motifs a trouver ? Si k = -1, alors on cherche les plus longs motifs possibles.\n\n")))
    
    if c in {1, 3} :
        if k > 2 and math.log2 (k) / int (math.log2 (k)) != 1 :
            k = 0
    
    while k < 2 and k != -1 :
        print ("\nTaille invalide. Recommencez.")
        print ("Taille des motifs a trouver ? Si k = -1, alors on cherche les plus longs motifs possibles.")
        k = int (float (input ()))
        
        if c in {1, 3} :
            if k > 2 and math.log2 (k) / int (math.log2 (k)) != 1 :
                k = 0
    
    if c in {3, 4, 5} :
        print ("\nQuorum ? Sachant que le nombre de .pdb lu est :", nb)
        quorum = int (float (input ()))
        
        while quorum <= 0 or quorum > len (l) :
            print ("\nQuorum invalide. Recommencez.")
            print ("Quorum ? Sachant que le nombre de .pdb lu est :", nb)
            quorum = int (float (input ()))
            
    if c == 5 :
        seuil = float (input ("Seuil de distance pour fusion ?\n\n"))
        
        while seuil < 0.0 or seuil > 1.0 :
            print ("\nSeuil invalide. Recommencez.")
            print ("Seuil de distance pour fusion ?\n\n")
            seuil = float (input ())
    
    
    if c == 1 :
        (t, maxk) = main (l[0], discretisationa, degenerescencea, k)
    elif c == 2 :
        (t, maxk) = main2 (l[0], discretisationa, degenerescencea, discretisationd, degenerescenced, k)
    elif c == 3 :
        (t, maxk) = main3 (l, discretisationa, degenerescencea, k, quorum)
    elif c == 4 :
        (t, maxk) = main4 (l, discretisationa, degenerescencea, discretisationd, degenerescenced, k, quorum)
    else :
        (t, maxk) = main_final (l, discretisationa, degenerescencea, discretisationd, degenerescenced, k, quorum, seuil)
        
    if k == -1 :
        print ("\nLes motifs les plus longs possibles sont de taille :", maxk)
    
    print ("\nLe nombre de motifs obtenus :")
    return t


### Fonction de suppression des motifs
    

def clear_motif_files () :
    
    """ void -> void
    
        Supprime tous les fichiers motifs. """
        
    l = glob.glob ("*motif*")
    
    for c in l:
            os.remove (c)
    
    return
