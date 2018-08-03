from fenics import *
import os
import errno
from datetime import datetime
from copy import copy
import numpy as np

__this_files_path = os.path.realpath(__file__)
__this_files_dir  = os.path.dirname(__this_files_path)

# Pfad zum Ordner mit den Gitter Dateien
DATA_DIR = os.path.abspath(os.path.join(__this_files_dir, 'Meshes'))


def get_vertex_normal(meshdata, facet_index):
    '''
    Berechnet die aeusseren (?) Normalen und gibt diese nach Numerierung der Facetten
    nur 2D moeglich bis jetzt
    Format: 1d array: Koordinaten des i'ten Normalenvektor sind 2i und 2i+1
    verschieden Interpolationen der Normalen auf Vertex durch Normale der umliegenden
    Facetten denkbar. (Wie finite diffrences stencils)
    '''

    normal = []

    # erste Normalen fuer erste Facette
    f = Facet(meshdata.mesh, facet_index[-1])
    f_norm_1 = f.normal().array()
    f = Facet(meshdata.mesh, facet_index[0])
    f_norm_2 = f.normal().array()

    # Interpolation der Facettennormalen
    normal_coords = 0.5 * (f_norm_1 + f_norm_2)
    normal.append(normal_coords[0])
    normal.append(normal_coords[1])

    for i in range(len(facet_index)-1):

        # hole Normalenkoordinaten der umliegenden Facetten von Vertex i+1
        # .normal gibt durch FEniCS einen 3D Vektor zurueck
        f = Facet(meshdata.mesh, facet_index[i])
        f_norm_1 = f.normal().array()
        f = Facet(meshdata.mesh, facet_index[i+1])
        f_norm_2 = f.normal().array()

        # Interpolation der Facettennormalen
        normal_coords = 0.5 * (f_norm_1 + f_norm_2)
        normal.append(normal_coords[0])
        normal.append(normal_coords[1])

    # korrigiere die Vorzeichen der Normalenvektoren, um aeussere Normalen zu erhalten


    return np.array(normal)


def order_boundary_facets(meshdata, bound_index):
    '''
     Gibt geordnete Liste der Boundary-Facetten an nach FEniCS-Facet Indizierung)
     (Ordnung haengt von bound_index ab!
    '''

    bound_facet_index = []

    for i in range(len(bound_index)-1):
        v = Vertex(meshdata.mesh, bound_index[i])

        # waehle Facette, die naechsten Vertex enthaelt
        for f in facets(v):
            if bound_index[i+1] in f.entities(0):
                bound_facet_index.append(f.index())


    # fuer den letzten Index
    v = Vertex(meshdata.mesh, bound_index[len(bound_index)-1])
    for f in facets(v):
        if bound_index[0] in f.entities(0):
            bound_facet_index.append(f.index())


    return bound_facet_index


def order_boundary_verts(meshdata):
    '''
     Gibt geordnete Indexliste der Boundary-Punkte am Interace zurueck
     Muessen nur einmal berechnet werden, da invariant bei Deformation (testen?)
    '''


    # Vertex-Indizes des Interface bekommen, bleiben invariant bei Deformation
    ind_boundary_facets = []
    for i in range(0, len(meshdata.boundaries)):
        if meshdata.boundaries[i] > 4:
            ind_boundary_facets.append(i)

    ind_boundary_facets = list(set(ind_boundary_facets))

    ind_boundary_vertices = []
    for c in cells(meshdata.mesh):
        for f in facets(c):
            if f.index() in ind_boundary_facets:
                for v in vertices(f):
                    ind_boundary_vertices.append(v.index())

    interface_indizes = list(set(ind_boundary_vertices))


    # ordne Interface Indizes der Connectivity nach, d.h. Nachbarn im Array sind Nachbarn im Interface
    ordered_interface_indizes = []
    ordered_interface_indizes.append(interface_indizes[0])

    i = 0
    while(len(ordered_interface_indizes) < len(interface_indizes)):
        v = Vertex(meshdata.mesh, ordered_interface_indizes[i])
        local_connectivity = []

        for f in facets(v):
            # kann hoechstens zwei Elemente in 2D haben, da 1dim Unterman.
            # entferne Vertex, an welchem wir uns bei Iterationsschritt befinden

            if ordered_interface_indizes[i] in f.entities(0):
                f_erased = f.entities(0).tolist()
                del f_erased[f_erased.index(ordered_interface_indizes[i])]
                local_connectivity.append(f_erased[0])

            if len(local_connectivity) >= 2: break

        if(i == 0): ordered_interface_indizes.append(local_connectivity[0])
            # im ersten Knoten hat man 2 moegliche Nachbarn

        else:
            # ein Nachbar ist bereits in der Liste, waehle den anderen
            if local_connectivity[0] in ordered_interface_indizes:
                ordered_interface_indizes.append(local_connectivity[1])

            else: ordered_interface_indizes.append(local_connectivity[0])

        i += 1

    return ordered_interface_indizes


def input_signal(meshdata, ampl, omega, eps, perimeter_val,  analytic = False):
    """
    Erzeugt eine Schwingung auf dem Interface, diese ist Element der Orthonormalbasis auf
    dem Shape-Tangentialraum. Man sieht nur endlich viele aufgeloest, je nach Gitterfeinheit;
    also genau die Anzahl an Knoten viel ONV

    kann auf jede periodische Funktion auf der geschlossenen Boundary verallgemeinert werden

    omega: frequency; ampl: amplitude
    
    Was ist der Zsmmhg zwischen diesen ONV und Eigenvektoren bzgl. Shape-Hessian?
    
    """

    # Calculate Perimeter of Interface
    V = FunctionSpace(meshdata.mesh, 'P', 1)
    ones = Function(V)
    ones.vector()[:] = 1.

    # if(analytic):
    #     perimeter = perimeter_val
    #
    # else:
    #     dS = Measure('dS', subdomain_data = meshdata.boundaries)
    #     perimeter = assemble(ones*dS(5) + ones*dS(6))


    # berechne Indizes der Vertizes, Facetten und Koordinaten der aeusseren Normalen auf Interface
    interface_indizes = order_boundary_verts(meshdata)
    interface_facet_inds = order_boundary_facets(meshdata, interface_indizes)
    normal_values = get_vertex_normal(meshdata, interface_facet_inds)



    V_vec  = VectorFunctionSpace(meshdata.mesh, 'P', 1, dim = 2)
    signal = Function(V_vec)
    vtod_vec = dolfin.vertex_to_dof_map(V_vec)

    # Koordinaten der aeusseren Normalen fuer FEniCS Funktion
    func_vals = np.zeros(len(signal.vector().get_local()))

    # Berechnung der Koeffizienten des Signalvektorfeldes in Richtung auesserer Normalen
    # hier wird ausgenutzt, dass die Gitterpunkte alle aequidistant sind!

    ###########################
    # korrekter waere: bei allgemeinem Gitter berechne Abstand von einem Punkt zum benachbarten,
    # starte bei 0; addiere naechsten Abstand (genauer: Weg auf Interface --> Randintegral, eig. egal, da Punkte
    # benachbart!)
    # und werte in der mit der Parametrisierung des
    # Interfaces transformierten Signalfunktion Koeffizient aus und weise auf Punkt zu
    ###########################

    vals = (1./len(interface_indizes))*np.array(list(range(0,len(interface_indizes))))
    signal_vals = ampl*np.sin(2.*pi*omega*vals)

    for i in range(0, len(interface_indizes)):
        func_vals[vtod_vec[2*interface_indizes[i]  ]] = eps*normal_values[i]
        func_vals[vtod_vec[2*interface_indizes[i]+1]] = eps*normal_values[i+1]

    print(signal_vals)
    signal.vector()[:] = func_vals



    return signal