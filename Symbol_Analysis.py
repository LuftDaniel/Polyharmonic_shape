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


def order_boundary_verts(meshdata):

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

    print(ordered_interface_indizes[0])
    i = 0
    while(len(ordered_interface_indizes) <= len(interface_indizes)):
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
            # im ersten Knoten hat man 2 Nachbarn moeglich

        else:
            # ein Nachbar ist bereits in der Liste, waehle den anderen
            if local_connectivity[0] in ordered_interface_indizes:
                ordered_interface_indizes.append(local_connectivity[1])

            else: ordered_interface_indizes.append(local_connectivity[0])

        i += 1

    return ordered_interface_indizes


def input_signal(meshdata, omega, ampl):
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

    dS = Measure('dS', subdomain_data = meshdata.boundaries)
    # kann man so umbauen, dass man auch Donut und Hufeisen integrieren kann

    perimeter = assemble(ones*dS(5) + ones*dS(6))

    interface_indizes = order_boundary_verts(meshdata)

    testfunc = Function(V)
    vals = list(range(0, len(interface_indizes)))

    func_vals = np.zeros(len(testfunc.vector().get_local()))
    vtod_scalar = dolfin.vertex_to_dof_map(V)

    for i in range(0, len(interface_indizes)):
        func_vals[vtod_scalar[interface_indizes[i]]] = float(vals[i])

    testfunc.vector()[:] = func_vals


    V_vec  = VectorFunctionSpace(meshdata.mesh, 'P', 1, dim = 2)
    signal = Function(V_vec)
    n      = FacetNormal(meshdata.mesh)


    return testfunc