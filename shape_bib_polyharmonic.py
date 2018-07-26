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


class MeshData:

    # Objekt mit allen Daten des Gitters
    def __init__(self, mesh, subdomains, boundaries, ind):

        # FEniCS Mesh
        self.mesh = mesh

        # FEniCS Subdomains
        self.subdomains = subdomains

        # FEniCS Boundaries
        self.boundaries = boundaries

        # Indizes der Knotenpunkte mit Traeger nicht am inneren Rand
        self.indNotIntBoundary = ind


class bfgs_memory:
    """
    Klasse, welche alle Information fuer das L-BFGS-Verfahren enthaelt.
    Besteht aus length Gradienten- und Deformationsvektoren, welche jeweils
    in einem Array der Historie nach absteigend sortiert sind.
    Besitzt die Funktion, Gradienten und Deformationen zu updaten, wobei der aelteste
    Eintrag verworfen wird. step_nr ist ein counter fuer das Verfahren.
    """

    # Objekt mit gespeicherten Gradienten und Deformationen der letzten l Schritte in Form von Arrays
    def __init__(self, gradient, deformation, length, step_nr):

        # Liste von Gradientenvektorfeldern
        if (len(gradient) == length): self.gradient = gradient
        else: raise SystemExit("Fehler: Anzahl der Gradienten passt nicht zur Memorylaenge!")

        # Liste von Deformationsvektorfeldern
        if (len(deformation) == length): self.deformation = deformation
        else: raise SystemExit("Fehler: Anzahl der Deformationen passt nicht zur Memorylaenge!")

        # Anzahl der gespeicherten letzten Schritte
        self.length = length

        # Anzahl der bereits ausgefuehrten l-BFGS-Schritte
        if (step_nr >= 0 and isinstance(step_nr, int)): self.step_nr = step_nr
        else: raise SystemExit("Fehler: step_nr muss Integer groesser gleich 0 sein!")

    # macht ein Update der Memory; neueste Elemente bekommen Index 0
    def update_grad(self, upd_grad):

        for i in range(self.length-1): self.gradient[-(i+1)]    = self.gradient[-(i+2)]

        self.gradient[0] = upd_grad

    def update_defo(self, upd_defo):

        for i in range(self.length-1): self.deformation[-(i+1)] = self.deformation[-(i+2)]

        self.deformation[0] = upd_defo

    def initialize_grad(self, meshData, i):
        # erzeugt eine FEniCS-Funktion des Gradientenfeldes i auf einem Mesh aus den gespeicherten Arrays
        # ermoeglicht dadurch Transport; i entspricht Index in Speichermatrix self.gradient, aufsteigend im Alter

        if isinstance(meshData, MeshData): pass
        else: raise SystemExit("initialize_grad benoetigt Objekt der MeshData-Klasse als Input!")

        V = VectorFunctionSpace(meshData.mesh, "P", 1, dim=2)
        f = Function(V)
        f.vector()[:] = self.gradient[i]
        return f

    def initialize_defo(self, meshData, i):
        # erzeugt eine FEniCS-Funktion des Deformationsfeldes i auf einem Mesh aus den gespeicherten Arrays
        # ermoeglicht dadurch Transport; i entspricht Index in Speichermatrix self.deformation, aufsteigend im Alter

        if isinstance(meshData, MeshData): pass
        else: raise SystemExit("initialize_defo benoetigt Objekt der MeshData-Klasse als Input!")

        V = VectorFunctionSpace(meshData.mesh, "P", 1, dim=2)
        f = Function(V)
        f.vector()[:] = self.deformation[i]
        return f


def create_outputfolder():
    '''
    Erstelt einen Outputordner, falls nicht vorhanden
    '''

    try:
        os.mkdir('Output')
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass

    # Erstelle Ordner fuer jeden Durchlauf nach Datum und Zeit
    outputfolder = os.path.join(__this_files_dir,
                                'Output',
                                datetime.now().strftime('%Y%m%d_%H%M%S'))
    os.mkdir(outputfolder)

    # Gibt den Pfad zum speziellen Order im Output Ordner zurueck
    return outputfolder


def load_mesh(name):
    '''
    Initialisiert ein Objekt der MeshData-klasse mittels per GMesh und Dolfin erzeugter
    .xml Dateien
    '''
    # Pfad zur speziellen Gitter Datei
    path_meshFile = os.path.join(DATA_DIR, name)

    # Erstelle FEniCS Mesh mit Subdomains und Boundaries
    mesh = Mesh(path_meshFile + ".xml")
    subdomains = MeshFunction("size_t", mesh,
                              path_meshFile + "_physical_region.xml")
    boundaries = MeshFunction("size_t", mesh,
                              path_meshFile + "_facet_region.xml")

    # Berechne Indizes mit Traeger nicht am inneren Rand
    ind = __get_index_not_interior_boundary(mesh, subdomains, boundaries)

    # Rueckgabe als MeshData Objekt
    return MeshData(mesh, subdomains, boundaries, ind)


def __get_index_not_interior_boundary(mesh, subdomains, boundaries, interior = True):
    """
    Gibt Indizes der Elemente ohne Traeger am inneren Rand zurueck.
    Falls interior = False, so werden die Indizes der Vertices
    des Inneren Randes zurueckgegeben. Diese nicht verwechseln mit DOF-Indizes!
    """

    # Facetten Indizes des inneren Randes bestimmen
    ind_interior_boundary_facets = []
    for i in range(0,len(boundaries)):
        if boundaries[i] > 4:
            ind_interior_boundary_facets.append(i)

    # Knoten Indizes des inneren Randes bestimmen
    ind_interior_boundary_vertices = []
    for c in cells(mesh):
        for f in facets(c):
            if f.index() in ind_interior_boundary_facets:
                for v in vertices(f):
                    ind_interior_boundary_vertices.append(v.index())

    if(interior == False): return list(set(ind_interior_boundary_vertices))

    ind_interior_boundary_vertices = list(set(ind_interior_boundary_vertices))

    # Element Indizes des inneren Randes bestimmen
    ind_around_interior_boundary_cells = []
    for c in cells(mesh):
        ind = False
        for v in vertices(c):
            if v.index() in ind_interior_boundary_vertices:
                ind = True
        if ind:
            ind_around_interior_boundary_cells.append(c.index())

    # Als neue Subdomain definieren
    new_sub = MeshFunction("size_t", mesh, 2)
    new_sub.set_all(0)
    for i in ind_around_interior_boundary_cells:
        if subdomains[i] == 1:
            new_sub[i] = 1
        else:
            new_sub[i] = 2

    # Indizes berechenen mit Traeger nicht am inneren Rand ueber Testproblem
    V = VectorFunctionSpace(mesh, "P", 1, dim=2)

    dx_int = Measure('dx',
                     domain=mesh,
                     subdomain_data=new_sub)

    v = TestFunction(V)

    dummy_y = Constant((1.0, 1.0))

    f_elas_int_1 = inner(dummy_y, v)*dx_int(1)
    F_elas_int_1 = assemble(f_elas_int_1)
    f_elas_int_2 = inner(dummy_y, v)*dx_int(2)
    F_elas_int_2 = assemble(f_elas_int_2)

    # Indizes setzen durch alle Punkte mit Wert 0, die keinen Einfluss haben
    ind1 = (F_elas_int_1.get_local() == 0.0)
    ind2 = (F_elas_int_2.get_local() == 0.0)

    ind = ind1 | ind2

    if(interior): return ind


def mesh_distance(mesh1, mesh2):
    """
    Berechnet einen integrierten Abstand der Minima aller Punkte zweier Formen.
    mesh2 dient dabei als Ausgangsform, d.h. von dort aus wird Abstand gemessen.
    Input sind Objekte der MeshData-Klasse.
    """

    # Berechne Indizes der Vertices der Boundaries
    boundary_index_mesh1 = __get_index_not_interior_boundary(mesh1.mesh, mesh1.subdomains, mesh1.boundaries, interior = False)
    boundary_index_mesh2 = __get_index_not_interior_boundary(mesh2.mesh, mesh2.subdomains, mesh2.boundaries, interior = False)

    # Berechne die Abstaende alle Vertices von mesh1 zu jeweils einem festen Vertex aus mesh2,
    # dann bilde das Minimum und fuege in Liste hinzu
    distance_list_vertex = np.zeros(mesh2.mesh.num_vertices())

    for i in boundary_index_mesh2:
        temp_distance_list = []

        for j in boundary_index_mesh1:
            local_dist = np.linalg.norm(mesh2.mesh.coordinates()[i] - mesh1.mesh.coordinates()[j])
            temp_distance_list.append(local_dist)


        dist = np.amin(temp_distance_list)
        distance_list_vertex[i] = dist

    # definiere eine Funktion auf mesh2 mit Abstandswerten auf Boundaries der Form
    V = FunctionSpace(mesh2.mesh, "P", 1)
    distance_function = Function(V)
    vtod = dolfin.vertex_to_dof_map(V)

    # uebersetze Vertexindizes in zugehoerige DOF-Indizes (der fenicsfunction), da Indizes nicht gleich
    distance_list_dof = np.zeros(len(distance_function.vector().get_local()))
    for i in boundary_index_mesh2: distance_list_dof[vtod[i]] = distance_list_vertex[i]

    # definiere Funktion auf der Form
    distance_function.vector()[:] = distance_list_dof

    # Berechne das zugehoerige Integral
    dS = Measure('dS', subdomain_data=mesh2.boundaries)
    distance_integral = distance_function('+')*dS(5) + distance_function('+')*dS(6)
    value = assemble(distance_integral)

    return value


def targetfunction(meshData, deformation, y_z, fValues, nu):
    """
    Berechnet den Wert des Zielfunktionals nach Verschiebung
    """

    # erzeuge lokale Gitterkopie
    msh = Mesh(meshData.mesh)
    sbd = MeshFunction("size_t", msh, 2)
    sbd.set_values(meshData.subdomains.array())
    bnd = MeshFunction("size_t", msh, 1)
    bnd.set_values(meshData.boundaries.array())
    ind = __get_index_not_interior_boundary(msh, sbd, bnd)
    local_mesh = MeshData(msh, sbd, bnd, ind)

    # verschiebe Gitterkopie
    ALE.move(local_mesh.mesh, deformation)

    # Berechne Zustand in verschobenem Gitter
    y = solve_state(local_mesh, fValues)

    # Assembliere Zielfunktional
    V = FunctionSpace(local_mesh.mesh, 'P', 1)
    z_1 = project(y_z[0], V)
    z_2 = project(y_z[1], V)
    z   = [z_1, z_2]

    j = 1./2.*norm(project(y[1]-z[1], V), 'L2', local_mesh.mesh)**2

    ds = Measure('dS', subdomain_data=local_mesh.boundaries)
    ones = Function(V)
    ones.vector()[:] = 1.
    j_reg_integral = ones*ds(5) + ones*ds(6)
    j_reg = nu*assemble(j_reg_integral)

    J = j + j_reg

    return J


def shape_deriv(meshData, p, y, z, fValues, nu, V):
    """
    Berechnet die Formableitung in Richtung V; V ist FEniCS Funktion, z ist Projektion von y_z
    """

    dx = Measure('dx',
                 domain=meshData.mesh,
                 subdomain_data=meshData.subdomains)
    dS = Measure('dS', subdomain_data=meshData.boundaries)

    f1 = Constant(fValues[0])
    f2 = Constant(fValues[1])
    epsilon_V = sym(nabla_grad(V))
    n = FacetNormal(meshData.mesh)

    # Gradient f ist 0, da f stueckweise konstant
    Dj = ((-inner(grad(y[0]), dot(epsilon_V*2, grad(p[0])))
        -inner(grad(y[1]), dot(epsilon_V * 2, grad(p[1])))) * dx
        + nabla_div(V) * (1 / 2 * (y[1] - z[1]) ** 2 + inner(grad(y[0]), grad(p[0]))
        + inner(grad(y[1]), grad(p[1]))) * dx
        - nabla_div(V)*(f1*p[0] + y[0]*p[1])*dx(1) - nabla_div(V)*(f2*p[0] + y[0]*p[1])*dx(2))

    Dj_reg = nu*((nabla_div(V('+'))
                  - inner(dot(nabla_grad(V('+')), n('+')), n('+')))*dS(5)
                  + (nabla_div(V('+'))
                     - inner(dot(nabla_grad(V('+')), n('+')), n('+')))*dS(6))

    deriv = assemble(Dj + Dj_reg)

    return deriv


def solve_state(meshData, fValues):
    """
    Loest Zustandsgleichung ohne Variationsungleichung
    """

    # Funktionen Raum
    V = FunctionSpace(meshData.mesh, "P", 1)

    # Randbedingungen
    y_out = Constant(0.0)
    bcs = [ DirichletBC(V, y_out, meshData.boundaries, i) for i in range(1, 5) ]

    # Problem definieren
    dx = Measure('dx',
                 domain=meshData.mesh,
                 subdomain_data=meshData.subdomains)

    f1 = Constant(fValues[0])
    f2 = Constant(fValues[1])

    # Loesung der ersten Gleichung berechnen
    y_trial = TrialFunction(V)
    v = TestFunction(V)

    a = lhs(inner(grad(y_trial), grad(v))*dx('everywhere'))
    b = rhs(f1*v*dx(1) + f2*v*dx(2))

    y_1 = Function(V, name="state_sol")
    solve(a == b, y_1, bcs)

    # Loesung der zweiten Gleichung berechnen
    y_trial = TrialFunction(V)
    v = TestFunction(V)

    a = lhs(inner(grad(y_trial), grad(v))*dx('everywhere'))
    b = rhs(y_1*v*dx('everywhere'))

    y_2 = Function(V, name="state_sol")
    solve(a == b, y_2, bcs)

    # Rueckgabe der Loesungen
    y = [y_1, y_2]
    #y.rename("state_sol", "label")

    return y


def solve_adjoint(meshData, y, z):
    """
    Loest Adjungierte Gleichung ohne Variationsungleichung
    """

    # Funktionen Raum
    V = FunctionSpace(meshData.mesh, "P", 1)

    # Randbedingungen
    p_out = Constant(0.0)
    bcs = [ DirichletBC(V, p_out, meshData.boundaries, i) for i in range(1, 5) ]

    # Variationsproblem definieren
    dx = Measure('dx',
                 domain=meshData.mesh,
                 subdomain_data=meshData.subdomains)

    # Loesung der ersten Gleichung berechnen
    p_trial = TrialFunction(V)
    v = TestFunction(V)

    a = lhs(inner(grad(p_trial), grad(v))*dx)
    l = rhs(-y[1]*v*dx + z[1]*v*dx)

    p_1 = Function(V)
    solve(a == l, p_1, bcs)

    # Loesung der zweiten Gleichung berechnen
    p_trial = TrialFunction(V)
    v = TestFunction(V)

    a = lhs(inner(grad(p_trial), grad(v))*dx)
    l = rhs(p_1*v*dx)

    p_2 = Function(V)
    solve(a == l, p_2, bcs)

    # Rueckgabe der Loesung
    # ACHTUNG: ERSTER EINTRAG VON P ENTSPRICHT ZWEITEM EINTRAG VON Y UND VICE VERSA (muss noch geaendert werden)
    p = [p_2, p_1]
    return p


def calc_lame_par(meshData, mu_min_value, mu_max_value):
    """
    Berechnet die lokal variierenden Lame-Parameter mu_elas
    """

    # Funktionen Raum
    V = FunctionSpace(meshData.mesh, "P", 1)

    # Randbedingungen
    mu_min = Constant(mu_min_value)
    mu_max = Constant(mu_max_value)
    bcs = ([ DirichletBC(V, mu_min, meshData.boundaries, i) for i in range(1, 5) ]
           + [ DirichletBC(V, mu_max, meshData.boundaries, i) for i in range(5, 7) ])

    # Variationsproblem definieren
    dx = Measure('dx',
                 domain=meshData.mesh,
                 subdomain_data=meshData.subdomains)

    mu_elas = TrialFunction(V)
    v = TestFunction(V)

    f = Expression("0.0", degree=1)

    a = inner(grad(mu_elas), grad(v))*dx
    l = f*v*dx

    # Loesung berechnen
    mu_elas = Function(V, name="lame_par")
    solve(a == l, mu_elas, bcs)

    # Rueckgabe des Lame-Parameters
    return mu_elas


def solve_linelas(meshData, p, y, z, fValues, mu_elas, nu, zeroed = True):
    """
    Loest lineare Elastizitaetsgleichung ohne Variationsungleichung
    """

    # Funktionenraum
    V = VectorFunctionSpace(meshData.mesh, "P", 1, dim=2)

    # Randbedingungen
    u_out = Constant((0.0, 0.0))
    bcs = [ DirichletBC(V, u_out, meshData.boundaries, i) for i in range (1, 5) ]


    U = TrialFunction(V)
    v = TestFunction(V)

    LHS = bilin_a(meshData, U, v, mu_elas)

    F_elas = shape_deriv(meshData, p, y, z, fValues, nu, v)

    # Alle die keine Traeger am innerend Rand haben auf 0 setzen
    if(zeroed): F_elas[meshData.indNotIntBoundary] = 0.0

    for bc in bcs:
        bc.apply(LHS)
        bc.apply(F_elas)

    # Norm der assemblierten rechten Seite wegen Abbruchbedingung
    nrm_f_elas = norm(F_elas, 'L2', meshData.mesh)

    # Berechne Loesung
    U = Function(V, name="deformation_vec")
    solve(LHS, U.vector(), F_elas)

    # Rueckgabe des Gradientenvektorfeldes U und der Abbruchbedingung
    return U, nrm_f_elas


def bilin_a(meshData, U, V, mu_elas):
    """
    Berechnet den Wert der Bilinearform (lin. El.) fuer gegebene Vektorfelder U, V
    Beide Vektorfelder muessen auf dem selben Mesh definiert sein
    Lame parameter lambda = 0
    """

    dx = Measure("dx", domain=meshData.mesh, subdomain_data=meshData.subdomains)

    epsilon_V = (1./2.)*sym(nabla_grad(V))
    sigma_U   = mu_elas*sym(nabla_grad(U))

    a     = inner(sigma_U, epsilon_V)*dx('everywhere')
    value = assemble(a)

    return value


def bfgs_step(meshData, memory, mu_elas, q_target):
    """
    berechnet aus einer BFGS-memory eine Mesh-Deformation q mittels double-loop-L-BFGS-Verfahren, welche zu memory.grad[0] gehoert
    benoetigt memory.grad[0] als aktuellen Gradienten, memory.deformation[0] als aktuell neueste Deformation
    Output q ist eine Fenics-Funktion der Art Function(V), V=VectorFunctionSpace(mesh, "P", 1, dim=2)
    """

    if isinstance(meshData, MeshData): pass
    else: raise SystemExit("bfgs_step benoetigt Objekt der MeshData-Klasse als Input!")

    if isinstance(memory, bfgs_memory): pass
    else: raise SystemExit("bfgs_step benoetigt Objekt der  BFGS-Memory-Klasse als Input!")

    V = VectorFunctionSpace(meshData.mesh, "P", 1, dim=2)
    q = Function(V)
    q.vector()[:] = q_target
    alpha = np.zeros(memory.length)

    diff_grad = Function(V)
    first_diff_grad = Function(V)


    if(memory.step_nr + 1 >= memory.length):
        # bei voll besetzter memory werden alle Eintraege verwendet

        for i in range(memory.length-1):
            # Vorwaertsschleife
            i         = i+1
            diff_grad.vector()[:] = (memory.initialize_grad(meshData, i-1).vector() - memory.initialize_grad(meshData, i).vector())
            alpha[i]              = bilin_a(meshData, memory.initialize_defo(meshData, i-1), q, mu_elas) / bilin_a(meshData, diff_grad, memory.initialize_defo(meshData, i-1), mu_elas)
            q.vector()[:]         = q.vector() - float(alpha[i])*diff_grad.vector()

        # Reskalierung von q
        first_diff_grad.vector()[:] = (memory.initialize_grad(meshData, 0).vector() - memory.initialize_grad(meshData, 1).vector())
        gamma                       = bilin_a(meshData, first_diff_grad, memory.initialize_defo(meshData, 0), mu_elas) / bilin_a(meshData, first_diff_grad, first_diff_grad, mu_elas)
        q.vector()[:]               = gamma*q.vector()

        for i in range(memory.length-1):
            # Rueckwaertsschleife
            i                     = i+1
            diff_grad.vector()[:] = (memory.initialize_grad(meshData, -(i+1)).vector() - memory.initialize_grad(meshData, -i).vector())
            beta                  = bilin_a(meshData, diff_grad, q, mu_elas) / bilin_a(meshData, diff_grad, memory.initialize_defo(meshData, -(i+1)), mu_elas)
            q.vector()[:]         = q.vector() + (float(alpha[-i]) - beta)*memory.initialize_defo(meshData, -(i+1)).vector()

    elif(memory.step_nr == 0):
            # der erste BFGS-Schritt ist ein Gradientenschritt, U Gradient ist in negativer Richtung
            q.vector()[:] = -1.*q.vector().get_local()
            return q

    else:
        # bei nicht voll besetzter memory werden lediglich die besetzten Eintraege verwendet

        for i in range(memory.step_nr):
            # Vorwaertsschleife
            i         = i+1
            diff_grad.vector()[:] = (memory.initialize_grad(meshData, i-1).vector() - memory.initialize_grad(meshData, i).vector())
            alpha[i]              = bilin_a(meshData, memory.initialize_defo(meshData, i-1), q, mu_elas) / bilin_a(meshData, diff_grad, memory.initialize_defo(meshData, i-1), mu_elas)
            q.vector()[:]         = q.vector() - float(alpha[i])*diff_grad.vector()

        # Reskalierung von q
        first_diff_grad.vector()[:] = (memory.initialize_grad(meshData, 0).vector() - memory.initialize_grad(meshData, 1).vector())
        gamma                       = bilin_a(meshData, first_diff_grad, memory.initialize_defo(meshData, 0), mu_elas) / bilin_a(meshData, first_diff_grad, first_diff_grad, mu_elas)
        q.vector()[:]               = gamma*q.vector()

        for i in range(memory.step_nr):
            # Rueckwaertsschleife
            shift = (memory.length-1) - memory.step_nr
            i                     = i+1
            diff_grad.vector()[:] = (memory.initialize_grad(meshData, -(i+1)-shift).vector() - memory.initialize_grad(meshData, -i-shift).vector())
            beta                  = bilin_a(meshData, diff_grad, q, mu_elas) / bilin_a(meshData, diff_grad, memory.initialize_defo(meshData, -(i+1)-shift), mu_elas)
            q.vector()[:]         = q.vector() + (float(alpha[-i-shift]) - beta)*memory.initialize_defo(meshData, -(i+1)-shift).vector()

    q.vector()[:] = -1. * q.vector()
    return q

