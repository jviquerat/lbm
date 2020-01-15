# Generic imports
import os
import numpy as np

# Imports with probable installation required
try:
    import meshio
except ImportError:
    print('*** Missing required packages, I will install them for you ***')
    os.system('pip3 install meshio')
    import meshio

### ************************************************
### Mesh class
class Mesh():
    def __init__(self, name):
        # Basic values
        self.name      = name
        self.extension = ''
        self.n_nodes   = 0
        self.n_elts    = 0
        self.n_edges   = 0
        self.n_tris    = 0
        self.n_tets    = 0
        self.node_dim  = 0
        self.elt_dim   = 0

        # Arrays
        self.nodes   = np.array([])
        self.edges   = np.array([])
        self.tris    = np.array([])
        self.tets    = np.array([])

### ************************************************
### Function to read mesh in cemef format
def read_cemef_mesh(mesh):
    # Check that file exists
    if (not os.path.isfile(mesh.name+'.t')):
        print('Input file does not exist')
        quit()

    # Open file and handle header line
    with open(mesh.name+'.t') as file:
        header          = file.readline().split()
        mesh.n_nodes    = int(header[0])
        mesh.n_elts     = int(header[2])

        # Handle dimensionality
        mesh.node_dim    = int(header[1])
        mesh.elt_dim     = int(header[3])

        # Allocate arrays to max size
        mesh.nodes      = np.zeros([mesh.n_nodes,3])
        mesh.edges      = np.zeros([mesh.n_elts, 2], dtype=int)
        mesh.tris       = np.zeros([mesh.n_elts, 3], dtype=int)
        mesh.tets       = np.zeros([mesh.n_elts, 4], dtype=int)

        # Read vertices
        for i in range(0,mesh.n_nodes):
            line      = file.readline().split()
            dim_nodes = len(line)
            for j in range(0,dim_nodes):
                mesh.nodes[i,j] = float(line[j])
                for j in range(dim_nodes,3):
                    mesh.nodes[i,j] = 0.0

        # Possible empty line
        pos = file.tell()
        line = file.readline().split()
        if (len(line) != 0):
            file.seek(pos)

        # Read elements
        for i in range(0,mesh.n_elts):
            line = file.readline().split()

            # Tetrahedra
            if ((mesh.elt_dim == 4) and (int(line[3]) != 0)):
                mesh.tets[mesh.n_tets,0] = int(line[0]) - 1
                mesh.tets[mesh.n_tets,1] = int(line[1]) - 1
                mesh.tets[mesh.n_tets,2] = int(line[2]) - 1
                mesh.tets[mesh.n_tets,3] = int(line[3]) - 1
                mesh.n_tets             += 1

            # Triangles
            if (   ((mesh.elt_dim == 3) and (int(line[2]) != 0))
                or ((mesh.elt_dim == 4) and (int(line[3]) == 0))):
                mesh.tris[mesh.n_tris,0] = int(line[0]) - 1
                mesh.tris[mesh.n_tris,1] = int(line[1]) - 1
                mesh.tris[mesh.n_tris,2] = int(line[2]) - 1
                mesh.n_tris             += 1

            # Edges
            if (int(line[2]) == 0):
                mesh.edges[mesh.n_edges,0] = int(line[0]) - 1
                mesh.edges[mesh.n_edges,1] = int(line[1]) - 1
                mesh.n_edges              += 1

        # Resize arrays
        mesh.edges = mesh.edges[0:mesh.n_edges,:]
        mesh.tris  = mesh.tris [0:mesh.n_tris, :]
        mesh.tets  = mesh.tets [0:mesh.n_tets, :]

### ************************************************
### Function to write mesh in cemef format
def write_cemef_mesh(mesh):
    # Handle dimensionality and format-specific features
    if ((mesh.node_dim == 0) or (mesh.elt_dim == 0)):
        if (mesh.n_tets == 0):
            mesh.node_dim = 2
            mesh.elt_dim  = 3
        else:
            mesh.node_dim = 3
            mesh.elt_dim  = 4
            mesh.n_edges  = 0
            mesh.edges    = np.array([])
    mesh.n_elts = mesh.n_edges + mesh.n_tris + mesh.n_tets

    # Open file
    with open(mesh.name+'.t','w') as file:

        # Handle header line
        file.write('{} {} {} {}\n'.format(
            mesh.n_nodes, mesh.node_dim, mesh.n_elts, mesh.elt_dim))

        # Write vertices
        for i in range(0,mesh.n_nodes):
            for j in range(0,mesh.node_dim):
                file.write('{} '.format(mesh.nodes[i,j]))
            file.write('\n')

        # Write tetrahedra
        for i in range(0,mesh.n_tets):
            file.write('{} {} {} {}'.format(mesh.tets[i,0] + 1,
                                            mesh.tets[i,1] + 1,
                                            mesh.tets[i,2] + 1,
                                            mesh.tets[i,3] + 1))
            file.write('\n')

        # Write triangles
        for i in range(0,mesh.n_tris):
            file.write('{} {} {} '.format(mesh.tris[i,0] + 1,
                                          mesh.tris[i,1] + 1,
                                          mesh.tris[i,2] + 1))
            for j in range(3,mesh.elt_dim):
                file.write('{} '.format(0))
            file.write('\n')

        # Write edges
        for i in range(0,mesh.n_edges):
            file.write('{} {} '.format(mesh.edges[i,0] + 1,
                                       mesh.edges[i,1] + 1))
            for j in range(2,mesh.elt_dim):
                file.write('{} '.format(0))
            file.write('\n')

### ************************************************
### Function to read mesh with meshio
def read_meshio_mesh(mesh):
    mesh_copy    = meshio.read(mesh.name+'.'+mesh.extension)
    mesh.nodes   = mesh_copy.points
    mesh.n_nodes = len(mesh.nodes)

    if ('line'     in mesh_copy.cells):
        mesh.edges   = mesh_copy.cells.get('line')
        mesh.n_edges = len(mesh.edges)

    if ('triangle' in mesh_copy.cells):
        mesh.tris    = mesh_copy.cells.get('triangle')
        mesh.n_tris  = len(mesh.tris)

    if ('tetra'    in mesh_copy.cells):
        mesh.tets    = mesh_copy.cells.get('tetra')
        mesh.n_tets  = len(mesh.tets)

### ************************************************
### Function to write mesh with meshio
def write_meshio_mesh(mesh):
    nodes = mesh.nodes
    cells = {}
    if ((mesh.n_edges > 0) and (mesh.extension != 'xml')):
        cells['line']     = mesh.edges
    if (mesh.n_tris  > 0):
        cells['triangle'] = mesh.tris
    if (mesh.n_tets  > 0):
        cells['tetra']    = mesh.tets
    mesh_copy = meshio.Mesh(nodes, cells)
    meshio.write(mesh.name+'.'+mesh.extension, mesh_copy)

### ************************************************
### Read mesh file in any format
def read_mesh(filename, extension=''):
    # Handle extension
    if (extension == ''):
        name, extension = filename.split('.')[:]
        if (extension == ''):
            print('I could not determine mesh extension')
            print('Please provide a valid mesh file with extension')
            quit()
    if (not extension in ('t','mesh','xml','xdmf')):
        print('The mesh format you provided is unknwon')
        quit()

    # Read mesh
    mesh           = Mesh(name)
    mesh.extension = extension

    if (extension == 't'):
        read_cemef_mesh(mesh)
    if (extension in ('mesh','xml','xdmf')):
        read_meshio_mesh(mesh)

    return mesh

### ************************************************
### Write mesh in any format
def write_mesh(mesh, extension=''):
    # Handle the case of an input like '.extension'
    trash, extension = extension.split('.')[:]

    # Handle extension
    if (extension == ''):
        if (mesh.extension == ''):
            print('I could not determine which mesh extension to use')
            print('By default I will use medit format')
            extension='mesh'
        extension = mesh.extension
    if (not extension in ('t','mesh','xml','xdmf')):
        print('The mesh format you provided is unknwon')
        print('By default I will use medit format')
        extension='mesh'

    # Reset mesh extension to output extension
    mesh.extension = extension

    # Write mesh
    if (extension == 't'):
        write_cemef_mesh(mesh)
    if (extension in ('mesh','xml','xdmf')):
        write_meshio_mesh(mesh)
