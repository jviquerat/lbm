# Generic imports
import PIL
import copy
import numpy             as np
import matplotlib        as mplt
import matplotlib.pyplot as plt

### ************************************************
### Class defining lattice object
class Lattice:
    ### ************************************************
    ### Constructor
    def __init__(self,
                 name = None,
                 xmin = None,
                 xmax = None,
                 ymin = None,
                 ymax = None,
                 nx   = None,
                 ny   = None):

        if (name is None): name = 'lattice'
        if (xmin is None): xmin =-10.0
        if (xmax is None): xmax = 20.0
        if (ymin is None): ymin = -5.0
        if (ymax is None): ymax = -5.0
        if (nx   is None): nx   = 100
        if (ny   is None): ny   = 200

        self.name = name
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.nx   = nx
        self.ny   = ny

    ### ************************************************
    ### Generate lattice
    def generate(self, polygon):

        # Because we loop on the lattice left-right and top-down,
        # we need to flip the polygon up-down
        poly       = copy.deepcopy(polygon)
        poly[:,1] *= -1.0

        # Declare lattice array
        self.lattice = np.zeros((self.ny, self.nx), dtype=bool)

        # Fill lattice
        for i in range(self.nx):
            for j in range(self.ny):
                pt           = self.lattice_coords(j, i)
                inside       = self.is_inside(poly, pt)
                self.lattice[j,i] = inside

    ### ************************************************
    ### Get lattice coordinates from integers
    def lattice_coords(self, j, i):

        # Compute and return the coordinates of the lattice node (i,j)
        dx = (self.xmax - self.xmin)/(self.nx - 1)
        dy = (self.ymax - self.ymin)/(self.ny - 1)
        x  = self.xmin + i*dx
        y  = self.ymin + j*dy

        return [x, y]

    ### ************************************************
    ### Determine if a pt is inside or outside a closed polygon
    def is_inside(self, poly, pt):

        # Initialize
        j         = len(poly) - 1
        odd_nodes = False

        # Check if point is inside or outside
        # This is a valid algorithm for any non-convex polygon
        for i in range(len(poly)):
            if (((poly[i,1] < pt[1] and poly[j,1] >= pt[1])  or
                 (poly[j,1] < pt[1] and poly[i,1] >= pt[1])) and
                 (poly[i,0] < pt[0] or  poly[j,0]  < pt[0])):

                # Compute slope
                slope = (poly[j,0] - poly[i,0])/(poly[j,1] - poly[i,1])

                # Check side
                if ((poly[i,0] + (pt[1] - poly[i,1])*slope) < pt[0]):
                    odd_nodes = not odd_nodes

            # Increment
            j = i

        return odd_nodes

    ### ************************************************
    ### Generate lattice image
    def generate_image(self):
        plt.axis('off')
        plt.imshow(self.lattice,
                   cmap=mplt.cm.inferno)
        plt.savefig(self.name, dpi=200, bbox_inches='tight')
        plt.close()
        self.trim_white(self.name+'.png')

    ### ************************************************
    ### Crop white background from image
    def trim_white(self, filename):
        im   = PIL.Image.open(filename)
        bg   = PIL.Image.new(im.mode, im.size, (255,255,255))
        diff = PIL.ImageChops.difference(im, bg)
        bbox = diff.getbbox()
        cp   = im.crop(bbox)
        cp.save(filename)
