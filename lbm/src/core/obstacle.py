### ************************************************
### Class defining an obstacle in the lattice
class obstacle:
    def __init__(self, name, n_pts, n_spts, type, size, pos):

        self.name   = name
        self.n_pts  = n_pts
        self.n_spts = n_spts
        self.type   = type
        self.size   = size
        self.pos    = pos

    def set_polygon(self, polygon):

        self.polygon = polygon

    def set_tag(self, tag):

        self.tag = tag

    def fill(self, area, boundary, ibb):

        self.area     = area
        self.boundary = boundary
        self.ibb      = ibb
