###############################################
### A very basic factory
class factory:
    def __init__(self):
        self.keys = {}

    def register(self, key, creator):
        self.keys[key] = creator

    def create(self, key, **kwargs):
        creator = self.keys.get(key)
        if not creator: raise ValueError(key)
        return creator(**kwargs)
