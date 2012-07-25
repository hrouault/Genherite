cimport Cellule

cdef class PyCellule:
    cdef Cellule.Cellule *pcell
    def __cinit__(self):
        self.pcell = new Cellule.Cellule()
    def __dealloc__(self):
        del self.pcell
    def evolution(self):
        return self.pcell.evolution()
