cdef extern from "Cellule.hpp":
    cdef cppclass Cellule:
        Cellule()

        double score

        void evolution()
        void optievolution()

        Cellule *copycellule()

