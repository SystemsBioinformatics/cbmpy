"""
CBMPy: fluxmodules sparse rationals module
=====================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Arne C. Reimers
Contact email: arne.c.reimers@gmail.com
"""


from fractions import Fraction
#from .myfraction import Fraction
import numpy


class Matrix:
    """ sparse matrix of rational numbers (fractions)

    The matrix grows by demand.
    """


    def __init__(self, listRep=None):
        """" Creates empty matrix.
        """
        self.__cols = []
        self.__nrows = 0

        if listRep != None:
            self.addListRep(listRep)

    def rows(self):
        return self.__nrows

    def cols(self):
        return len(self.__cols)

    def getCol(self, idx):
        """ returns the column at index idx as a frozen set of key-value pairs

        The returned dictionary does not contain any superfluous zeros.
        """
        out = set()
        # sanitize column idx
        if len(self.__cols) <= idx:
            return frozenset()
        else:
            c = self.__cols[idx]
            for k,v in c.items():
                if v == 0:
                    del c[k]
                else:
                    out.add((k,v))
            return frozenset(out)

    def __getitemIdx(self, source, rowidx, out):
        for c in source:
            col = {}
            for r,v in c.items():
                i = rowidx(r)
                if i != None:
                    col[i] = v
            out.__cols.append(col)

    def transpose(self, scalar=1):
        """ creates a transpose of the matrix

        Optionally, the values can be multiplied by a scalar.
        This is for example used to compute the dual matroid.
        """
        assert self.checkRows()
        out = Matrix()
        for i in range(self.rows()):
            out.__cols.append({})
        for ci, c in enumerate(self.__cols):
            for r,v in c.items():
                out.__cols[r][ci] = v*scalar
#        out.__cols[r] = [self.__cols[c][r] for c in range(self.cols()) if self.__cols[c].has_key(r)]
#        [out.__cols[r][ci]=v for  ci, c in enumerate(self.__cols) for r,v in c.items()]
        out.__nrows = len(self.__cols)
        assert out.checkRows()
        return out

    def copy(self):
        """ copies this matrix """
        assert self.checkRows()
        out = Matrix()
        for c in self.__cols:
            c2 = c.copy()
            out.__cols.append(c2)
        out.__nrows = self.__nrows
        assert out.checkRows()
        return out

    def __getitem__(self, indices):
        """ fetch entries of this matrix

        indices must be a 2-tuple indexing rows and columns.
        This returns a matrix.

        If the index is a 2-tuple of ints, a Fraction is returned.
        Otherwise, a Matrix is returned.
        """
        assert len(indices) == 2
        assert self.checkRows()
        rows = indices[0]
        cols = indices[1]

        if isinstance(cols, int):
            if isinstance(rows, int):
                if cols < self.cols():
                    return self.__cols[cols].get(rows,0)
                else:
                    return 0
            else:
                cols = [cols]
        elif isinstance(rows, int):
            rows = [rows]

        if isinstance(rows, slice):
            rlen = len(xrange(*rows.indices(self.rows())))
        else:
            rlen = len(rows)

        out = Matrix()
        out.__nrows = rlen

        if isinstance(rows, slice):
            ind = rows.indices(self.rows())
            idx = lambda x: revslice(x, ind)
        else:
            ind = {}
            for i,v in enumerate(rows):
                ind[v] = i
            idx = lambda x: ind.get(x, None)

        if isinstance(cols, slice):
            subcols = self.__cols[cols]
        else:
            subcols = [self.__cols[i] for i in cols if i < self.cols()]

        self.__getitemIdx(subcols, idx, out)
        assert out.checkRows()
        return out

    def __resetSlice(self, source, rslice):
        assert self.checkRows()
        for c in source:
            for r in c.keys():
                if (r >= rslice.start
                        and r < rslice.stop
                        and r - rslice.start % rslice.step == 0):
                    del c[r]
        assert self.checkRows()

    def __resetTest(self, source, rTest):
        assert self.checkRows()
        for c in source:
            for r in c.keys():
                if r in rTest:
                    del c[r]
        assert self.checkRows()

    def reset(self, rows, cols):
        """ reset all elements in rows x __cols to zero

        This is equivalent (but a bit faster) to self[rows, __cols] = 0
        """
        assert self.checkRows()
        if isinstance(cols, slice):
            subcols = self.__cols[cols]
        else:
            subcols = [self.__cols[i] for i in cols if i < self.cols()]

        if isinstance(rows, slice):
            self.__resetSlice(subcols, rows)
        else:
            self.__resetTest(subcols, rows)
        assert self.checkRows()

    def __iadd__(self, other):
        """ += operation """
        assert self.checkRows()
        assert other.checkRows()
        mcol = other.cols()
        while self.cols() < mcol:
            self.__cols.append({})
        # copy new values
        for i,c in enumerate(other.__cols):
            mycol = self.__cols[i]
            for r,v in c.items():
                nv = mycol.get(r,0) + v
                if nv != 0:
                    mycol[r] = nv
                elif mycol.has_key(r):
                    del mycol[r]

        self.__nrows = max(self.__nrows, other.__nrows)
        assert self.checkRows()
        return self

    def __isub__(self, other):
        """ -= operation """
        assert self.checkRows()
        assert other.checkRows()
        mcol = other.cols()
        while self.cols() < mcol:
            self.__cols.append({})
        # copy new values
        for i,c in enumerate(other.__cols):
            mycol = self.__cols[i]
            for r,v in c.items():
                nv = mycol.get(r,0) - v
                if nv != 0:
                    mycol[r] = nv
                elif mycol.has_key(r):
                    del mycol[r]

        self.__nrows = max(self.__nrows, other.__nrows)
        assert self.checkRows()
        return self

    def __idiv__(self, scalar):
        """ /= operation """
        assert self.checkRows()
        for c in self.__cols:
            for r,v in c.items():
                c[r] = v / scalar
        return self

    def __itruediv__(self, scalar):
        return self.__idiv__(scalar)

    def __imul__(self, v):
        """ *= operation

        works for scalars and matrices.
        """
        assert self.checkRows()
        if isinstance(v, Matrix) and (v.rows() > 1 or v.cols() > 1):
            # proper matrix multiplication
            # we do not care about matching matrix dimensions, because
            # we can always fill up with zeros
            out = Matrix()
            for tr in range(self.rows()):
                for ci,c in enumerate(v.__cols):
                    x = 0
                    for ri, rv in c.items():
                        x += self[tr, ri] * rv
                    out[tr,ci] = x
            # copy result to self
            self.__cols = out.__cols
            # the number of rows stays the same so no update needed
        else:
            # just mutliplication with a scalar
            if isinstance(v, Matrix):
                scalar = v[0,0]
            else:
                scalar = v
            for c in self.__cols:
                for r,val in c.items():
                    c[r] = val * scalar
        return self

    def __add__(self, other):
        """ + operation """
        assert self.checkRows()
        out = Matrix()
        out.__iadd__(self)
        out.__iadd__(other)
        return out

    def __sub__(self, other):
        """ - operation """
        assert self.checkRows()
        out = Matrix()
        out.__iadd__(self)
        out.__isub__(other)
        return out

    def __mul__(self, v):
        """ * operation (for multiplication with scalars)"""
        assert self.checkRows()
        if isinstance(v, Matrix) and (v.rows() > 1 or v.cols() > 1):
            # proper matrix multiplication
            # we do not care about matching matrix dimensions, because
            # we can always fill up with zeros
            out = Matrix()
            for tr in range(self.rows()):
                for ci,c in enumerate(v.__cols):
                    x = 0
                    for ri, rv in c.items():
                        x += self[tr, ri] * rv
                    out[tr,ci] = x
        else:
            # just mutliplication with a scalar
            if isinstance(v, Matrix):
                scalar = v[0,0]
            else:
                scalar = v
            out = Matrix()
            for c in self.__cols:
                oc = {}
                for r,val in c.items():
                    oc[r] = val * scalar
                out.__cols.append(oc)
            out.__nrows = self.__nrows
#        assert self.checkRows()
#        out = Matrix()
#        out.__iadd__(self)
#        out.__imul__(v)
        return out

    def __div__(self, scalar):
        """ / operation (for division by scalars)"""
        assert self.checkRows()
        out = Matrix()
        out.__iadd__(self)
        out.__idiv__(scalar)
        return out

    def coladd(self, col, v, scalar=1):
        """ computes self[:,col] += scalar * v

        col can only be a single index
        """
        assert self.checkRows()
        assert v.checkRows()
        assert v.cols() >= 1

        while self.cols() <= col:
            self.__cols.append({})

        c = self.__cols[col]
        oc = v.__cols[0]
        for r,val in oc.items():
            nv = c.get(r,0) + scalar * val
            if nv != 0:
                c[r] = nv
                if r >= self.__nrows:
                    self.__nrows = r+1
            else:
                del c[r]

        assert self.checkRows()

    def __setitem__(self, indices, v):
        assert len(indices) == 2
        assert self.checkRows()
        rows = indices[0]
        cols = indices[1]

        if isinstance(cols, int):
            cols = [cols]
        if isinstance(rows, int):
            rows = [rows]

        if isinstance(rows, slice):
            if rows.start == None:
                rstart = 0
            else:
                rstart = rows.start
            if rows.step == None:
                maxrow = max(self.rows(), rstart + v.rows())
            elif rows.step > 0:
                maxrow = max(self.rows(), rstart + v.rows()*rows.step)
            else:
                maxrow = self.rows()
            rows = xrange(*rows.indices(maxrow))
        if isinstance(cols, slice):
            if cols.start == None:
                cstart = 0
            else:
                cstart = cols.start
            if cols.step == None:
                maxcol = max(self.cols(), cstart + v.cols())
            elif cols.step > 0:
                maxcol = max(self.cols(), cstart + v.cols()*cols.step)
            else:
                maxcol = self.cols()
            cols = xrange(*cols.indices(maxcol))

        if len(rows) == 1 and len(cols) == 1:
            row = rows[0]
            col = cols[0]
            # just set single value
            if v == 0:
                if self.cols() > col:
                    if self.__cols[col].has_key(row):
                        del self.__cols[col][row]
            else:
                while self.cols() <= col:
                    self.__cols.append({})
                self.__nrows = max(self.__nrows, row+1)
                if isinstance(v, Matrix):
                    assert v.cols() <= 1
                    assert v.rows() <= 1
                    if v.cols() == 0 or v.rows() == 0:
                        if self.__cols[col].has_key(row):
                            del self.__cols[col][row]
                    else:
                        self.__cols[col][row] = v[0,0]
                elif isinstance(v, Fraction):
                    self.__cols[col][row] = v
                else:
                    self.__cols[col][row] = Fraction(v) # make sure its a fraction

        else:
            assert isinstance(v, Matrix)
            assert self.checkRows()
            # empty matrix
            self.reset(rows, cols)
            # add empty columns if necessary
            mcol = max(cols)
            while self.cols() <= mcol:
                self.__cols.append({})
            # copy new values
            for i,c in enumerate(v.__cols):
                for r,val in c.items():
                    assert i < len(cols)
                    assert r < len(rows), "%d < %d %d" % (r, len(rows), v.rows())
                    self.__cols[cols[i]][rows[r]] = val
                    self.__nrows = max(self.__nrows, rows[r]+1)
        assert self.checkRows()

    def __neg__(self):
        """ return -self """
        out = Matrix()
        out.__isub__(self)
        return out

    def __eq__(self, other):
        if isinstance(other, Matrix):
            for i in range(max(len(self.__cols), len(other.__cols))):
                if i >= len(self.__cols):
                    # check if the other column only contains zeros
                    for v in other.__cols[i].values():
                        if v != 0:
                            return False
                elif i >= len(other.__cols):
                    # check if the own column only contains zeros
                    for v in self.__cols[i].values():
                        if v != 0:
                            return False
                else:
                    # is a valid index for both matrices
                    sc = self.__cols[i]
                    oc = other.__cols[i]
                    for k,v in sc.items():
                        if oc[k] != v:
                            return False
                    for k,v in oc.items():
                        if sc[k] != v:
                            return False
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)


    def addMetabolicNetwork(self, cmod):
        """ adds the stoichiometric matrix of a metabolic network

        This is the default method to load the stoichiometric matrix into
        ths matrix. Just create an empty matrix and call this function.

        A list of column names is returned.
        """
        var_spec_id = [s.getId() for s in cmod.species if not s.is_boundary]
        num_col = len(cmod.reactions)
        while self.cols() < num_col:
            self.__cols.append({})
        labels = []
        for c in range(num_col):
            labels.append(cmod.reactions[c].getId())
            for reag in cmod.reactions[c].getStoichiometry():
                if reag[1] in var_spec_id:
                    r = var_spec_id.index(reag[1])
                    self.__cols[c][r] = self.__cols[c].get(r,0) + Fraction(reag[0])
        self.__nrows = len(var_spec_id)
        assert self.checkRows()
        return labels, var_spec_id

    def addListRep(self, matrix):
        """ adds matrix in list representation (ordered by rows) """
        self.__nrows = max(self.__nrows, len(matrix))
        for ri, r in enumerate(matrix):
            while len(r) > self.cols():
                self.__cols.append({})
            for ci, v in enumerate(r):
                if v != 0:
                    col = self.__cols[ci]
                    col[ri] = col.get(ri,0)+Fraction(v)

        assert self.checkRows()

    def toNumpy(self):
        """ returns floating point numpy representation.

        The numpy matrix is a full matrix!
        """
        assert self.checkRows()
        # find max row
        out = numpy.zeros((self.rows(), self.cols()))
        for i,c in enumerate(self.__cols):
            for r,v in c.items():
                out[r,i] = v
        return out

    def checkRows(self):
#        for c in self.__cols:
#            for r in c.keys():
#                assert r < self.rows()
#                if r >= self.rows():
#                    return False
        return True

    def isZero(self):
        for c in self.__cols:
            for r in c.values():
                if r != 0:
                    return False
        return True

def revslice(x, sl):
    if sl[2] == 0:
        return None
    elif sl[2] > 0:
        if x >= sl[1]:
            return None
        if x < sl[0]:
            return None
    else:
        # sl[2] < 0
        if x >= sl[0]:
            return None
        if x < sl[1]:
            return None
    y = x - sl[0]
    if y % sl[2] != 0:
        return None
    return y / sl[2]