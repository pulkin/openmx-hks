# ======================================================================
# This file is under construction.
# ======================================================================

import numpy
from scipy import linalg
from warnings import warn

def greens_function(
    W0,
    W_positive,
    W_negative,
    tolerance,
    maxiter,
):

    es0 = W0
    e00 = W0
    alp = -W_positive
    bet = -W_negative
    gr00 = linalg.inv(es0)
    gt = gr00

    iterations = 0

    while True:
        gr02 = linalg.inv(e00)
        gr01 = numpy.dot(gr02, bet)
        gr00 = numpy.dot(alp, gr01)
        es0 = es0 - gr00
        gr01 = numpy.dot(gr02, alp)
        gr00 = numpy.dot(gr02, bet)
        gr02 = numpy.dot(bet, gr01)
        e00 = e00 - gr02
        gr02 = numpy.dot(alp, gr00)
        e00 = e00 - gr02
        gr02 = numpy.dot(alp, gr01)
        alp = gr02
        gr02 = numpy.dot(bet, gr00)
        bet = gr02
        gr00 = linalg.inv(es0)
        rms = numpy.abs(gt-gr00).max()
        iterations += 1
        if rms>tolerance and iterations<maxiter:
            gt = gr00
        else:
            break

    if rms>tolerance:
        raise Exception("Green's function iteration error: after {:d} iterations the error is {:e} (required {:e})".format(iterations, rms, tolerance))

    return gr00

class LinearOperator(object):
    """
    A class defining basic algebra and commuting rules.
    """

    def zeros_like(self):
        """
        A compaitable "0" operator.
        """
        raise NotImplementedError

    def eye_like(self):
        """
        A compaitable "1" operator.
        """
        raise NotImplementedError

    def tr(self):
        """
        Transpose.
        """
        raise NotImplementedError

    def cc(self):
        """
        Complex conjugate.
        """
        raise NotImplementedError

    def hc(self):
        """
        Hermitian conjugate.
        """
        return self.tr().cc()

    def absmax(self):
        """
        A measure for the deviation from zero.
        """
        raise NotImplementedError

    def __neg__(self):
        raise NotImplementedError

    def __abs__(self):
        raise NotImplementedError

    def __eq__(self, other):
        return (self-other).absmax() == 0

    def __add__(self, other):

        if isinstance(other, (int,float,complex)):
            return self + self.eye_like()*other

        else:
            raise NotImplementedError

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return - self + other

    def __mul__(self, other):
        raise NotImplementedError

    def __rmul__(self, other):
        return self*other

    def __truediv__(self, other):

        if isinstance(other, (int,float,complex)):
            return self*(1./other)

        else:
            raise NotImplementedError

class TightBinding(LinearOperator):

    def __init__(self, m, vectors = None, dimensions = None, shape = None):
        """
        A class representing tight-binding (periodic) matrix.

        Args:

            m (dict or list): a tight-binding matrix. In the case of a
            dict, the dict keys represent integers corresponding to
            matrix block location while dict values are block matrices.
            Otherwise the list items are tight-binding block matrices and
            the keyword argument ``vectors'' providing block indices has
            to be supplied separately.

        Kwargs:

            vectors (list): a list of tight-binding block indices for the
            case of ``m'' being a list;

            dimensions (int): number of dimensions if no blocks
            provided or, equivalently, this TightBinding is zero;

            shape (tuple of ints): shape of the blocks if no blocks
            provided or, equivalently, this TightBinding is zero.
        """

        if len(m) == 0:

            if dimensions is None or shape is None:
                raise ValueError("Empty input: could not determine number of dimensions and/or shape of tight-binding blocks")

            if not isinstance(dimensions,int):
                raise ValueError("Dimensions keyword argument should be integer")

            a,b = shape
            if not isinstance(a,int) or not isinstance(b,int):
                raise ValueError("Shape keyword argument should be a tuple of two integers")

            self.dims = dimensions
            self.shape = (a,b)
            self.__m__ = {}

        else:
            #TODO: Move all this mess to __setitem__

            if isinstance(m,dict):
                if not vectors is None:
                    warn("Vectors keyword argument is ignored")

            else:
                if vectors is None:
                    m = {tuple():m}
                else:
                    m = dict(zip((tuple(i) for i in vectors), m))

            t = None

            for k,v in m.items():

                if isinstance(v, (list,tuple)):
                    tt = numpy.ndarray
                else:
                    tt = type(v)

                if t is None:
                    t = tt

                if not t == tt:
                    raise ValueError("Inconsistent types along the input dict: {} and {}".format(str(t),str(tt)))

            self.__m__ = {}

            if t == numpy.ndarray or t == numpy.core.memmap:

                for k,v in m.items():

                    if isinstance(k,int):
                        k = (k,)

                    self.__m__[k] = numpy.array(v, dtype = numpy.complex)
            elif t == TightBinding:

                for k,v in m.items():

                    if not isinstance(k, int):
                        raise ValueError("The keys in the input should be ints, found {} instead".format(str(type(k))))

                    for kk, vv in v.__m__.items():
                        self.__m__[(k,)+kk] = numpy.array(vv, dtype = numpy.complex)

            else:
                raise ValueError("Unknown type: {}".format(str(t)))

            self.dims = None
            self.shape = None
            d1 = None

            for k,v in self.__m__.items():

                shape = tuple(v.shape)
                if not len(shape) == 2:
                    raise ValueError("{} is not a 2D matrix: shape = {}".format(str(k), str(shape)))

                if self.shape is None:
                    self.dims = len(k)
                    d1 = k
                    self.shape = shape

                elif not self.dims == len(k):
                    raise ValueError("Inconsistent dimensions: {} vs {}".format(str(d1),str(k)))

                elif not self.shape == shape:
                    raise ValueError("Inconsistent matrix size: {} in {} vs {} in {}".format(str(self.shape), str(d1), str(shape), str(k)))

            if not dimensions is None and not dimensions == self.dims:
                raise ValueError("Dimensions keyword argument = {:d} does not correspond to the input with {:d} dimensions".format(dimensions, self.dims))

            if not shape is None and not shape == self.shape:
                raise ValueError("Shape keyword argument = {} does not correspond to the input with shape {}".format(repr(shape), repr(self.shape)))

    # ==================================================================
    # Serialization
    # ==================================================================

    def __getstate__(self):
        return dict(
            dims=self.dims,
            shape=self.shape,
            __m__=self.__m__,
        )

    def __setstate__(self, data):
        self.dims = data["dims"]
        self.shape = data["shape"]
        self.__m__ = data["__m__"]

    # ==================================================================
    # Constructors
    # ==================================================================

    @staticmethod
    def zeros(dims,shape):
        """
        Prepares a zero TightBinding of the given geometry.

        Args:

            dims (int): number of dimensions;

            shape (tuple of 2 ints): shape of matrix blocks;

        Returns:

            A zero TightBinding.
        """
        return TightBinding({}, dimensions = dims, shape = shape)

    @staticmethod
    def eye(dims,shape):
        """
        Prepares an "eye" TightBinding of the given geometry.

        Args:

            dims (int): number of dimensions;

            shape (tuple of 2 ints): shape of matrix blocks;

        Returns:

            An "eye" TightBinding.
        """
        return TightBinding({
            (0,)*dims : numpy.eye(shape[0], M = shape[1]),
        })

    # ==================================================================
    # Helpers
    # ==================================================================

    def foreach(self, p, **kwargs):
        """
        Performs array-wise operation on this TightBinding.

        Args:

            p (function): a function accepting key and value and returning
            a new key and a new value.

        Kwargs:

            All keyword arguments are passed to the constructor.

        Returns:

            A new TightBinding with a given operation applied.
        """
        result = {}
        for k,v in self.__m__.items():
            r = p(k,v)
            if not r is None and not r[0] is None and not r[1] is None:
                result[tuple(r[0])] = r[1]
        return TightBinding(result, **kwargs)

    @staticmethod
    def __tr_i__(i):
        return tuple(-ii for ii in i)

    def test_compatible(self, other):
        """
        Tests whether this and the other TightBinding are compatible to
        perform simple arithmetics. Raises an exception otherwise.
        """
        if self.dims != other.dims:
            raise ValueError("Dimension mismatch: {:d} vs {:d}".format(self.dims, other.dims))

        if self.shape != other.shape:
            raise ValueError("Block shape mismatch: {} vs {}".format(self.shape, other.shape))

    # ==================================================================
    # Required by parent
    # ==================================================================

    def zeros_like(self):
        """
        Prepares a zero TightBinding of the same geometry as self.

        Returns:

            A zero TightBinding.
        """
        return TightBinding.zeros(self.dims, self.shape)

    def eye_like(self):
        """
        Prepares an "eye" TightBinding of the same geometry as self.

        Returns:

            An "eye" TightBinding.
        """
        return TightBinding.eye(self.dims, self.shape)

    def tr(self):
        """
        Calculates transposed TightBinding.

        Returns:

            A transposed of the TightBinding.
        """
        return self.foreach(
            lambda k,v: (TightBinding.__tr_i__(k), numpy.transpose(v)),
            dimensions = self.dims,
            shape = self.shape,
        )

    def cc(self):
        """
        Calculates complex conjugate TightBinding.

        Returns:

            A complex conjugate of the TightBinding.
        """
        return self.foreach(
            lambda k,v: (k, numpy.conj(v)),
            dimensions = self.dims,
            shape = self.shape,
        )

    def absmax(self, location = False):
        """
        Retrieves the maximum value by modulus.

        Kwargs:

            location (bool): returns location of maximum as well

        Returns:

            A float with the maximum value.
        """
        mx = None
        block = None
        index = None
        for k,v in self.__m__.items():
            if mx is None:
                block = k
                index = numpy.unravel_index(numpy.argmax(abs(v)),v.shape)
                mx = abs(v[index])
            else:
                i = numpy.unravel_index(numpy.argmax(abs(v)),v.shape)
                if abs(v[i]) > mx:
                    block = k
                    index = i
                    mx = abs(v[index])

        if location:
            return mx, block, index

        else:
            return mx

    def __neg__(self):
        return self.foreach(lambda k,v: (k,-v))

    def __abs__(self):
        return self.foreach(lambda k,v: (k,abs(v)))

    def __add__(self, other):

        if isinstance(other, TightBinding):

            self.test_compatible(other)

            keys = set(self.__m__.keys()) | set(other.__m__.keys())
            result = {}
            for k in keys:
                result[k] = self.__m__.get(k,0) + other.__m__.get(k,0)
            return TightBinding(result)

        else:
            return super(TightBinding,self).__add__(other)

    def __mul__(self, other):
        return self.foreach(lambda k,v: (k,v*other))

    # ==================================================================
    # Matrix blocks access
    # ==================================================================

    def __getattr__(self, name):

        if name == "diagonal":
            key = (0,)*self.dims
            if key in self.__m__:
                return self[key]
            else:
                return numpy.zeros(self.shape)

        else:
            raise AttributeError("{} object has no attribute {}".format(self.__class__.__name__, repr(name)))

    def __setattr__(self, name, value):

        if name == "diagonal":
            key = (0,)*self.dims
            self[key] = value

        else:
            super(TightBinding,self).__setattr__(name,value)

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)

        if key in self.__m__:
            return self.__m__[key]

        else:
            if not len(key) == self.dims:
                raise ValueError("Argument number mismatch: found {:d}, required {:d}".format(len(key), self.dims))

            return numpy.zeros(self.shape, dtype = numpy.complex)

    def __setitem__(self, key, item):
        if not isinstance(key, tuple):
            key = (key,)

        if not len(key) == self.dims:
            raise ValueError("Argument number mismatch: found {:d}, required {:d}".format(len(key), self.dims))

        for k in key:
            if not isinstance(k, (int, numpy.integer)):
                raise ValueError("The keys should be ints, found {} instead".format(str(type(k))))

        item = numpy.array(item, dtype = numpy.complex)

        if not len(item.shape) == 2:
            raise ValueError("Not a 2D matrix: shape = {}".format(str(item.shape)))

        if not tuple(item.shape) == self.shape:
            raise ValueError("Wrong dimensions: shape = {}".format(str(item.shape)))

        self.__m__[key] = item

    # ==================================================================
    # Other arithmetic operations preserving geometry
    # ==================================================================

    def nan(self, check_inf = True):
        """
        Checks if NaN are present in the tight binding.

        Kwargs:

            check_inf (bool): checks infinities as well.

        Returns:

            A TightBinding with empty matrices where corresponding
            (NaN, inf) elements are set to one.
        """
        return self.foreach(
            lambda k,v: (k, numpy.isnan(v) + (check_inf == True) * numpy.isinf(v)),
            dimensions = self.dims,
            shape = self.shape,
        )

    # ==================================================================
    # Arithmetic operations affecting geometry
    # ==================================================================

    def squeezed(self,*dims):
        """
        Squeezes dimensions from this TightBinding.

        Args:

            dimension(s) to be squeezed or "all".

        Returns:

            A TightBinding with dimensions squeezed out.
        """
        if len(dims) == 1 and dims[0] == "all":

            return self.squeezed(*tuple(i for i in range(self.dims) if self.is_isolated(i)))

        else:

            for i in dims:
                if not self.is_isolated(i):
                    raise ValueError("This TightBinding is not isolated along dimension {:d}".format(i))

            return self.foreach(
                lambda k,v: (tuple(j for i,j in enumerate(k) if not i in dims),v),
                dimensions = self.dims-len(dims),
                shape = self.shape,
            )

    def insert_dimension(self, position):
        """
        Adds a dimension to this TightBinding.

        Args:

            position (int): a new index of the dimension.

        Returns:

            A new TightBinding with the dimension added.
        """
        return self.foreach(
            lambda k,v: (k[:position] + (0,) + k[position:], v),
            dimensions = self.dims+1,
            shape = self.shape,
        )

    def fourier(self, k, index = None):
        """
        Performs a Fourier transform.

        Args:

            k (float): the wave number or a wave vector.

        Kwargs:

            index (int): index to transform; if None transforms along
            first ``len(k)'' indeces.

        Returns:

            Transformed tight binding matrix.
        """
        if self.dims == 1 and index is None:
            index = 0

        if not index is None:

            new_data = {}

            for key, value in self.__m__.items():

                key2 = tuple(key[:index] + key[index+1:])
                matrix = value*numpy.exp(2j*numpy.pi*key[index]*k)

                if key2 in new_data:
                    new_data[key2] += matrix

                else:
                    new_data[key2] = matrix

            return TightBinding(new_data)

        else:

            result = self
            for i in k:
                result = result.fourier(i, index = 0)
            return result

    # ==================================================================
    # Classification
    # ==================================================================

    def is_square(self):
        """
        Checks if it is a square matrix.

        Returns:

            True if it is a square matrix.
        """
        return self.shape[0] == self.shape[1]

    def is_1DNN(self):
        """
        Determines if it is a valid 1D nearest neighbour matrix.

        Returns:

            True if it is a valid one.
        """
        if not self.dims == 1:
            return False
        if not set(self.__m__.keys()) <= set(((0,),(1,),(-1,))):
            return False

        return True

    def is_hermitian(self, eps = 0):
        """
        Determines if Hermitian down to a given precision.

        Kwargs:

            eps (float): epsilon

        Returns:

            True if Hermitian.
        """
        return (self - self.hc()).absmax() <= eps

    def is_trivial(self):
        """
        Determines if this is a trivial (0-dim) TightBinding.

        Returns:

            True if trivial.
        """
        return self.dims == 0

    def is_isolated(self, d = None):
        """
        Determines if this is an isolated TightBinding. If ``d'' is
        provided checks only apply to the specified dimension.

        Kwargs:

            d (int): dimension to check.

        Returns:

            True if is isolated.
        """
        if d is None:
            return len(self.__m__) == 0 or (len(self.__m__) == 1 and (0,)*self.dims in self.__m__)

        else:
            for k in self.__m__.keys():
                if not k[d] == 0:
                    return False

            return True

    # ==================================================================
    # Utility
    # ==================================================================

    def __repr__(self):
        return "TightBinding dim {:d} [{:d}x{:d}]".format(self.dims, self.shape[0], self.shape[1])

    def copy(self):
        """
        Calculates a deep copy.

        Returns:

            A copy.
        """
        return TightBinding(dict((k,v.copy()) for k,v in self.__m__.items()), dimensions = self.dims, shape = self.shape)

    def size(self):
        """
        Retrieves the number of tight-binding matrices.

        Returns:

            The number of tight-binding matrices.
        """
        return len(self.__m__)

    def simplify(self):
        """
        Simplifies this TightBinding in-place by removing all-zero submatrices.
        """
        new = {}
        for k,v in self.__m__.items():
            if numpy.any(v!=0):
                new[k] = v
        self.__m__ = new

    def simplified(self):
        """
        Simplifies this TightBinding by taking out all zero matrices.

        Returns:

            A simplified TightBinding.
        """
        return self.foreach(
            lambda k,v: (k,v) if numpy.count_nonzero(v)>0 else None,
            dimensions = self.dims,
            shape = self.shape,
        )

    def inverted(self, dim = 0):
        """
        Creates a space-inverted version of self.

        Kwargs:

            dim (int): dimension to invert;

        Returns:

            An inverted version of self.
        """
        return self.foreach(
            lambda k,v: (k[:dim] + (-k[dim],) + k[dim+1:],v),
            dimensions = self.dims,
            shape = self.shape,
        )

    def __s2m__(self, s, d):
        """ Formats the input into array of indexes."""
        if isinstance(s, (list, tuple, numpy.ndarray, set)):
            if isinstance(s, set):
                s = sorted(s)
            s = numpy.array(s)
            if s.dtype == numpy.bool or s.dtype == numpy.int:
                return numpy.arange(self.shape[d])[s]
            else:
                raise ValueError("Unknown type: {}".format(s.dtype))
        elif isinstance(s, slice):
            return numpy.arange(self.shape[d])[s]
        else:
            raise ValueError("Cannot format input: {}".format(s))

    def subsystem(self, rows, columns):
        """
        Creates a subsystem of this TightBinding by applying a given
        selection to each matrix in the Hamiltonian.

        Args:

            rows (array): row selection;

            columns (array): columns selection;

        Returns:

            A subsystem TightBinding selected.
        """
        x = self.__s2m__(rows,0)
        y = self.__s2m__(columns,1)
        target = numpy.ix_(x,y)
        return self.foreach(
            lambda k,v: (k, v[target]),
            dimensions = self.dims,
            shape = (len(x), len(y)),
        )

    def eig_path(self, pts, b = None):
        """
        Calculates eigenvalues along a path in k space.

        Args:

            pts (array): an array containing k-points to calculate at;

        Kwargs:

            b (TightBinding): the rhs of a generailzed eigenvalue problem;

        Returns:

            Eigenvalues in a multidimensional array.
        """
        pts = numpy.array(pts)
        result = []
        for i in pts:

            if isinstance(i, (float, numpy.float)):
                if self.dims != 1:
                    raise ValueError("Plain 1D k-point arrays are supported only for 1D TightBindings, this one is {:d}D".format(self.dims))

            elif isinstance(i, numpy.ndarray):
                if len(i) != self.dims:
                    raise ValueError("The number of actual dimensions {:d} is different from the number of k-point dimensions {:d}".format(self.dims, len(i)))

            else:
                raise ValueError("Unknown k-point: {}".format(repr(i)))

            hh = self.fourier(i).diagonal

            if not b is None:
                bb = b.fourier(i).diagonal

            else:
                bb = None

            try:
                result.append(linalg.eigvalsh(hh, b = bb))

            except numpy.linalg.linalg.LinAlgError:

                # Workaround for ill-conditioned overlap
                bb_val, bb_vec = linalg.eig(bb)
                bb_val = numpy.abs(bb_val)
                bb = bb_vec.dot(numpy.diag(bb_val)).dot(bb_vec.conj().T)

                warn("Failed solving generalized eigenvalue equation, using a workaround")
                result.append(linalg.eigvalsh(hh, b = bb))

        return numpy.array(result)

    def shift_subblock(self, rows, columns, vector):
        """
        Shifts TightBinding submatrices along a vector. May be considered
        as moving one or more atoms into a different unit cell.

        Args:

            rows (array): submatrix rows mask;

            columns (array): submatrix columns mask;

            vector (array): an array of integers specifying the submatrix
            shift.
        """
        move = self.subsystem(rows, columns).simplified()

        if move != 0:

            target = numpy.ix_(self.__s2m__(rows,0),self.__s2m__(columns,1))

            # Set moving blocks to zero
            for v in self.__m__.values():
                v[target] = 0
            self.simplify()

            for k,v in move.__m__.items():
                new_k = tuple(i+j for i,j in zip(k,vector))
                m = self[new_k]
                m[target] = v
                self[new_k] = m

    def hermitian_shift_subblock(self, selection, vector):
        """
        Performs a shift of a subblock of a square Hermitian matrix. This
        operation leaves eigenvalues intact while amplitudes in eigenvectors
        acquire a phase.

        Args:

            selection (array): diagonal submatrix mask;

            vector (array): an array of integers specifying the submatrix
            shift.
        """
        if not self.is_square():
            raise ValueError("Not a square matrix")

        vector = numpy.array(vector)
        s = numpy.zeros(self.shape[0], dtype=bool)
        s[self.__s2m__(selection,0)] = True
        self.shift_subblock(s, ~s, vector)
        self.shift_subblock(~s, s, -vector)

    def super(self, size, dim):
        """
        Creates a supercell version of this tight binding.

        Args:

            size (int): a factor by which matrix blocks are increased;

            dim (int): dimension along which a supercell is formed.

        Returns:

            A supercell TightBinding.
        """
        if not 0 <= dim < self.dims:
            raise ValueError("Wrong dimension: {:d} for {:d}-dim TightBinding".format(dim, self.dims))

        result = {}
        for k,v in self.__m__.items():
            for i in range(size):
                j = k[dim]+i
                new_k = k[:dim] + (j // size,) + k[dim+1:]

                if not new_k in result:
                    result[new_k] = numpy.zeros(numpy.array(self.shape)*size, dtype = numpy.complex)

                j = j % size
                result[new_k][i*self.shape[0]:(i+1)*self.shape[0],j*self.shape[1]:(j+1)*self.shape[1]] = v

        return TightBinding(result)

    # ==================================================================
    # Methods returning devices.
    # ==================================================================

    def subdevice(self, center, leads, leads2):
        """
        Creates a device from this TightBinding.

        Args:

            center (array): center indexes;

            leads (list): list of leads indexes inside center region;

            leads2 (list): list of next lead replica indexes;

        Returns:

            A MultiterminalDevice specified by the input indexes.
        """

        center = self.__s2m__(center,0)

        center_t = self.subsystem(center,center).simplified()
        leads_t = []
        connections_t = []

        for l_i,(l,l2) in enumerate(zip(leads,leads2)):

            l = self.__s2m__(l,0)
            l2 = self.__s2m__(l2,0)

            # Lead
            l_t = self.subsystem(l,l).simplified()

            # Lead replica
            l2_t = self.subsystem(l2,l2).simplified()

            # Hoppings
            l_b = self.subsystem(l,l2).simplified()
            l_f = self.subsystem(l2,l).simplified()

            # Check if the lead belongs to the scattering region
            if not set(l).issubset(set(center)):
                raise ValueError("Lead {:d} does not completely belong to the central region".format(l_i))

            # Check absence of overlap
            center_without_lead = set(center) - set(l)
            if self.subsystem(center_without_lead,l2).absmax() != 0 or self.subsystem(l2,center_without_lead).absmax() != 0:
                raise ValueError("Lead replica {:d} interacts with the central region".format(l_i))

            leads_t.append(TightBinding({
                0: l_t,
                1: l_f,
                -1: l_b,
            }))

            c = numpy.eye(self.shape[0], dtype = numpy.int)
            connections_t.append(c[numpy.ix_(l, center)])

        return MultiterminalDevice(center_t, leads_t, connections_t)

    def periodic_device(self, size = 0):
        """
        Initializes a periodic device along dimension 0 with 2 identical
        leads and no scattering region.

        Kwargs:

            size (int): number of units to include into a scattering
            region;

        Returns:

            A periodic 2-terminal device.
        """
        d = self.shape[0]
        s = self.super(size+4, 0)
        return s.subdevice(
            slice(d,-d),
            [slice(d,2*d),slice(d*(2+size),d*(3+size))],
            [slice(0,d),slice(d*(3+size),d*(4+size))],
        ).squeezed(0)

class MultiterminalDevice(LinearOperator):

    def __init__(self, center, leads, connections):
        """
        Describes a multiterminal device.

        Args:

            center (TightBinding): a matrix of the center part of device;

            leads (array): semi-infinite leads connected to a device,
            an array where each element is a TightBinding;

            connections (array): corresponding connections of the leads,
            each element is an array leads size by center size.
        """

        if not center.is_square():
            raise ValueError("Center is not a square TightBinding")

        self.center = center.copy()
        self.dims = center.dims
        self.leads = []
        self.connections = []

        for l,c in zip(leads,connections):
            self.append_lead(l,c)

    # ==================================================================
    # Serialization
    # ==================================================================

    def __getstate__(self):
        return dict(
            dims=self.dims,
            center=self.center,
            leads=self.leads,
            connections=self.connections,
        )

    def __setstate__(self, data):
        self.dims = data["dims"]
        self.center = data["center"]
        self.leads = data["leads"]
        self.connections = data["connections"]

    # ==================================================================
    # Constructors
    # ==================================================================

    @staticmethod
    def zeros(dims,shape_c,shape_l,connections):
        """
        Prepares a zero MTD of the given geometry.

        Args:

            dims (int): number of dimensions of the central part;

            shape_c (tuple of 2 ints): shape of the center matrix;

            shape_l (array of tuples of 2 ints): shapes of leads'
            matrices;

            connections (array): corresponding connections of the leads,
            each element is an array leads size by center size.

        Returns:

            A zero MTD.
        """
        return MultiterminalDevice(
            TightBinding.zeros(dims,shape_c),
            tuple(TightBinding.zeros(dims+1,i) for i in shape_l),
            connections,
        )

    @staticmethod
    def eye(dims,shape_c,shape_l,connections):
        """
        Prepares an "eye" MTD of the given geometry.

        Args:

            dims (int): number of dimensions of the central part;

            shape_c (tuple of 2 ints): shape of the center matrix;

            shape_l (array of tuples of 2 ints): shapes of leads'
            matrices;

            connections (array): corresponding connections of the leads,
            each element is an array leads size by center size.

        Returns:

            An "eye" MTD.
        """
        return MultiterminalDevice(
            TightBinding.eye(dims,shape_c),
            tuple(TightBinding.eye(dims+1,i) for i in shape_l),
            connections,
        )

    # ==================================================================
    # Helpers
    # ==================================================================

    def test_compatible(self, other):
        """
        Tests whether this and other device are compatible to perform
        simple arithmetics. Raises an exception otherwise.
        """
        if not self.center.dims == other.center.dims:
            raise ValueError("Devices have different dimensionality: {:d} vs {:d}".format(self.center.dims, other.center.dims))

        if not self.center.shape == other.center.shape:
            raise ValueError("The shape of the central part does not match: {} vs {}".format(self.center.shape, other.center.shape))

        if not len(self.leads) == len(other.leads):
            raise ValueError("The number of leads is different: {:d} vs {:d}".format(len(self.leads), len(other.leads)))

        for n, (i,j) in enumerate(zip(self.leads, other.leads)):
            if not i.shape == j.shape:
                raise ValueError("The shape of a lead {:d} does not match: {} vs {}".format(n,i.shape,j.shape))

        for n, (i,j) in enumerate(zip(self.connections, other.connections)):
            if not numpy.array_equal(i,j):
                raise ValueError("The connections arrays for lead {:d} are not equal".format(n))

    # ==================================================================
    # Required by parent
    # ==================================================================

    def zeros_like(self):
        """
        Prepares a zero MTD of the same geometry as self.

        Returns:

            A zero MTD.
        """
        return MultiterminalDevice.zeros(
            self.dims,
            self.center.shape,
            tuple(i.shape for i in self.leads),
            self.connections,
        )

    def eye_like(self):
        """
        Prepares an "eye" MTD of the same geometry as self.

        Returns:

            An "eye" MTD.
        """
        return MultiterminalDevice.eye(
            self.dims,
            self.center.shape,
            tuple(i.shape for i in self.leads),
            self.connections,
        )

    def tr(self):
        """
        Calculates transposed MultiterminalDevice.

        Returns:

            A transposed of the MultiterminalDevice.
        """
        return MultiterminalDevice(
            self.center.tr(),
            list(i.tr() for i in self.leads),
            list(i.copy() for i in self.connections),
        )

    def cc(self):
        """
        Calculates complex conjugate MultiterminalDevice.

        Returns:

            A complex conjugate of the MultiterminalDevice.
        """
        return MultiterminalDevice(
            self.center.cc(),
            list(i.cc() for i in self.leads),
            list(i.conj() for i in self.connections),
        )

    def absmax(self):
        """
        Retrieves the maximum value by modulus.

        Returns:

            A float with the maximum value.
        """
        result = max(list(i.absmax() for i in self.leads))
        result = max(result, self.center.absmax())
        for l, c in zip(self.leads, self.connections):
            for k,v in l.__m__.items():
                result = max(result, abs(v.dot(c)).max())
        return result

    def __neg__(self):
        return MultiterminalDevice(
            -self.center,
            list(-i for i in self.leads),
            list(i.copy() for i in self.connections),
        )

    def __abs__(self):
        return MultiterminalDevice(
            abs(self.center),
            list(abs(i) for i in self.leads),
            list(i.copy() for i in self.connections),
        )

    def __add__(self, other):

        if isinstance(other, MultiterminalDevice):

            self.test_compatible(other)
            return MultiterminalDevice(
                self.center + other.center,
                list(i+j for i,j in zip(self.leads, other.leads)),
                list(i.copy() for i in self.connections),
            )

        else:

            return super(MultiterminalDevice, self).__add__(other)

    def __mul__(self, other):
        return MultiterminalDevice(
            self.center*other,
            list(i*other for i in self.leads),
            list(i.copy() for i in self.connections),
        )

    # ==================================================================
    # Lead operations
    # ==================================================================

    def append_lead(self, lead, connection):
        """
        Appends a lead to the device.

        Args:

            lead (TightBinding): a lead to add;

            connection (matrix): connection matrix of the lead. None is
            accepted to assume a zero matrix;
        """

        if not lead.dims == self.dims+1:
            raise ValueError("Inconsistent dimensions of the lead: {:d} instead of {:d}".format(lead.dims, self.dims+1))

        if not lead.is_square():
            raise ValueError("Lead is not a square TightBinding: {}".format(str(lead.shape)))

        if connection is None:
            connection = numpy.zeros((lead.shape[0], self.center.shape[1]), dtype = numpy.int)

        else:
            if connection.shape != (lead.shape[0], self.center.shape[1]):
                raise ValueError("Expected an array of size {:d}x{:d} for the connection, found {}".format(lead.shape[0], self.center.shape[1], repr(connection.shape)))

        self.leads.append(lead)
        self.connections.append(connection)

    def remove_lead(self, lead):
        """
        Removes one of the leads from the device.

        Args:

            lead (int): a lead to remove;
        """
        self.leads = self.leads[:lead] + self.leads[lead+1:]
        self.connections = self.connections[:lead] + self.connections[lead+1:]

    def consume_lead(self, lead):
        """
        Increases the device by one unit of one of the leads.

        Args:

            lead (int): a lead to include.
        """
        l = self.leads[lead]
        c = self.connections[lead]

        center = self.center.foreach(lambda k,v: (k, numpy.pad(v,((0,l.shape[0]),(0,l.shape[1])), mode = 'constant')))
        for k,v in center.__m__.items():
            v[self.center.shape[0]:,self.center.shape[1]:] = l[(0,)+k]
            v[self.center.shape[0]:,:self.center.shape[1]] = l[(1,)+k].dot(c)
            v[:self.center.shape[0],self.center.shape[1]:] = c.T.dot(l[(-1,)+k])
        self.center = center

        self.connections = list(numpy.pad(i, ((0,0),(0,center.shape[1] - i.shape[1])), mode = 'constant') for i in self.connections)
        self.connections[lead] = numpy.eye(l.shape[0], M = center.shape[0], k = center.shape[0]-l.shape[0])

    # ==================================================================
    # Arithmetic operations affecting geometry
    # ==================================================================

    def squeezed(self, *dims):
        """
        Squeezes dimensions from this MTD.

        Args:

            dimension(s) to be squeezed or "all".

        Returns:

            An MTD with dimensions squeezed out.
        """
        if len(dims) == 1 and dims[0] == "all":

            return self.squeezed(*tuple(i for i in range(self.dims) if self.is_isolated(i)))

        else:

            for i in dims:
                if not self.is_isolated(i):
                    raise ValueError("This MTD is not isolated along dimension {:d}".format(i))

            keep = numpy.ones(self.dims, dtype = bool)
            keep2 = numpy.ones(self.dims+1, dtype = bool)
            dims = numpy.array(dims, dtype = int)
            keep[dims] = False
            keep2[dims+1] = False

            return MultiterminalDevice(
                self.center.foreach(lambda k,v: (numpy.array(k)[keep], v)),
                tuple(i.foreach(lambda k,v: (numpy.array(k)[keep2], v)) for i in self.leads),
                self.connections,
            )

    def fourier(self, k, index = None):
        """
        Performs Fourier transform.

        Args:

            k (float,array): the wave number or a wave vector.

        Kwargs:

            index (int): index to transform; if None transforms along
            first ``len(k)'' indeces.

        Returns:

            Transformed MultiterminalDevice.
        """
        if self.dims == 1 and index is None:
            index = 0

        if not index is None:

            return MultiterminalDevice(
                self.center.fourier(k, index = index),
                list(i.fourier(k, index = index+1) for i in self.leads),
                list(i.copy() for i in self.connections),
            )

        else:

            result = self
            for i in k:
                result = result.fourier(i, index = 0)
            return result

    # ==================================================================
    # Classification
    # ==================================================================

    def is_isolated(self, d):
        """
        Determines if this is an isolated device along the dimension of
        the central region specified.

        Args:

            d (int): dimension to check.

        Returns:

            True if is isolated.
        """
        if not self.center.is_isolated(d):
            return False

        for i in self.leads:
            if not i.is_isolated(d+1):
                return False

        return True

    # ==================================================================
    # Utility
    # ==================================================================

    def __repr__(self):
        return "MultiterminalDevice dim {:d} leads {:d}".format(self.dims, len(self.leads))

    def copy(self):
        """
        Calculates a deep copy.

        Returns:

            A copy.
        """
        return MultiterminalDevice(
            self.center.copy(),
            list(i.copy() for i in self.leads),
            list(i.copy() for i in self.connections),
        )

#class BlochCalculator(object):

    #def __init__(self, lead):
        #"""
        #A Bloch function based calculator for the transport properties of
        #the lead.

        #Args:

            #lead (TightBinding): lead to calculate.
        #"""
        #self.lead = lead

        #if not self.lead.is_1DNN():
            #raise ValueError("Not a 1D NN tight binding")
        #if not self.lead.is_square():
            #raise ValueError("Not a square TightBinding")

        #self.states_r = None
        #self.states_ri = None
        #self.w_r = None
        #self.w_ri = None
        #self.v_r = None

        #self.states_l = None
        #self.states_li = None
        #self.w_l = None
        #self.w_li = None
        #self.v_l = None

    #def calculate_bloch_states(self):
        #"""
        #Calculates bloch states. Stores right-going modes
        #into 'self.states_r' and corresponding bloch coefficients |w|<1
        #into 'self.w_r' (inverse in 'self.w_ri'). The left-going modes
        #are stored in 'self.states_l', their bloch coefficients |w|>1 in
        #'self.w_l' and their inverse in 'self.w_li'. Their velocities
        #are stored in
        #"""

        #N = self.lead.shape[0]
        #lhs = numpy.zeros((2*N,2*N), dtype = numpy.complex)
        #rhs = numpy.zeros((2*N,2*N), dtype = numpy.complex)

        #lhs[:N,N:] = numpy.eye(N)
        #lhs[N:,:N] = self.lead[-1]
        #lhs[N:,N:] = self.lead[0]

        #rhs[:N,:N] = numpy.eye(N)
        #rhs[N:,N:] = -self.lead[1]

        #w, r = linalg.eig(lhs, b = rhs)

        ## Select
        #r[:N,numpy.isinf(abs(w))] = r[N:,numpy.isinf(abs(w))]
        #r = r[:N,:]

        ## Normalize
        #r /= ((abs(r)**2).sum(axis = 0)**.5)[numpy.newaxis,:]

        ## Inverse
        #w_i = numpy.where(numpy.isinf(w), 0, w**-1)
        ## Note: if w == 0 then corresponding r is an eigenstate of
        ## self.lead[-1] with zero eigenvalue; if w == inf then r is
        ## an eigenstate of self.lead[1] with a zero eigenvalue.

        ## Group velocity
        #v = []

        #for i in range(len(w)):

            #if w[i] == 0:
                #v.append(1j)

            #elif w_i[i] == 0:
                #v.append(-1j)

            #else:
                ##m = 1j*( - w[i]*self.lead[1] + w_i[i]*self.lead[-1] )
                #m = 1j*( - w_i[i]*self.lead[-1] + w[i]*self.lead[1] )
                #v.append(numpy.conj(r[:,i]).dot(m).dot(r[:,i]))

        #v = numpy.array(v)

        #left = abs(w)>1
        #right = abs(w)<1

        #if not sum(left) == N:
            #raise Exception("Could not identify left-going modes (found: {:d}, expected: {:d}). Increase the imaginary part.".format(sum(left), N))

        #if not sum(right) == N:
            #raise Exception("Could not identify right-going modes (found: {:d}, expected: {:d}). Increase the imaginary part.".format(sum(right), N))

        #self.states_r, self.w_r, self.w_ri, self.v_r = r[:,right], w[right], w_i[right], v[right]
        #self.states_l, self.w_l, self.w_li, self.v_l = r[:,left], w[left], w_i[left], v[left]

        #self.states_ri = numpy.linalg.inv(self.states_r)
        #self.states_li = numpy.linalg.inv(self.states_l)

    #def bloch_matrix(self, kind, power):
        #"""
        #Calculates the Bloch matrix.

        #Args:

            #kind (str): either "right", "+", or "left", "-";

            #power (int): power of the Bloch matrix.

        #Returns:

            #The power of a bloch matrix.
        #"""
        #if self.states_r is None:
            #self.calculate_bloch_states()

        #r,ri,w,w_i = {
            #"right": (self.states_r, self.states_ri, self.w_r, self.w_ri),
            #"left":  (self.states_l, self.states_li, self.w_l, self.w_li),
        #}[{
            #"left": "left",
            #"right": "right",
            #"-": "left",
            #"+": "right",
        #}[kind]]

        #if power>0:

            #if sum(w_i == 0)>0:
                #raise ValueError("Requested Bloch matrix is infinite")

            #return r.dot(numpy.diag(w**power)).dot(ri)

        #elif power<0:

            #if sum(w == 0)>0:
                #raise ValueError("Requested Bloch matrix is infinite")

            #return r.dot(numpy.diag(w_i**-power)).dot(ri)

        #else:

            #return r.dot(ri)

    #def self_energy(self, conj = False):
        #"""
        #Calculates self-energy.

        #Kwargs:

            #conj (bool): whether to calculate conjugated self-energy

        #Returns:

            #Self-energy.
        #"""

        #if not conj:
            #return -self.lead[-1].dot(self.bloch_matrix("-",-1))

        #else:
            #return -self.lead[1].dot(self.bloch_matrix("-",1))
            ##return self.lead[-1].dot(self.bloch_matrix("+",-1))

    #def gamma(self):
        #"""
        #Calculates gamma matrix.

        #Returns:

            #Self-energy.
        #"""
        ##return 1j*( self.self_energy() - self.self_energy(conj = True) ).dot(self.states_r).dot(numpy.diag(abs(abs(self.w_r)-1)>1e-5)).dot(self.states_ri)
        #se = self.self_energy()
        #return 1j*( se - se.conj().T )

class GreensFunctionCalculator(object):

    def __init__(self, lead):
        """
        A Greens function based calculator for the transport properties of
        the lead.

        Args:

            lead (TightBinding): lead to calculate.
        """
        self.lead = lead

        if not self.lead.is_1DNN():
            raise ValueError("Not a 1D NN tight binding")
        if not self.lead.is_square():
            raise ValueError("Not a square TightBinding")

        self.gf_r = None

    def gf(self, tolerance = None, maxiter = 1000):
        """
        Calculates the retarded Green's function matrix and stores it into
        'self.gf_r'.

        Kwargs:

            tolerance (float): tolerance for Green's function iterations;

            maxiter (int): maximum number of iterations.

        Returns:

            A surface Green's function matrix.
        """
        w0 = self.lead[0]
        w1 = self.lead[-1]
        w2 = self.lead[1]

        if tolerance is None:
            tolerance = 1e-10 * max(max(abs(w0.max()), abs(w1).max()), abs(w2).max())

        self.gf_r = greens_function(
            w0,
            w1,
            w2,
            tolerance,
            maxiter,
        )

        return self.gf_r

    def bloch_matrix(self):
        """
        Calculates the Bloch matrix.

        Returns:

            The Bloch matrix.
        """
        if self.gf_r is None:
            self.gf()

        return -self.gf_r.dot(self.lead[1])

    def self_energy(self):
        """
        Calculates self-energy.

        Returns:

            Self-energy calculated.
        """
        if self.gf_r is None:
            self.gf()

        return self.lead[-1].dot(self.gf_r).dot(self.lead[1])

    def gamma(self):
        """
        Calculates Gamma-function.

        Returns:

            Gamma function calculated.
        """
        se = self.self_energy()
        return 1j*(se - se.conj().T)

    def __(self):
        if gamma:

            gmm = w1.dot(gf).dot(w2)
            gmm = 1j*(gmm - gmm.conj().T)
            rtn.append(gmm)

        if gamma_:

            gf2 = numpy.linalg.inv(gf.dot(w1).dot(w2))
            rtn.append(1j*w1.dot(gf-gf2).dot(w2))

            #bm = -gf.dot(w2)
            #val, vec = numpy.linalg.eig(bm)
            #vec_i = numpy.linalg.inv(vec)
            #val_i = numpy.where(abs(val) <1e-14,0,val**-1)
            #bm2 = vec.dot(numpy.diag(val_i)).dot(vec_i)
            #rtn.append(1j*(w1.dot(bm) - w2.dot(bm2)))
            ##rtn.append(numpy.linalg.inv(gf.dot(bm)))
            ##rtn.append(1j*(w1.dot(bm) - w1.dot(bm2)))
            ##rtn.append(1j*(w1.dot(bm) - numpy.linalg.inv(gf)))
            ##se = w1.dot(gf).dot(w2)
            ##rtn.append(1j*(-se - w2.dot(w1).dot(numpy.linalg.inv(se))))
            if direction == "positive":
                rtn[-1] = -rtn[-1]

        if len(rtn) == 1:
            return rtn[0]

        else:
            return tuple(rtn)

class MTDCalculator(object):

    def __init__(self, device, calculator = GreensFunctionCalculator):
        """
        A calculator for the transport properties of a multiterminal
        device.

        Args:

            device (MultiterminalDevice): device to calculate.
        """

        if not device.dims == 0:
            raise ValueError("{:d} periodic dimension(s) have to be eliminated before using this calculator".format(device.dims))

        self.device = device
        self.gf_r = None
        self.calculators = list(calculator(i) for i in device.leads)

    def lead_self_energy(self,lead,**kwargs):
        """
        Returns a lead self energy adapted to the size of the scattering
        region matrix.

        Kwargs:

            the kwargs are passed to the leads' calculators

        Returns:

            A matrix with the lead self energy.
        """
        return self.device.connections[lead].conj().T.dot(self.calculators[lead].self_energy()).dot(self.device.connections[lead])

    def gf(self, **kwargs):
        """
        Calculates the Green's function matrix and stores it into
        'self.gf_r'.

        Kwargs:

            the kwargs are passed to the leads' calculators

        Returns:

            A matrix with the Green's function.
        """

        gi = self.device.center.diagonal.copy()

        for i in range(len(self.device.leads)):

            gi -= self.lead_self_energy(i, **kwargs)

        self.gf_r = linalg.inv(gi)

        return self.gf_r

    def transmission(self, source, drain):
        """
        Calculates a two-terminal transmission.

        Args:

            source (int): source of Bloch waves;

            drain (int): drain of Bloch waves;

        Returns:

            The total transmission.
        """
        if self.gf_r is None:
            self.gf()

        gamma_source = self.calculators[source].gamma()
        gamma_drain = self.calculators[drain].gamma()
        G_sd = self.device.connections[drain].dot(self.gf_r).dot(self.device.connections[source].conj().T)

        return numpy.trace(gamma_source.dot(G_sd.conj().T).dot(gamma_drain).dot(G_sd))

    def tmatrix(self, source, drain):
        """
        Calculates a transmission matrix.

        Args:

            source (int): source of Bloch waves;

            drain (int): drain of Bloch waves;

        Returns:

            A transmission matrix.
        """
        if self.gf_r is None:
            self.gf()

        G_sd = self.device.connections[drain].dot(self.gf_r).dot(self.device.connections[source].conj().T)

        gamma_source = self.calculators[source].gamma()
        gamma_drain = self.calculators[drain].gamma()

        vals, vecs = linalg.eig(gamma_source)
        vald, vecd = linalg.eig(gamma_drain)

        return numpy.diag((-self.calculators[drain].v_r.real)**.5).\
            dot(self.calculators[drain].states_ri).dot(G_sd).\
            dot(self.calculators[source].states_li.conj().T).\
            dot(numpy.diag(self.calculators[source].v_l.real**.5))

