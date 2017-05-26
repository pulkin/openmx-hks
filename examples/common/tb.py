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

class TightBinding(object):
    
    def __init__(self, m, vectors = None):
        """
        A class representing tight binding (periodic) matrix.
        
        Args:
        
            m (dict or list): a tight binding matrix. In the case of a
            dict, the dict keys represent integers corresponding to
            matrix block location while dict values are block matrices.
            Otherwise the list items are tight-binding block matrices and
            the keyword argument ``vectors'' providing block indices has
            to be supplied separately.
            
        Kwargs:
        
            vectors (list): a list of tight binding block indices for the
            case of ``m'' being a list.
        """
        if len(m) == 0:
            raise ValueError("Empty input")
        
        if isinstance(m,dict):
            if not vectors is None:
                warn("The vectors keyword argument is ignored")
        
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
        
        if t == numpy.ndarray:
            
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
    
    def copy(self):
        """
        Calculates a copy.
        
        Returns:
        
            A (shallow) copy.
        """
        return TightBinding(dict((k,v.copy()) for k,v in self.__m__.items()))
        
    def __repr__(self):
        return "TightBinding dim {:d} [{:d}x{:d}]".format(self.dims, self.shape[0], self.shape[1])
        
    @staticmethod
    def __tr_i__(i):
        return tuple(-ii for ii in i)
        
    def __eq__(self, other):
        return (self-other).absmax() == 0
        
    def __neg__(self):
        return self.foreach(lambda k,v: (k,-v))
        
    def __abs__(self):
        return self.foreach(lambda k,v: (k,abs(v)))
            
    def __add__(self, other):
        if isinstance(other, TightBinding):
            
            keys = set(self.__m__.keys()) | set(other.__m__.keys())
            result = {}
            for k in keys:
                result[k] = self.__m__.get(k,0) + other.__m__.get(k,0)
            return TightBinding(result)
            
        else:
            result = self.copy()
            result.diagonal += other*numpy.eye(self.shape[0], M = self.shape[1])
            return result
        
    def __radd__(self, other):
        return self + other
        
    def __sub__(self, other):
        return self + (-other)
        
    def __rsub__(self, other):
        return - self + other
        
    def __mul__(self, other):
        return self.foreach(lambda k,v: (k,v*other))
        
    def __rmul__(self, other):
        return self*other
        
    def __div__(self, other):
        if isinstance(other, TightBinding):
            keys = set(self.__m__.keys()) | set(other.__m__.keys())
            result = {}
            for k in keys:
                result[k] = self.__m__.get(k,0) / other.__m__.get(k,0)
            return TightBinding(result)
        else:
            return self.foreach(lambda k,v: (k,v/other))
            
    def __getattr__(self, name):
        
        if name == "diagonal":
            key = (0,)*self.dims
            if key in self.__m__:
                return self[key]
            else:
                return numpy.zeros(self.shape)
            
        else:
            raise AttributeError
            
    def __setattr__(self, name, value):
        
        if name == "diagonal":
            key = (0,)*self.dims
            self.__m__[key] = value
            
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
            if not isinstance(k, int):
                raise ValueError("The keys should be ints, found {} instead".format(str(type(k))))
            
        item = numpy.array(item, dtype = numpy.complex)
        
        if not len(item.shape) == 2:
            raise ValueError("Not a 2D matrix: shape = {}".format(str(item.shape)))
            
        if not tuple(item.shape) == self.shape:
            raise ValueError("Wrong dimensions: shape = {}".format(str(item.shape)))
                
        self.__m__[key] = item
    
    def __s2m__(self, s, d):
        result = numpy.zeros(self.shape[d], dtype = bool)
        result[s] = True
        return result
    
    def foreach(self, p):
        """
        Performs array-wise operation on this TightBinding.
        
        Args:
        
            p (function): a function accepting key and value and returning
            a new key and a new value.
            
        Returns:
        
            A new tightbinding with a given operation applied.
        """
        result = {}
        for k,v in self.__m__.items():
            r = p(k,v)
            if not r is None and not r[0] is None and not r[1] is None:
                result[r[0]] = r[1]
        return TightBinding(result)
        
    def size(self):
        """
        Calculates total number of entries in the TightBinding.
        
        Returns:
        
            A number of entries in the TightBinding.
        """
        return self.diagonal.size*len(self.__m__)
        
    def isnan(self, check_inf = True):
        """
        Checks if NaN are present in the tight binding.
        
        Kwargs:
        
            check_inf (bool): checks infinities as well.
            
        Returns:
        
            A Tight binding with empty matrices where corresponding
            (NaN, inf) elements are set to one.
        """
        return self.foreach(lambda k,v: (k, numpy.isnan(v) + (check_inf == True) * numpy.isinf(v)))
        
    def is_square(self):
        """
        Checks if it is a square matrix.
        
        Returns:
        
            True if it is a square matrix.
        """
        return self.shape[0] == self.shape[1]
        
    def tr(self):
        """
        Calculates transposed matrix.
        
        Returns:
        
            A transposed of the tight binding matrix.
        """
        return self.foreach(lambda k,v: (TightBinding.__tr_i__(k), numpy.transpose(v)))
        
    def cc(self):
        """
        Calculates complex conjugate.
        
        Returns:
        
            A complex conjugate of the tight binding matrix.
        """
        return self.foreach(lambda k,v: (k, numpy.conj(v)))
        
    def simplified(self):
        """
        Simplifies this TightBinding by taking out all zero matrices.
        
        Returns:
        
            A simplified TightBinding.
        """
        return self.foreach(lambda k,v: (k,v) if numpy.count_nonzero(v)>0 or k == (0,)*self.dims else None)

    def simplify(self):
        """
        Simplifies this TightBinding in-place by removing all-zero submatrices.
        """
        new = {}
        for k,v in self.__m__.items():
            if numpy.any(v!=0):
                new[k] = v
        self.__m__ = new
        
    def squeezed(self):
        """
        Squeezes decoupled dimensions from this tight-binding.
        
        Returns:
        
            A tight-binding with squeezed dimensions.
        """
        v = numpy.array(self.__m__.keys())
        unique = tuple(tuple(numpy.unique(x).tolist()) for x in v.T)
        keep = tuple(i != (0,) for i in unique)
        return self.foreach(lambda k,v: (
            tuple(i for i,j in zip(k,keep) if j),
            v
        ))
        
    def hc(self):
        """
        Calculates Hermitian conjugate.
        
        Returns:
        
            A Hermitian conjugate of the tight binding matrix.
        """
        return self.tr().cc()
        
    def absmax(self, location = False):
        """
        Retrieves maximum value by modulus.
        
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
            
    def hermitian_error(self):
        """
        Calculates largest non-hermitian element of the Hamiltonian.
        
        Returns:
        
            A float giving a measure of the matrix non-hermitianity.
        """
        if not self.is_square():
            raise ValueError("Not a square matrix")
            
        return (self - self.hc()).absmax()
        
    def fourier(self, k, index = None):
        """
        Performs Fourier transform.
        
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
        
    def is_trivial(self):
        """
        Determines if this is a trivial (0-dim) TightBinding.
        
        Returns:
        
            True if trivial.
        """
        return self.dims == 0

    def is_zero(self):
        """
        Determines if zero.

        Returns:

            True if contains zero matrices only.
        """
        for v in self.__m__.values():
            if numpy.any(v!=0):
                return False
        return True
        
    def is_isolated(self, dim = None):
        """
        Determines if this is an isolated TightBinding. If ``dim'' is
        provided checks only apply to the specified dimension.
        
        Kwargs:
        
            dim (int): dimension to check.
            
        Returns:
        
            True if is isolated.
        """
        if dim is None:
            return len(self.__m__) == 0 or (len(self.__m__) == 1 and (0,)*self.dims in self.__m__)
            
        else:
            for k in self.__m__.keys():
                if not k[dim] == 0:
                    return False
                    
            return True
        
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
        result = []
        for i in pts:
            
            hh = self.fourier(i).__m__.values()[0]
            
            if not b is None:
                bb = b.fourier(i).__m__.values()[0]
                
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
        
    def eye(self):
        """
        Generates an eye tight binding matrix.
            
        Returns:
        
            A TightBinding of the same structure as input.
        """
        return TightBinding({
            (0,)*self.dims : numpy.eye(self.shape[0], M = self.shape[1]),
        })
        
    def zeros(self):
        """
        Generates a zeros tight binding matrix.
            
        Returns:
        
            A TightBinding of the same structure as input.
        """
        return TightBinding({
            (0,)*self.dims : numpy.zeros(self.shape),
        })
        
    def super(self, size, dim):
        """
        Creates a supercell version of this tight binding.
        
        Args:
        
            size (int): a multiple to increase the tight binding;
            
            dim (int): dimension along which to increase the tight binding.
            
        Returns:
        
            A supercell tight binding.
        """
        if not dim < self.dims:
            raise ValueError("Wrong dimension: {:d} for {:d}-dim TightBinding".format(dim, self.dims))
            
        result = {}
        for k,v in self.__m__.items():
            for i in range(size):
                j = k[dim]+i
                new_k = k[:dim] + (j/size,) + k[dim+1:]
                
                if not new_k in result:
                    result[new_k] = numpy.zeros(numpy.array(self.shape)*size, dtype = numpy.complex)
                
                j = j % size
                result[new_k][i*self.shape[0]:(i+1)*self.shape[0],j*self.shape[1]:(j+1)*self.shape[1]] = v
                
        return TightBinding(result)
        
    def inverted(self, dim = 0):
        """
        Creates a space-inverted version of self.
        
        Kwargs:
        
            dim (int): dimension to invert;
            
        Returns:
        
            An inverted version of self.
        """
        return self.foreach(lambda k,v: (k[:dim] + (-k[dim],) + k[dim+1:],v))
        
    def subsystem(self, rows, columns):
        """
        Creates a subsystem of this TightBinding by applying a given
        selection to each matrix in the Hamiltonian.
        
        Args:
        
            rows (array): row selection;
            
            columns (array): columns selection;
            
        Returns:
        
            A subsystem TghtBinding selected.
        """
        target = numpy.ix_(self.__s2m__(rows,0),self.__s2m__(columns,1))
        return self.foreach(lambda k,v: (k, v[target]))
               
    def shift_subblock(self, rows, columns, vector):
        """
        Shifts TightBinding submatrices along a vector.
        
        Args:
        
            rows (array): submatrix rows mask;
            
            columns (array): submatrix columns mask;
            
            vector (array): an array of integers specifying the submatrix
            shift.
        """
        move = self.subsystem(rows, columns)
        if not move.is_zero():
            move = move.simplified()
        else:
            move = None
        target = numpy.ix_(self.__s2m__(rows,0),self.__s2m__(columns,1))
        if not move is None:

            for v in self.__m__.values():
                v[target] = 0

            self.simplify()

            for k,v in move.__m__.items():
                if numpy.any(v!=0):
                    new_k = tuple(i+j for i,j in zip(k,vector))
                    m = self[new_k]
                    m[target] = v
                    self[new_k] = m
            
    def hermitian_shift_subblock(self, selection, vector):
        """
        Performs a shift of a subblock of a square Hermitian matrix such
        that effective TightBinding remains the same.
        
        Args:
        
            selection (array): diagonal submatrix mask;
            
            vector (array): an array of integers specifying the submatrix
            shift.
        """
        if not self.is_square():
            raise ValueError("Not a square matrix")
            
        vector = numpy.array(vector)
        self.shift_subblock(selection, ~selection, vector)
        self.shift_subblock(~selection, selection, -vector)
        
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
        
        center_t = self.subsystem(center,center).simplified()
        leads_t = []
        connections_t = []
        
        for l_i,(l,l2) in enumerate(zip(leads,leads2)):
            
            # Lead
            l_t = self.subsystem(l,l).simplified()
            
            # Lead replica
            l2_t = self.subsystem(l2,l2).simplified()
            
            # Hoppings
            l_b = self.subsystem(l,l2).simplified()
            l_f = self.subsystem(l2,l).simplified()
            
            # Check absence of overlap
            for another_i, another in enumerate(leads2):
                if not l_i == another_i and not (self.subsystem(l,another).absmax() == 0 and self.subsystem(another,l).absmax() == 0):
                    raise ValueError("Leads {:d} and {:d} have non-zero hopping".format(l_i, another_i))
                    
            leads_t.append(TightBinding({
                0: l_t,
                1: l_f,
                -1: l_b,
            }))
            
            c = numpy.zeros((l_t.shape[0], center_t.shape[0]), dtype = numpy.int)
            indexes = numpy.zeros(self.shape[0], dtype = numpy.bool)
            indexes[l] = True
            indexes = indexes[center]
            c[numpy.arange(l_t.shape[0]),indexes] = 1
            connections_t.append(c)
            
        return MultiterminalDevice(center_t, leads_t, connections_t)
        
    def periodic_device(self, size = 0, dim = 0):
        """
        Initializes a periodic device with 2 identical leads and no
        scattering region.
        
        Kwargs:
        
            size (int): number of units to include into a scattering
            region;
            
            dim (int): dimension to repeat along;
        
        Returns:
        
            A periodic 2-terminal device.
        """
        d = self.shape[0]
        return self.super(size+4, dim).subdevice(
            slice(d,-d),
            [slice(d,2*d),slice(d*(2+size),d*(3+size))],
            [slice(0,d),slice(d*(3+size),d*(4+size))],
        ).remove_dim(dim)
        
class BlochCalculator(object):
    
    def __init__(self, lead):
        """
        A Bloch function based calculator for the transport properties of
        the lead.
        
        Args:
        
            lead (TightBinding): lead to calculate.
        """
        self.lead = lead
                
        if not self.lead.is_1DNN():
            raise ValueError("Not a 1D NN tight binding")
        if not self.lead.is_square():
            raise ValueError("Not a square TightBinding")
        
        self.states_r = None
        self.states_ri = None
        self.w_r = None
        self.w_ri = None
        self.v_r = None
        
        self.states_l = None
        self.states_li = None
        self.w_l = None
        self.w_li = None
        self.v_l = None
        
    def calculate_bloch_states(self):
        """
        Calculates bloch states. Stores right-going modes
        into 'self.states_r' and corresponding bloch coefficients |w|<1
        into 'self.w_r' (inverse in 'self.w_ri'). The left-going modes
        are stored in 'self.states_l', their bloch coefficients |w|>1 in
        'self.w_l' and their inverse in 'self.w_li'. Their velocities
        are stored in 
        """
        
        N = self.lead.shape[0]    
        lhs = numpy.zeros((2*N,2*N), dtype = numpy.complex)
        rhs = numpy.zeros((2*N,2*N), dtype = numpy.complex)
        
        lhs[:N,N:] = numpy.eye(N)
        lhs[N:,:N] = self.lead[-1]
        lhs[N:,N:] = self.lead[0]
        
        rhs[:N,:N] = numpy.eye(N)
        rhs[N:,N:] = -self.lead[1]
        
        w, r = linalg.eig(lhs, b = rhs)
        
        # Select
        r[:N,numpy.isinf(abs(w))] = r[N:,numpy.isinf(abs(w))]
        r = r[:N,:]
        
        # Normalize
        r /= ((abs(r)**2).sum(axis = 0)**.5)[numpy.newaxis,:]
        
        # Inverse
        w_i = numpy.where(numpy.isinf(w), 0, w**-1)
        # Note: if w == 0 then corresponding r is an eigenstate of
        # self.lead[-1] with zero eigenvalue; if w == inf then r is
        # an eigenstate of self.lead[1] with a zero eigenvalue.
        
        # Group velocity
        v = []
        
        for i in range(len(w)):
            
            if w[i] == 0:
                v.append(1j)
            
            elif w_i[i] == 0:
                v.append(-1j)
                
            else:
                #m = 1j*( - w[i]*self.lead[1] + w_i[i]*self.lead[-1] )
                m = 1j*( - w_i[i]*self.lead[-1] + w[i]*self.lead[1] )
                v.append(numpy.conj(r[:,i]).dot(m).dot(r[:,i]))
        
        v = numpy.array(v)
        
        left = abs(w)>1
        right = abs(w)<1
        
        if not sum(left) == N:
            raise Exception("Could not identify left-going modes (found: {:d}, expected: {:d}). Increase the imaginary part.".format(sum(left), N))
            
        if not sum(right) == N:
            raise Exception("Could not identify right-going modes (found: {:d}, expected: {:d}). Increase the imaginary part.".format(sum(right), N))
        
        self.states_r, self.w_r, self.w_ri, self.v_r = r[:,right], w[right], w_i[right], v[right]
        self.states_l, self.w_l, self.w_li, self.v_l = r[:,left], w[left], w_i[left], v[left]
        
        self.states_ri = numpy.linalg.inv(self.states_r)
        self.states_li = numpy.linalg.inv(self.states_l)
        
    def bloch_matrix(self, kind, power):
        """
        Calculates the Bloch matrix.
        
        Args:
        
            kind (str): either "right", "+", or "left", "-";
            
            power (int): power of the Bloch matrix.
            
        Returns:
        
            The power of a bloch matrix.
        """
        if self.states_r is None:
            self.calculate_bloch_states()
            
        r,ri,w,w_i = {
            "right": (self.states_r, self.states_ri, self.w_r, self.w_ri),
            "left":  (self.states_l, self.states_li, self.w_l, self.w_li),
        }[{
            "left": "left",
            "right": "right",
            "-": "left",
            "+": "right",
        }[kind]]
        
        if power>0:
            
            if sum(w_i == 0)>0:
                raise ValueError("Requested Bloch matrix is infinite")
                
            return r.dot(numpy.diag(w**power)).dot(ri)
            
        elif power<0:
            
            if sum(w == 0)>0:
                raise ValueError("Requested Bloch matrix is infinite")
                
            return r.dot(numpy.diag(w_i**-power)).dot(ri)
            
        else:
            
            return r.dot(ri)
            
    def self_energy(self, conj = False):
        """
        Calculates self-energy.
        
        Kwargs:
        
            conj (bool): whether to calculate conjugated self-energy
            
        Returns:
        
            Self-energy.
        """
        
        if not conj:
            return -self.lead[-1].dot(self.bloch_matrix("-",-1))
            
        else:
            return -self.lead[1].dot(self.bloch_matrix("-",1))
            #return self.lead[-1].dot(self.bloch_matrix("+",-1))
        
    def gamma(self):
        """
        Calculates gamma matrix.
            
        Returns:
        
            Self-energy.
        """
        #return 1j*( self.self_energy() - self.self_energy(conj = True) ).dot(self.states_r).dot(numpy.diag(abs(abs(self.w_r)-1)>1e-5)).dot(self.states_ri)
        se = self.self_energy()
        return 1j*( se - se.conj().T )
                
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

class MultiterminalDevice(object):
    
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
        
        self.dims = center.dims
        
        if not center.is_square():
            raise ValueError("Center is not a square TightBinding")
        
        for i_i,i in enumerate(leads):
            if not i.dims == self.dims+1:
                raise ValueError("Inconsistent dimensions of lead {:d}: {:d} instead of {:d}".format(i_i, i.dims, self.dims+1))
            if not i.is_square():
                raise ValueError("Lead is not a square TightBinding: {}".format(str(i.shape)))
        
        self.center = center
        self.leads = leads
            
        for i_c, c in enumerate(connections):
            if not c.shape == (leads[i_c].shape[0], center.shape[1]):
                raise ValueError("Expected an array of size {:d}x{:d} for connection {:d}, found {:d}x{:d}".format(leads[i_c].shape[0], center.shape[1], i_c, c.shape[0], c.shape[1]))
                
        self.connections = list(numpy.array(i, dtype = numpy.complex) for i in connections)
        
        self.gf_r = None
        self.lead_calculators = None
        
    def test_compatible(self, other):
        """
        Tests whether this and other device are arithmetically
        compatible. Raises exception otherwise.
        """
        if not self.center.dims == other.center.dims:
            raise Exception("Devices have different dimensionality")
            
        if not self.center.shape == other.center.shape:
            raise Exception("The shape of the central part does not match")
            
        for n, (i,j) in enumerate(zip(self.leads, other.leads)):
            if not i.shape == j.shape:
                raise Exception("The shape of a lead {:d} does not match".format(n))
                
        for n, (i,j) in enumerate(zip(self.connections, other.connections)):
            if not numpy.array_equal(i,j):
                raise Exception("The connections arrays for lead {:d} are not equal".format(n))
                
    def foreach(self, p, include_connections = False):
        """
        Performs a matrix-wise operation on this MultiterminalDevice.
        
        Args:
        
            p (function): a function accepting value and returning
            a new value.
            
        Kwargs:
        
            include_connections (bool): includes connections.
            
        Returns:
        
            A new MultiterminalDevice with a given operation applied.
        """
        f = lambda k,v: (k,p(v))
        
        center = self.center.foreach(f)
        
        leads = []
        for l in self.leads:
            leads.append(l.foreach(f))
            
        if include_connections:
            connections = []
            for a in self.connections:
                connections.append(p(a))
                
        else:
            connections = list(i.copy() for i in self.connections)
            
        return MultiterminalDevice(center,leads,connections)
                
    def __eq__(self, other):
        
        try:
            self.test_compatible(other)
        except:
            return False
            
        if not self.center == other.center:
            return False
            
        for i,j in zip(self.leads, other.leads):
            if not i==j:
                return False
                
        return True
        
    def __neg__(self):
        return self.foreach(lambda v: -v)
        
    def __abs__(self):
        return self.foreach(lambda v: abs(v))
            
    def __add__(self, other):
        if isinstance(other, MultiterminalDevice):
            self.test_compatible(other)
            return MultiterminalDevice(
                self.center + other.center,
                list(i+j for i,j in zip(self.leads, other.leads)),
                list(i.copy() for i in self.connections),
            )
            
        else:
            return self + other * self.eye()
        
    def __radd__(self, other):
        return self + other
        
    def __sub__(self, other):
        return self + (-other)
        
    def __rsub__(self, other):
        return - self + other
        
    def __mul__(self, other):
        return self.foreach(lambda v: v*other)
        
    def __rmul__(self, other):
        return self*other
        
    def __div__(self, other):
        return self.foreach(lambda v: v/other)
            
    def add_lead(self, lead):
        """
        Includes a particular lead into the scattering region.
        
        Args:
        
            lead (int): a lead to include.
            
        Returns:
        
            A new instance with a larger central region.
        """
        l = self.leads[lead]
        c = self.connections[lead]
        
        center = self.center.foreach(lambda k,v: (k, numpy.pad(v,((0,l.shape[0]),(0,l.shape[1])), mode = 'constant')))
        for k,v in center.__m__.items():
            v[self.center.shape[0]:,self.center.shape[1]:] = l[(0,)+k]
            v[self.center.shape[0]:,:self.center.shape[1]] = l[(1,)+k].dot(c)
            v[:self.center.shape[0],self.center.shape[1]:] = c.T.dot(l[(-1,)+k])
        
        connections = list(numpy.pad(i, ((0,0),(0,center.shape[1] - i.shape[1])), mode = 'constant') for i in self.connections)
        connections[lead] = numpy.eye(l.shape[0], M = center.shape[0], k = center.shape[0]-l.shape[0])
        
        return MultiterminalDevice(center, self.leads, connections)
    
    def is_isolated(self, dim):
        """
        Determines if this is an isolated device along the dimension of
        the central region specified.
        
        Args:
        
            dim (int): dimension to check.
            
        Returns:
        
            True if is isolated.
        """
        if not self.center.is_isolated(dim = dim):
            return False
            
        for i in self.leads:
            if not i.is_isolated(dim = dim+1):
                return False
                
        return True
        
    def remove_dim(self, dim):
        """
        Removes the dimension if the device is isolated along it.
        
        Args:
        
            dim (int): a dimension of the central region to get rid of;
            
        Returns:
        
            A new device with the dimension removed.
        """
        if not self.is_isolated(dim):
            raise ValueError("Not isolated along dim {:d}".format(dim))
            
        takeaway_center = lambda k,v: (k[:dim]+k[dim+1:], v)
        takeaway_lead = lambda k,v: (k[:dim+1]+k[dim+2:], v)
        
        return MultiterminalDevice(
            self.center.foreach(takeaway_center),
            list(i.foreach(takeaway_lead) for i in self.leads),
            list(self.connections),
        )
        
    def eye(self):
        """
        Generates an eye multiterminal device.
            
        Returns:
        
            A MultiterminalDevice of the same structure as self.
        """
        return MultiterminalDevice(
            self.center.eye(),
            list(i.eye() for i in self.leads),
            list(i.copy() for i in self.connections),
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
        print vals, vald
        
        return numpy.diag((-self.calculators[drain].v_r.real)**.5).\
            dot(self.calculators[drain].states_ri).dot(G_sd).\
            dot(self.calculators[source].states_li.conj().T).\
            dot(numpy.diag(self.calculators[source].v_l.real**.5))
    
    def __(self):
        if True:
            pass
        else:
            print "---------------"
            print "G"
            print G
            print "---------------"
            print "Caroli"
            print numpy.trace(leads_gamma[source].dot(G_sd.conj().T).dot(leads_gamma[drain]).dot(G_sd))
            
            w1, r1, v1 = self.leads[source].bloch_states(return_velocities = True)
            w2, r2, v2 = self.leads[drain].bloch_states(return_velocities = True)
            #propagating1 = numpy.abs(v1) > 1e-6
            #propagating2 = numpy.abs(v2) > 1e-6
            #positive1 = numpy.logical_or(
                #numpy.logical_and(
                    #numpy.logical_not(propagating1),
                    #numpy.abs(w1) > 1
                #),
                #numpy.logical_and(
                    #propagating1,
                    #numpy.real(v1) > 0
                #)
            #)
            #positive2 = numpy.logical_or(
                #numpy.logical_and(
                    #numpy.logical_not(propagating2),
                    #numpy.abs(w2) > 1
                #),
                #numpy.logical_and(
                    #propagating2,
                    #numpy.real(v2) > 0
                #)
            #)
            positive1 = numpy.abs(w1)<1
            positive2 = numpy.abs(w2)>1
            #print "V1", numpy.real(v1)
            #print propagating1
            print "W1", numpy.abs(w1)
            print "--", positive1
            #print "V2", numpy.real(v2)
            #print propagating2
            print "W2", numpy.abs(w2)
            print "--", positive2
            w1, r1, v1 = w1[positive1], r1[:,positive1], v1[positive1]
            w2, r2, v2 = w2[positive2], r2[:,positive2], v2[positive2]
            
            r1d = linalg.inv(r1)
            r2d = linalg.inv(r2)
            
            print "Reliable"
            rel = -abs(r2d.dot(G_sd).dot(leads_gamma[source]).dot(r1))**2*v2[:,numpy.newaxis].real/v1[numpy.newaxis,:].real
            print rel
            print "Reliable2"
            result = -abs(r2d.dot(G_sd).dot(r1d.conj().T))**2*v2[:,numpy.newaxis].real*v1[numpy.newaxis,:].real
            print result
            print "Reliable3"
            t = numpy.diag(-v2**.5).dot(r2d).dot(G_sd).dot(r1d.conj().T).dot(numpy.diag(v1**.5))
            rel = abs(t)**2
            print rel
            print numpy.trace(t.dot(t.conj().T))
            print "G_L n"
            print self.leads[source].gf(energy, b = b.leads[source])
            print "V_L"
            print v1
            print "TEST V_L"
            print r1.conj().T.dot(leads_gamma[source]).dot(r1)
            print "Gamma_L"
            print leads_gamma[source]
            print "TEST Gamma_L"
            tgl = r1d.conj().T.dot(numpy.diag(v1)).dot(r1d)
            tgr = r2d.conj().T.dot(numpy.diag(v2)).dot(r2d)
            print tgl
            print leads_gamma[source]/tgl
            print "V_R"
            print v2
            print "TEST V_R"
            print r2.conj().T.dot(leads_gamma[drain]).dot(r2)
            print "Gamma_R"
            print leads_gamma[drain]
            print "TEST Gamma_R"
            print tgr
            print leads_gamma[drain]/tgr
            print "TEST Caroli"
            print numpy.trace(tgl.dot(G_sd.conj().T).dot(tgr).dot(G_sd))
            print "TEST"
            test = r2d.dot(G_sd).dot(r1d.conj().T)
            test = test.dot(numpy.diag(abs(v1))).dot(test.conj().T).dot(numpy.diag(abs(v2)).conj())
            print numpy.trace(test)
            print "components"
            print r2d.dot(G_sd).dot(r1d.conj().T)
            print v1
            print v2
            print "---------------"
            
            #print r1d.dot(r1d.conj().T)
            #vv = v1[numpy.newaxis,:]*v2[:,numpy.newaxis]
            #print "VV"
            #print vv
            #print "<>"
            #print r2d.dot(G_sd).dot(r1d.conj().T)
            #print "ans"
            #print abs(r2d.dot(G_sd).dot(r1d.conj().T))**2*vv
            ##print v2
            
            return {
                "transmission-matrix": result,
                "incoming-modes": r1,
                "outgoing-modes": r2,
                "incoming-v": v1,
                "outgoing-v": v2,
            }
