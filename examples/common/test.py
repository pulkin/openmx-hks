#!/usr/bin/env python

# ======================================================================
# This file is under construction.
# ======================================================================

import unittest

import numpy
from numpy import testing, random
import pickle

from tb import *

class LOTests(unittest.TestCase):
    
    def test_method_list(self):
        """
        This is simply a list of methods to implement.
        """
        o = LinearOperator()
        with self.assertRaises(NotImplementedError):
            o.zeros_like()
        with self.assertRaises(NotImplementedError):
            o.eye_like()
        with self.assertRaises(NotImplementedError):
            o.tr()
        with self.assertRaises(NotImplementedError):
            o.cc()
        with self.assertRaises(NotImplementedError):
            o.absmax()
        with self.assertRaises(NotImplementedError):
            -o
        with self.assertRaises(NotImplementedError):
            abs(o)
        with self.assertRaises(NotImplementedError):
            o+object()
        with self.assertRaises(NotImplementedError):
            o.zeros_like()
        with self.assertRaises(NotImplementedError):
            o*3
        with self.assertRaises(NotImplementedError):
            o/object()
            
class TightBindingTests(unittest.TestCase):
    
    def setUp(self):
        self.r1 = TightBinding({
            (1,): [[1]],
            (-1,): [[1]],
        })
        
        self.c1 = TightBinding({
            (0,): [[0]],
            (1,): [[1j]],
            (-1,): [[-1j]],
        })
        
        self.r2 = TightBinding({
            (1,): [[1,0],[0,1]],
            (-1,): [[1,0],[0,1]],
        })
        
        self.c2 = TightBinding({
            (0,): [[0,0],[0,0]],
            (1,): [[1j,0],[0,1j]],
            (-1,): [[-1j,0],[0,-1j]],
        })
        
        self.r5 = TightBinding({
            (0,): numpy.diag(numpy.arange(5)),
            (1,): 0.1*numpy.eye(5),
            (-1,): 0.1*numpy.eye(5),
        })
        
        self.c5 = TightBinding({
            (0,): numpy.diag(numpy.arange(5)),
            (1,): 0.1j*numpy.eye(5),
            (-1,): -0.1j*numpy.eye(5),
        })
        
        self.a = [self.r1,self.c1,self.r2,self.c2,self.r5,self.c5]
        
    # ==================================================================
    # Serialization
    # ==================================================================
    
    def test_equality(self):
        for i in self.a:
            assert pickle.loads(pickle.dumps(i)) == i
        
    # ==================================================================
    # Constructors
    # ==================================================================
    
    def test_constructors(self):
        for dim in range(4):
            for s in range(3,5):
                
                t = TightBinding.zeros(dim, (s,s))
                assert t.dims == dim
                assert t.shape == (s,s)
                assert t.size() == 0
                
                t = TightBinding.eye(dim, (s,s))
                assert t.dims == dim
                assert t.shape == (s,s)
                assert t.size() == 1
                testing.assert_array_equal(t.diagonal, numpy.eye(s))
                
    # ==================================================================
    # Arithmetics
    # ==================================================================
        
    def test_arithmetics(self):
        for i in self.a:
            assert i.zeros_like() == 0
            assert i.eye_like() == TightBinding.eye(i.dims, i.shape)
            assert i.tr() == i.hc().cc()
            assert i.cc() == i.hc().tr()
            assert i.hc() == i.tr().cc()
            assert i.hc() == i
            assert i.absmax() == {
                self.r1: 1,
                self.c1: 1,
                self.r2: 1,
                self.c2: 1,
                self.r5: 4,
                self.c5: 4,
            }[i]
            
            assert i == i
            assert (-i) == (-i)
            assert i + (-i) == 0
            
            assert abs(i).absmax() == i.absmax()
            
            assert 2*i == i + i
            assert 0.5*i == i - 0.5*i
            assert 0.5*i == i/2
            assert i + 3.14*i.eye_like() == i + 3.14 == 3.14 + i
            assert i - 3.14*i.eye_like() == i - 3.14 == - (3.14 - i)

    def test_absmax(self):
        N = 10
        tb = TightBinding(dict(zip(range(N),([[i*0.1j]] for i in range(N)))))
        assert tb.absmax() == 0.1*(N-1)
        assert tb.absmax(location = True) == (0.1*(N-1), (N-1,), (0,0))
            
    # ==================================================================
    # Getters/setters and diagonal
    # ==================================================================
                
    def test_get_zero(self):
        for i in self.a:
            assert numpy.all(i[1000] == 0)
            with self.assertRaises(ValueError):
                i[3,5]
                
    def test_set(self):
        x = TightBinding({0:[[0]]})
        assert numpy.all(x[1] == 0)
        x[1] = [[2]]
        assert numpy.all(x[1] == 2)
        with self.assertRaises(ValueError):
            x[1,2] = [[2]]
        with self.assertRaises(ValueError):
            x['a'] = [[2]]
        with self.assertRaises(ValueError):
            x[1] = [2]
        with self.assertRaises(ValueError):
            x[1] = [[0,0],[0,0]]
            
    def test_diagonal(self):
        for i in self.a:
            testing.assert_array_equal(i.diagonal, i[0])
            if (0,) in i.__m__:
                assert i.diagonal is i[0]
            else:
                assert numpy.count_nonzero(i.diagonal) == 0
            i.diagonal = numpy.eye(i.shape[0])
            testing.assert_array_equal(i.diagonal, numpy.eye(i.shape[0]))
            
    def test_attr(self):
        with self.assertRaises(AttributeError):
            self.r1.__some_non_existing_attribute__
            
    # ==================================================================
    # Arithmetic operations preserving geometry
    # ==================================================================
    
    def test_nan(self):
        for i in self.a:
            for m in i.nan().__m__.values():
                assert numpy.all(m == 0)

    # ==================================================================
    # Arithmetic operations affecting geometry
    # ==================================================================
        
    def test_squeezed(self):
        t = TightBinding({
            (0, 1): [[1]],
            (0, 0): [[0]],
            (0,-1): [[1]],
        })
        assert t.dims == 2
        t1 = t.squeezed(0)
        assert t1.dims == 1
        t2 = t.squeezed("all")
        assert t2.dims == 1
        with self.assertRaises(ValueError):
            t.squeezed(1)
            
    def test_insert_dimension(self):
        t = self.r1.insert_dimension(0)
        assert t.dims == 2
        assert t.squeezed("all") == self.r1
            
    def test_fourier_gamma(self):
        for i in self.a:
            t = i.fourier(0)
            tt = 0
            for ii in i.__m__.values():
                tt += ii
            testing.assert_array_equal(t.diagonal, tt)
            
    def test_fourier_edge(self):
        for i in self.a:
            t = i.fourier(0.5)
            tt = 0
            for kk, ii in i.__m__.items():
                kk = kk[0]
                if kk % 2 == 0:
                    tt += ii
                else:
                    tt -= ii
            testing.assert_allclose(t.diagonal, tt, atol = 1e-15)
                
    # ==================================================================
    # Classification
    # ==================================================================
    
    def test_classification(self):
        for i in self.a:
            assert i.is_square()
            assert i.is_1DNN()
            assert i.is_hermitian()
            assert not i.is_trivial()
            assert not i.is_isolated()
            
            i = i.fourier(0)
            
            assert i.is_trivial()
            assert i.is_isolated()
            
    def test_1dnn(self):
        not_1d = self.r1.insert_dimension(0)
        assert not not_1d.is_1DNN()
        not_nn = self.r1.copy()
        not_nn[2] = [[0.1]]
        assert not not_nn.is_1DNN()
            
    # ==================================================================
    # Utility
    # ==================================================================
    
    def test_utility(self):
        for i in self.a:
            assert isinstance(repr(i), str)
            assert "TightBinding" in repr(i)
            assert i.copy() == i
            assert not i.copy() is i
            assert i.size() == {
                self.r1: 2,
                self.c1: 3,
                self.r2: 2,
                self.c2: 3,
                self.r5: 3,
                self.c5: 3,
            }[i]
            
    def test_simplify(self):
        assert self.c1.size() == 3
        assert self.c1.simplified().size() == 2
        self.c1.simplify()
        assert self.c1.size() == 2
        
    def test_inverted(self):
        assert self.c1.inverted() == self.c1.cc()
        
    def test_subsystem(self):
        for i in self.a:
            assert i.subsystem(numpy.arange(i.shape[0]),numpy.arange(i.shape[0])) == i
        assert self.c2.subsystem([0],[0]) == self.c1
        
    def test_subsystem_order(self):
        s = self.r5.subsystem((3, 2), (3, 2))
        testing.assert_equal(numpy.diag(s.diagonal), (3, 2))
            
    def test_eigpath(self):
        for tb in self.a:
            for i in (tb, tb.inverted()):
                p = numpy.linspace(-.5,.5,31)
                eig = i.eig_path(p)
                a = numpy.diag(i[0])[numpy.newaxis,:]
                b = numpy.diag(i[1])[numpy.newaxis,:]
                testing.assert_allclose(eig,
                    numpy.sort(a + 2*numpy.cos(2*numpy.pi*p)[:,numpy.newaxis]*b.real - 2*numpy.sin(2*numpy.pi*p)[:,numpy.newaxis]*b.imag, axis = 1)
                )
                eig2 = i.eig_path(p, b = i.eye_like())
                testing.assert_allclose(eig,eig2)
                
    def test_eigpath_2D(self):
        a = 1
        b = .2
        c = .3
        d = .4
        tb = TightBinding({
            (0,0):[[a]],
            (0,1):[[b+1j*c]],
            (0,-1):[[b-1j*c]],
            (1,0):[[d]],
            (-1,0):[[d]],
        })
        for i in (tb, tb.inverted()):
            p = (numpy.linspace(-.5,.5,31)[:,numpy.newaxis]*[[1,1]])
            eig = i.eig_path(p)
            testing.assert_allclose(eig,
                (a +\
                2*numpy.cos(2*numpy.pi*p[:,1])*b -\
                2*numpy.sin(2*numpy.pi*p[:,1])*c +\
                2*numpy.cos(2*numpy.pi*p[:,0])*d)[:,numpy.newaxis],
            )
            
    def test_shift_sb(self):
        x = self.r5.copy()
        x.shift_subblock([1,2],[1,2],(1,))
        ref = TightBinding({
            (0,): numpy.diag([0,0.1,0.1,3,4]),
            (1,): numpy.diag([0.1,1,2,0.1,0.1]),
            (2,): numpy.diag([0,0.1,0.1,0,0]),
            (-1,): numpy.diag([0.1,0,0,0.1,0.1]),
        })
        assert x == ref
        x.shift_subblock([False,True,True,False,False], [False,True,True,False,False], (-1,))
        assert self.r5 == x
            
    def test_herm_shift(self):
        for i in self.a:
            p = numpy.linspace(-.5,.5,31)
            eig = i.eig_path(p[:,numpy.newaxis])
            mask = numpy.zeros(i.shape[0], dtype = numpy.bool)
            mask[::2] = True
            i.hermitian_shift_subblock(mask, (3,))
            eig2 = i.eig_path(p[:,numpy.newaxis])
            testing.assert_allclose(eig, eig2)
        x = self.r5.subsystem([0,1,2],[0,1])
        x.shift_subblock([0,1],[0,1],(1,))
        with self.assertRaises(ValueError):
            x.hermitian_shift_subblock([0,1],(1,))

    def test_herm_shift_2(self):
        sample = TightBinding({
            0: [
                [1, 2, 3],
                [4, 5, 6],
                [7, 8, 9],
                ],
            1: [
                [.1, .2, .3],
                [.4, .5, .6],
                [.7, .8, .9],
                ]
            })
        sample.hermitian_shift_subblock([1, 2], (1,))
        reference = TightBinding({
            -1:[
                [0, 2, 3],
                [0, 0, 0],
                [0, 0, 0],
                ],
            0: [
                [1, .2, .3],
                [0, 5, 6],
                [0, 8, 9],
                ],
            1: [
                [.1, 0, 0],
                [4, .5, .6],
                [7, .8, .9],
                ],
            2: [
                [0, 0, 0],
                [.4, 0, 0],
                [.7, 0, 0],
                ],
        })
        assert sample == reference

    def test_super(self):
        x = self.c1.super(2,0)
        assert x == TightBinding({
            (0,): [[0,1j],[-1j,0]],
            (1,): [[0,0],[1j,0]],
            (-1,): [[0,-1j],[0,0]],
        })
        for i in self.a:
            s = i.super(2,0)
            assert i == i.hc()
            eig = i.eig_path(numpy.linspace(0,1,30,endpoint = False)[:,numpy.newaxis])
            eig2 = s.eig_path(numpy.linspace(0,1,15,endpoint = False)[:,numpy.newaxis])
            testing.assert_allclose(numpy.sort(numpy.concatenate((eig[:15],eig[15:]),axis = 1),axis = 1),eig2, atol = 1e-14)
        with self.assertRaises(ValueError):
            i.super(2,1)
            
    def test_subdevice(self):
        large = TightBinding([
                [1,.1,0,0,0],
                [.1,1,.1,0,0],
                [0,.1,0,.1j,0],
                [0,0,-.1j,-1,.1j],
                [0,0,0,-.1j,-1],
            ],
        )
        d = large.subdevice([1,2,3],[[1],[3]],[[0],[4]])
        assert d.center == TightBinding([
            [1,.1,0],
            [.1,0,.1j],
            [0,-.1j,-1],
        ])
        assert d.leads[0] == TightBinding({
            0:[[1]],
            1:[[.1]],
            -1:[[.1]],
        })
        assert d.leads[1] == TightBinding({
            0:[[-1]],
            1:[[-.1j]],
            -1:[[.1j]],
        })
        testing.assert_array_equal(d.connections[0],[[1,0,0]])
        testing.assert_array_equal(d.connections[1],[[0,0,1]])
        with self.assertRaises(ValueError):
            large.subdevice([2,3],[[1],[3]],[[0],[4]])
        large.diagonal[0,3] = 1
        with self.assertRaises(ValueError):
            large.subdevice([1,2,3],[[1],[3]],[[0],[4]])
            
    def test_subdevice_order(self):
        large = TightBinding([
            [0, 1, 0, 0, 0, 0, 0, 0],
            [1, 2, 3, 0, 0, 0, 0, 0],
            [0, 3, 0, 1, 0, 0, 0, 0],
            [0, 0, 1, 2, 4, 0, 0, 0],
            [0, 0, 0, 4, 5, 6, 0, 0],
            [0, 0, 0, 0, 6, 7, 8, 0],
            [0, 0, 0, 0, 0, 8, 7, 8],
            [0, 0, 0, 0, 0, 0, 8, 7],
        ])
        lead1 = [ (1, 0), (3, 2) ]
        center = (2, 3, 6, 5, 4)
        lead2 = [ (6,), (7,) ]
        d = large.subdevice(center, (lead1[1], lead2[0]), (lead1[0], lead2[1]))
        assert d.center == TightBinding([
            [0, 1, 0, 0, 0],
            [1, 2, 0, 0, 4],
            [0, 0, 7, 8, 0],
            [0, 0, 8, 7, 6],
            [0, 4, 0, 6, 5],
        ])
        
        assert d.leads[0] == TightBinding({
            0: [[2, 1], [1, 0]],
            1: [[0, 3], [0, 0]],
            -1: [[0, 0], [3, 0]],
        })
        
        assert d.leads[1] == TightBinding({
            0: [[7]],
            1: [[8]],
            -1: [[8]],
        })
        
        testing.assert_array_equal(d.connections[0], [
            [0, 1, 0, 0, 0],
            [1, 0, 0, 0, 0],
        ])
        testing.assert_array_equal(d.connections[1], [
            [0, 0, 1, 0, 0],
        ])
            
    def test_pd(self):
        x = self.c1.periodic_device()
        assert x.center == TightBinding([[0,1j],[-1j,0]])
        assert x.leads[0] == x.leads[1].inverted() == self.c1
        
        y = self.c1.insert_dimension(0).periodic_device()
        assert y.center == TightBinding({
            1:[[1j,0],[0,1j]],
            -1:[[-1j,0],[0,-1j]],
        })
        assert y.leads[0] == y.leads[1] == self.c1.insert_dimension(0)
            
class MTDTests(unittest.TestCase):
    
    def setUp(self):
        self.p1 = TightBinding({
            0 : [[0]],
            1 : [[1]],
            -1: [[1]],
        }).periodic_device()
        self.p2 = TightBinding({
            0: [[1]],
        }).periodic_device()
        self.p3 = TightBinding({
            0 : [[1,0],[0,-1]],
            1 : [[.1,0],[0,.1j]],
            -1: [[.1,0],[0,-.1j]],
        }).periodic_device()
        self.p4 = TightBinding({
            (1,0)  : [[1]],
            (-1,0) : [[1]],
            (0,1)  : [[.1j]],
            (0,-1) : [[-.1j]],
        }).periodic_device()
        self.p_all = (self.p1,self.p2,self.p3,self.p4)
    
    # ==================================================================
    # Serialization
    # ==================================================================
    
    def test_equality(self):
        for i in self.p_all:
            assert pickle.loads(pickle.dumps(i)) == i
            
    # ==================================================================
    # Constructors
    # ==================================================================
    
    def test_constructors(self):
        
        z = MultiterminalDevice.zeros(0,(5,5),[],[])
        assert z.dims == 0
        assert z.center.shape == (5,5)
        assert z.center == 0
        assert len(z.leads) == 0
        assert len(z.connections) == 0
        
        z = MultiterminalDevice.eye(1,(5,5),[(3,3),(2,2)],[
            numpy.eye(3,5),
            numpy.eye(2,5),
        ])
        assert z.dims == 1
        assert z.center.shape == (5,5)
        assert z.center == z.center.eye_like()
        assert len(z.leads) == 2
        for l in z.leads:
            assert l == l.eye_like()
        assert numpy.count_nonzero(z.connections[0]) == 3
        assert numpy.count_nonzero(z.connections[1]) == 2
                
    # ==================================================================
    # Arithmetics
    # ==================================================================
        
    def test_arithmetics(self):
        for i in self.p_all:
            assert i.zeros_like() == 0
            assert i.eye_like() == MultiterminalDevice.eye(
                i.dims,
                i.center.shape,
                tuple(j.shape for j in i.leads),
                i.connections,
            )
            assert i.tr() == i.hc().cc()
            assert i.cc() == i.hc().tr()
            assert i.hc() == i.tr().cc()
            assert i.hc() == i
            assert i.absmax() == 1
            
            assert i == i
            assert (-i) == (-i)
            assert i + (-i) == 0
            
            assert abs(i).absmax() == i.absmax()
            
            assert 2*i == i + i
            assert 0.5*i == i - 0.5*i
            assert 0.5*i == i/2
            assert i + 3.14*i.eye_like() == i + 3.14 == 3.14 + i
            assert i - 3.14*i.eye_like() == i - 3.14 == - (3.14 - i)
            
    # ==================================================================
    # Compatability exceptions
    # ==================================================================
    
    def test_compat(self):
        self.p1.test_compatible(self.p2)
        with self.assertRaises(ValueError):
            self.p1.test_compatible(self.p3)
        with self.assertRaises(ValueError):
            self.p1.test_compatible(self.p4)
        p1c = self.p1.copy()
        p1c.remove_lead(1)
        with self.assertRaises(ValueError):
            self.p1.test_compatible(p1c)
        p1c.append_lead(self.p3.leads[0],None)
        with self.assertRaises(ValueError):
            self.p1.test_compatible(p1c)
        p3c = p1c.copy()
        p3c.connections[1][0,0] = 1
        with self.assertRaises(ValueError):
            p3c.test_compatible(p1c)
            
    # ==================================================================
    # Lead operations
    # ==================================================================
    
    def test_lead_op(self):
        for p in self.p_all:
            p_bup = p.copy()
            p.append_lead(p.leads[0],p.connections[0])
            assert len(p.leads) == 3
            assert len(p.connections) == 3
            assert p.leads[0] == p.leads[2]
            testing.assert_array_equal(p.connections[0], p.connections[2])
            with self.assertRaises(ValueError):
                p.append_lead(p.leads[0].insert_dimension(0),p.connections[0])
            with self.assertRaises(ValueError):
                p.append_lead(TightBinding.eye(1,(3,2)),None)
            with self.assertRaises(ValueError):
                p.append_lead(p.leads[0],numpy.eye(5))
            p.remove_lead(2)
            assert p == p_bup
    
    def test_consume_lead(self):
        assert self.p1.center == TightBinding([
            [0,1],
            [1,0],
        ])
        self.p1.consume_lead(0)
        assert self.p1.center == TightBinding([
            [0,1,1],
            [1,0,0],
            [1,0,0],
        ])
        testing.assert_array_equal(self.p1.connections[0], [[0,0,1]])
        
    # ==================================================================
    # Arithmetic operations affecting geometry
    # ==================================================================
        
    def test_squeezed(self):
        x = TightBinding({
            (1,0,0):[[1]],
            (-1,0,0):[[1]],
            (0,1,0):[[.3]],
            (0,-1,0):[[.3]],
        }).periodic_device()
        with self.assertRaises(ValueError):
            x.squeezed(0)
        x1 = x.squeezed(1)
        assert x1.center.dims == 1
        x2 = x.squeezed("all")
        assert x1==x2
        
    #def test_isolated(self):
        #assert not self.p4.is_isolated(0)
        #self.p4.center = TightBinding.zeros(1,(1,1))
        #assert not self.p4.is_isolated(0)
        
    ##def test_gf(self):
        ##for i in self.a:
            ##for e in numpy.linspace(-1,1,30)+1e-1j:
                ##w = e-i
                ##gfc = GreensFunctionCalculator(w)
                ##g = gfc.gf()
                ##se = gfc.self_energy()
                ##testing.assert_allclose((w[0] - se).dot(g), numpy.eye(g.shape[0]))
                
    ##def test_gamma(self):
        ##for i in self.a:
            ##for e in numpy.linspace(-1,1,30)+1e-3j:
                ##g, gm = (e-i).gf(gamma_ = True)
                ##w,r,vel = (e-i).bloch_states(return_velocities = True)
                ##s = abs(w)>1
                ##testing.assert_allclose(vel[s],numpy.diag(r[:,s].conj().T.dot(gm).dot(r[:,s])))
            
#class LeadCalculatorTest(unittest.TestCase):
    
    #def setUp(self):
        
        #self.tb = [
            #TightBinding({
                #(0,): [[0,1+.2j],[1-.2j,0]],
                #(1,): [[0,0.5],[.1-.2j,0.3]],
                #(-1,): [[0,.1+.2j],[0.5,0.3]],
            #}),
            #TightBinding({
                #(0,): [[0,1],[1,0]],
                #(1,): [[0,0],[.1,0]],
                #(-1,): [[0,.1],[0,0]],
            #}),
        #]
        
        #self.energies = [.1j, 1+.1j, 1+1e-8j]
            
    #def test_bloch(self):
        
        #for tb in self.tb:
            #for e in self.energies:
            
                #tbb = e - tb
                #calc = BlochCalculator(tbb)
                #calc.calculate_bloch_states()
                
                #for r,w,vel in ((calc.states_l, calc.w_l, calc.v_l), (calc.states_r, calc.w_r, calc.v_r)):
                    
                    #aw = abs(w)
                    #testing.assert_allclose((abs(r)**2).sum(axis = 0),1)
                    
                    #for i in range(len(w)):
                        
                        #if aw[i] == float("inf"):
                            #pass
                        
                        #elif aw[i] != aw[i]:
                            #raise ValueError("NaN encountered")
                            
                        #elif aw[i] == 0:
                            #pass
                            
                        #else:
                            
                            #k = (-1j*numpy.log(abs(w)) + numpy.angle(w))/numpy.pi/2
                            #m = tb.fourier(k[i]).diagonal
                            #v = r[:,i]
                        
                            #testing.assert_allclose(tb[-1].dot(v) + w[i]*tb[0].dot(v) + w[i]**2*tb[1].dot(v),e*v*w[i], atol = 1e-14)
                            #testing.assert_allclose(m.dot(v),e*v, atol = 1e-14)
                    
    #def test_gf(self):
        
        #for tb in self.tb:
            #for e in self.energies:
            
                #g1 = GreensFunctionCalculator(e - tb).gf()
                #g2 = GreensFunctionCalculator(e - tb.super(2,0)).gf()
                #testing.assert_allclose(g1,g2[2:,2:])
                #w = e - tb
                #testing.assert_allclose(g1, numpy.linalg.inv(w[0] - w[-1].dot(g1).dot(w[1])))
        
    #def test_bloch_matrix(self):
        
        #for tb in self.tb:
            #for e in self.energies:
            
                #tbb = e - tb
                #calc = BlochCalculator(tbb)
                #calc2 = GreensFunctionCalculator(tbb)
                
                #testing.assert_allclose(calc.bloch_matrix("-",-1), calc2.bloch_matrix(), atol = 1e-14)
                
                #testing.assert_allclose(calc.bloch_matrix("-",0), numpy.eye(2), atol = 1e-14)
                #testing.assert_allclose(calc.bloch_matrix("+",0), numpy.eye(2), atol = 1e-14)
                
                #testing.assert_allclose(tbb[-1].dot(calc.bloch_matrix("-",-2)) + tbb[0].dot(calc.bloch_matrix("-",-1)),-tbb[1], atol = 1e-14)
                #testing.assert_allclose(tbb[1].dot(calc.bloch_matrix("+",2)) + tbb[0].dot(calc.bloch_matrix("+",1)),-tbb[-1], atol = 1e-14)
            
    #def test_self_energy(self):
            
        #for tb in self.tb:
            #for e in self.energies:
            
                #tbb = e - tb

                #c1 = GreensFunctionCalculator(tbb)
                #se1 = c1.self_energy()
                
                #c2 = BlochCalculator(tbb)
                #se2 = c2.self_energy()
                
                #testing.assert_allclose(se1,se2, atol = 1e-14)
                
    #def self_energy_conj_test(self):
        
        #for tb in self.tb:
            #for e in numpy.linspace(1,0,3)+1e-8j:
                
                #tbb = e - tb
                #c = BlochCalculator(tbb)
                #se = c.self_energy().dot(c.states_r).dot(numpy.diag(abs(abs(c.w_r)-1)>1e-5)).dot(c.states_ri)
                #sec = c.self_energy(conj = True).dot(c.states_r).dot(numpy.diag(abs(abs(c.w_r)-1)>1e-5)).dot(c.states_ri)
                
                #testing.assert_allclose(se,sec.conj().T,atol = 1e-14)
            
    #def test_gamma(self):
        
        #for tb in self.tb:
            #for e in self.energies:
            
                #tbb = e - tb
                
                #c1 = GreensFunctionCalculator(tbb)
                #c2 = BlochCalculator(tbb)
                
                #c2.calculate_bloch_states()
                #print c2.states_l.conj().T.dot(c2.gamma()).dot(c2.states_l)
                #print c2.states_l.conj().T.dot(c1.gamma()).dot(c2.states_l)
                #print c2.v_l
                #testing.assert_allclose(c1.gamma(), c2.gamma(), atol = 1e-14)
            
    #def test__(self):
        #e = 1+.3e-8j
        #d = self.tb.periodic_device()
        #gg = (e-d).gf()[:2,:2]
        #print "Device"
        #print gg
        #print (d.leads[0] - d.leads[1]).absmax()
        #print "Transmission"
        #print (e-d).transmission(0,1)
        
        #g, bm, gm = (e-self.tb).gf(direction = 'negative', bloch_matrix = True, gamma = True)
        #g2, bm2, gm2 = (e-self.tb).gf(direction = 'positive', bloch_matrix = True, gamma = True)
        
        #w,r,vel = (e-self.tb).bloch_states(return_velocities = True)
        #print "vel"
        #print vel
        #incoming = (abs(w)>1)
        #outgoing = (abs(w)<1)
        #print "T"
        #print abs(numpy.linalg.inv(r[:,outgoing]).dot(gg).dot(gm).dot(r[:,incoming]))
        #raise
        #w,r,vel = w[right_going], r[:,right_going], vel[right_going]
        #print "selected"
        #print vel
        #print "lambda"
        #print w
        #print "bm|>/|>"
        #print bm.dot(r)/r
        #print "G*Gamma"
        #print gg.dot(gm)
        #print "bm"
        #print bm2
        #print "ratio"
        #print bm2/gg.dot(gm)
        #raise
        
#class TightBinding1DLargeTest(unittest.TestCase):
    
    #def setUp(self):
        #N = 5
        #random.seed(0)
        #self.tb = TightBinding({
            #(0,): random.rand(N,N) + 1j*random.rand(N,N) - 0.5 - 0.5j,
            #(1,):  random.rand(N,N) + 1j*random.rand(N,N) - 0.5 - 0.5j,
            #(-1,):  random.rand(N,N) + 1j*random.rand(N,N) - 0.5 - 0.5j,
        #})
        #self.tb = self.tb + self.tb.hc()
            
    #def test_bloch(self):
        
        #for e in numpy.linspace(-1.1,1.1,17)+1e-6j:
            
            #tbb = e - self.tb
            #calc = BlochCalculator(tbb)
            #calc.calculate_bloch_states()
            
            #for r,w,vel in ((calc.states_l, calc.w_l, calc.v_l), (calc.states_r, calc.w_r, calc.v_r)):
                
                #k = (-1j*numpy.log(abs(w)) + numpy.angle(w))/numpy.pi/2
                #aw = abs(w)
                #testing.assert_allclose((abs(r)**2).sum(axis = 0),1)
                
                #for i in range(len(w)):
                    
                    #m = self.tb.fourier(k[i]).diagonal
                    #v = r[:,i]
                    
                    #if aw[i] == float("inf"):
                        #pass
                    
                    #elif aw[i] != aw[i]:
                        #raise ValueError("NaN encountered")
                        
                    #elif aw[i] == 0:
                        #pass
                        
                    #else:
                        #testing.assert_allclose(self.tb[-1].dot(v) + w[i]*self.tb[0].dot(v) + w[i]**2*self.tb[1].dot(v),e*v*w[i], atol = 1e-14)
                        #testing.assert_allclose(m.dot(v),e*v, atol = 1e-14)

#class MTDTest(unittest.TestCase):
    
    #def __model_tb__(self, m1, m2):
        #return TightBinding({
            #(0,): m1,
            #(1,): m2,
            #(-1,): numpy.array(m2).conj().T,
        #})
        
    #def __model__(self, m1, m2, size = 0):
        #return self.__model_tb__(m1,m2).periodic_device(size = size)
        
    #def test_gf(self):
        #e = .3j
        #d = self.__model__([[3,0],[0,4]],[[1,0],[0,2]])
        #g1 = MTDCalculator(e - d).gf()
        #g2 = MTDCalculator(e - d.add_lead(0)).gf()
        #testing.assert_allclose(g1,g2[:g1.shape[0],:g1.shape[1]])
        
    #def test_gf2(self):
        #e = 0.1+.3e-8j
        #d = self.__model__([[0,0],[0,0]],[[1,0],[0,-1]], size = 3).add_lead(0)
        #d.center.diagonal[2:4,2:4] = [[.3,0],[0,.4]]
        #g1 = MTDCalculator(e - d).gf()
        #d = d.add_lead(0)
        #g2 = MTDCalculator(e - d).gf()
        #testing.assert_allclose(g1,g2[:g1.shape[0],:g1.shape[1]])
        #d = d.add_lead(1)
        #g3 = MTDCalculator(e - d).gf()
        #testing.assert_allclose(g1,g3[:g1.shape[0],:g1.shape[1]])
        
    #def test_gf3(self):
        #tb = self.__model_tb__([[0]],[[1]])
        #d = tb.periodic_device()
        #d.center.diagonal[0,0] += 0.3
        #d2 = d.add_lead(0).add_lead(1)
        #d3 = d2.add_lead(0)
        #d5 = d3.add_lead(1)
        #for e in numpy.linspace(0,3,30)+1e-3j:
            #gf3 = MTDCalculator(e - d3).gf()
            #gf5 = MTDCalculator(e - d5).gf()
                
            #testing.assert_allclose(gf3,gf5[:gf3.shape[0],:gf3.shape[1]])
        
    #def __test_1__(self, a, b, size = 1):
        #m = self.__model__([[a]],[[b]], size = size)
        
        #for calc in (GreensFunctionCalculator, BlochCalculator):
            
            #for e in numpy.linspace(0.01,3,30)+1e-9j:
                
                #mtdc = MTDCalculator(e - m, calculator = calc)
                #actual = -numpy.imag(numpy.trace(mtdc.gf()))
                #expected = 4*abs(b)**2 - (e-a)**2
                #expected = 1/expected**.5 if expected>0 else 0
                #expected *= m.center.shape[0]
                #testing.assert_allclose(actual,expected,atol = 1e-6*size)
                
                #transmission = mtdc.transmission(0, 1)
                
                #if abs(e-a)<2*abs(b):
                    #testing.assert_allclose(transmission,1,atol = 1e-7)
                    
                #else:
                    #testing.assert_allclose(transmission,0, atol = 1e-7)
            
    #def test_1_simple(self):
        #self.__test_1__(0,1)
        
    #def test_1_diag(self):
        #self.__test_1__(0.1,1)
        
    #def test_1_arb(self):
        #self.__test_1__(0.1,0.5*(1+1j))
        
    #def test_1_2x(self):
        #self.__test_1__(0.1,0.5*(1+1j), size = 2)
        
    #def test_1_3x(self):
        #self.__test_1__(0.1,0.5*(1+1j), size = 3)
        
    #def test_2_random(self):
        #random.seed(0)
        #tb = self.__model_tb__(random.rand(2,2) + 1j*random.rand(2,2) - 0.5 - 0.5j, random.rand(2,2) + 1j*random.rand(2,2) - 0.5 - 0.5j)
        #tb = self.__model_tb__([[.1+.2j, .3-.4j], [-.5-.6j, .7+.8j]],[[.09+.1j, -.11-.12j], [-.13+.14j, .15+.16j]])
        #tb = tb + tb.hc()
        #assert tb == tb.hc()
        #d1 = tb.periodic_device()
        #assert d1.leads[0] == d1.leads[0].hc()
        #assert d1.leads[1] == d1.leads[1].hc()
        #d2 = tb.periodic_device(size = 2)
        
        #for e in numpy.linspace(0.1,3,30)+1e-8j:
            
            #gf1 = MTDCalculator(e - d1).gf()
            #gf2 = MTDCalculator(e - d2).gf()
            #testing.assert_allclose(gf1,gf2[:gf1.shape[0],:gf1.shape[1]])

    #def test_2_random_ovlp(self):
        #random.seed(0)
        #d1 = tb.periodic_device()
        #d1_o = tb_o.periodic_device()
        #d2 = tb.periodic_device(size = 2)
        #d2_o = tb_o.periodic_device(size = 2)
        #d3 = tb.periodic_device(size = 3)
        #d3_o = tb_o.periodic_device(size = 3)
        
        #numpy.set_printoptions(linewidth = 200, precision = 6)
        
        #for e in numpy.linspace(-3,3,30)+1e-2j:
            
            #gf1 = MTDCalculator(d1_o*e - d1).gf()
            #gf2 = MTDCalculator(d2_o*e - d2).gf()
            #gf3 = MTDCalculator(d3_o*e - d3).gf()
            #testing.assert_allclose(gf1,gf2[:gf1.shape[0],:gf1.shape[1]])
            #testing.assert_allclose(gf2,gf3[:gf2.shape[0],:gf2.shape[1]])
            
    #def test_transmission_reordered(self):
        #random.seed(0)
        #N = 5
        #tb = self.__model_tb__(random.rand(N,N) + 1j*random.rand(N,N) - 0.5 - 0.5j, random.rand(N,N) + 1j*random.rand(N,N) - 0.5 - 0.5j)
        #tb = tb + tb.hc()
        
        #d1 = tb.periodic_device()
        
        #d1_r = tb.periodic_device()
        
        #new_order = random.permutation(d1_r.center.shape[0])
        #d1_r.center = d1_r.center.foreach(lambda k,v: (k, v[new_order,:][:,new_order]))
        #for i in range(len(d1_r.connections)):
            #d1_r.connections[i] = d1_r.connections[i][:,new_order]
        
        #for e in numpy.linspace(0.01,3,30)+1e-8j:
            #t1 = MTDCalculator(e - d1).transmission(0,1)
            #t2 = MTDCalculator(e - d1_r).transmission(0,1)
            #testing.assert_allclose(abs(t1),round(abs(t1)),atol = 1e-6)
            #testing.assert_allclose(t1,t2)

    #def test_transmission_supercell(self):
        #tb = self.__model_tb__([[0]],[[1]])
        #d = tb.periodic_device()
        #d.center.diagonal[0,0] += 0.3
        #d1 = d.add_lead(0).add_lead(1)
        #d2 = d1.add_lead(1)
        #for e in numpy.linspace(0.01,3,30)+1e-8j:
            #t1 = MTDCalculator(e - d1).transmission(0,1)
            #t2 = MTDCalculator(e - d2).transmission(0,1)
            #testing.assert_allclose(t1,t2, atol = 1e-8)
            
    #def test_transmission_matrix_triv(self):
        #tb = self.__model_tb__([[0]],[[1]])
        #d = tb.periodic_device()
        #for e in numpy.linspace(0.01,3,30)+1e-8j:
            #t1 = MTDCalculator(e - d).transmission(0,1)
            #t2 = MTDCalculator(e - d, calculator = BlochCalculator).tmatrix(0,1)
            #testing.assert_allclose(t1,numpy.trace(t2.dot(t2.conj().T)), atol = 1e-14)
      
    #def test_transmission_matrix_defect(self):
        #tb = self.__model_tb__([[0]],[[1]])
        #d = tb.periodic_device()
        #d.center.diagonal[0,0] += 0.3
        #d1 = d.add_lead(0).add_lead(1)
        #for e in numpy.linspace(0.01,3,30)+1e-8j:
            #t1 = MTDCalculator(e - d1).transmission(0,1)
            #t2 = MTDCalculator(e - d1, calculator = BlochCalculator).tmatrix(0,1)
            #testing.assert_allclose(t1,numpy.trace(t2.dot(t2.conj().T)), atol = 1e-14)
            
if __name__ == '__main__':
    unittest.main()
