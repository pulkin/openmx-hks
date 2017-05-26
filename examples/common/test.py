#!/usr/bin/env python

# ======================================================================
# This file is under construction.
# ======================================================================

import unittest

import numpy
from numpy import testing, random

from tb import *

class TightBinding1DDiagTests(unittest.TestCase):
    
    def setUp(self):
        self.r1 = TightBinding({
            (0,): [[0]],
            (1,): [[1]],
            (-1,): [[1]],
        })
        
        self.c1 = TightBinding({
            (0,): [[0]],
            (1,): [[1j]],
            (-1,): [[-1j]],
        })
        
        self.r2 = TightBinding({
            (0,): [[0,0],[0,0]],
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
        
    def test_eq(self):
        for i in self.a:
            assert i.copy() == i
            assert (-i) == (-i)
            assert i.hc() == i
            assert 2*i == i + i
            assert 0.5*i == i - 0.5*i
            assert 0.5*i == i/2
            assert i.is_1DNN()
            assert i + 3.14*i.eye() == i + 3.14
            
    def test_eigpath(self):
        for i in self.a:
            p = numpy.linspace(-.5,.5,31)
            eig = i.eig_path(p[:,numpy.newaxis])
            a = numpy.diag(i[0])[numpy.newaxis,:]
            b = numpy.diag(i[1])[numpy.newaxis,:]
            testing.assert_allclose(eig,
                numpy.sort(a + 2*numpy.cos(2*numpy.pi*p)[:,numpy.newaxis]*b.real - 2*numpy.sin(2*numpy.pi*p)[:,numpy.newaxis]*b.imag, axis = 1)
            )
            
    def test_herm_shift(self):
        for i in self.a:
            p = numpy.linspace(-.5,.5,31)
            eig = i.eig_path(p[:,numpy.newaxis])
            mask = numpy.zeros(i.shape[0], dtype = numpy.bool)
            mask[::2] = True
            i.hermitian_shift_subblock(mask, (3,))
            eig2 = i.eig_path(p[:,numpy.newaxis])
            testing.assert_allclose(eig, eig2)

    def test_super(self):
        for i in self.a:
            s = i.super(2,0)
            assert i == i.hc()
            assert len(i.__m__) == 3
            eig = i.eig_path(numpy.linspace(0,1,30,endpoint = False)[:,numpy.newaxis])
            eig2 = s.eig_path(numpy.linspace(0,1,15,endpoint = False)[:,numpy.newaxis])
            testing.assert_allclose(numpy.sort(numpy.concatenate((eig[:15],eig[15:]),axis = 1),axis = 1),eig2, atol = 1e-14)
            
    def test_gf(self):
        for i in self.a:
            for e in numpy.linspace(-1,1,30)+1e-1j:
                w = e-i
                gfc = GreensFunctionCalculator(w)
                g = gfc.gf()
                se = gfc.self_energy()
                testing.assert_allclose((w[0] - se).dot(g), numpy.eye(g.shape[0]))
                
    #def test_gamma(self):
        #for i in self.a:
            #for e in numpy.linspace(-1,1,30)+1e-3j:
                #g, gm = (e-i).gf(gamma_ = True)
                #w,r,vel = (e-i).bloch_states(return_velocities = True)
                #s = abs(w)>1
                #testing.assert_allclose(vel[s],numpy.diag(r[:,s].conj().T.dot(gm).dot(r[:,s])))
            
class LeadCalculatorTest(unittest.TestCase):
    
    def setUp(self):
        
        self.tb = [
            TightBinding({
                (0,): [[0,1+.2j],[1-.2j,0]],
                (1,): [[0,0.5],[.1-.2j,0.3]],
                (-1,): [[0,.1+.2j],[0.5,0.3]],
            }),
            TightBinding({
                (0,): [[0,1],[1,0]],
                (1,): [[0,0],[.1,0]],
                (-1,): [[0,.1],[0,0]],
            }),
        ]
        
        self.energies = [.1j, 1+.1j, 1+1e-8j]
            
    def test_bloch(self):
        
        for tb in self.tb:
            for e in self.energies:
            
                tbb = e - tb
                calc = BlochCalculator(tbb)
                calc.calculate_bloch_states()
                
                for r,w,vel in ((calc.states_l, calc.w_l, calc.v_l), (calc.states_r, calc.w_r, calc.v_r)):
                    
                    aw = abs(w)
                    testing.assert_allclose((abs(r)**2).sum(axis = 0),1)
                    
                    for i in range(len(w)):
                        
                        if aw[i] == float("inf"):
                            pass
                        
                        elif aw[i] != aw[i]:
                            raise ValueError("NaN encountered")
                            
                        elif aw[i] == 0:
                            pass
                            
                        else:
                            
                            k = (-1j*numpy.log(abs(w)) + numpy.angle(w))/numpy.pi/2
                            m = tb.fourier(k[i]).diagonal
                            v = r[:,i]
                        
                            testing.assert_allclose(tb[-1].dot(v) + w[i]*tb[0].dot(v) + w[i]**2*tb[1].dot(v),e*v*w[i], atol = 1e-14)
                            testing.assert_allclose(m.dot(v),e*v, atol = 1e-14)
                    
    def test_gf(self):
        
        for tb in self.tb:
            for e in self.energies:
            
                g1 = GreensFunctionCalculator(e - tb).gf()
                g2 = GreensFunctionCalculator(e - tb.super(2,0)).gf()
                testing.assert_allclose(g1,g2[2:,2:])
                w = e - tb
                testing.assert_allclose(g1, numpy.linalg.inv(w[0] - w[-1].dot(g1).dot(w[1])))
        
    def test_bloch_matrix(self):
        
        for tb in self.tb:
            for e in self.energies:
            
                tbb = e - tb
                calc = BlochCalculator(tbb)
                calc2 = GreensFunctionCalculator(tbb)
                
                testing.assert_allclose(calc.bloch_matrix("-",-1), calc2.bloch_matrix(), atol = 1e-14)
                
                testing.assert_allclose(calc.bloch_matrix("-",0), numpy.eye(2), atol = 1e-14)
                testing.assert_allclose(calc.bloch_matrix("+",0), numpy.eye(2), atol = 1e-14)
                
                testing.assert_allclose(tbb[-1].dot(calc.bloch_matrix("-",-2)) + tbb[0].dot(calc.bloch_matrix("-",-1)),-tbb[1], atol = 1e-14)
                testing.assert_allclose(tbb[1].dot(calc.bloch_matrix("+",2)) + tbb[0].dot(calc.bloch_matrix("+",1)),-tbb[-1], atol = 1e-14)
            
    def test_self_energy(self):
            
        for tb in self.tb:
            for e in self.energies:
            
                tbb = e - tb

                c1 = GreensFunctionCalculator(tbb)
                se1 = c1.self_energy()
                
                c2 = BlochCalculator(tbb)
                se2 = c2.self_energy()
                
                testing.assert_allclose(se1,se2, atol = 1e-14)
                
    def self_energy_conj_test(self):
        
        for tb in self.tb:
            for e in numpy.linspace(1,0,3)+1e-8j:
                
                tbb = e - tb
                c = BlochCalculator(tbb)
                se = c.self_energy().dot(c.states_r).dot(numpy.diag(abs(abs(c.w_r)-1)>1e-5)).dot(c.states_ri)
                sec = c.self_energy(conj = True).dot(c.states_r).dot(numpy.diag(abs(abs(c.w_r)-1)>1e-5)).dot(c.states_ri)
                
                testing.assert_allclose(se,sec.conj().T,atol = 1e-14)
            
    def test_gamma(self):
        
        for tb in self.tb:
            for e in self.energies:
            
                tbb = e - tb
                
                c1 = GreensFunctionCalculator(tbb)
                c2 = BlochCalculator(tbb)
                
                c2.calculate_bloch_states()
                print c2.states_l.conj().T.dot(c2.gamma()).dot(c2.states_l)
                print c2.states_l.conj().T.dot(c1.gamma()).dot(c2.states_l)
                print c2.v_l
                testing.assert_allclose(c1.gamma(), c2.gamma(), atol = 1e-14)
            
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

class MTDTest(unittest.TestCase):
    
    def __model_tb__(self, m1, m2):
        return TightBinding({
            (0,): m1,
            (1,): m2,
            (-1,): numpy.array(m2).conj().T,
        })
        
    def __model__(self, m1, m2, size = 0):
        return self.__model_tb__(m1,m2).periodic_device(size = size)
        
    def test_gf(self):
        e = .3j
        d = self.__model__([[3,0],[0,4]],[[1,0],[0,2]])
        g1 = MTDCalculator(e - d).gf()
        g2 = MTDCalculator(e - d.add_lead(0)).gf()
        testing.assert_allclose(g1,g2[:g1.shape[0],:g1.shape[1]])
        
    def test_gf2(self):
        e = 0.1+.3e-8j
        d = self.__model__([[0,0],[0,0]],[[1,0],[0,-1]], size = 3).add_lead(0)
        d.center.diagonal[2:4,2:4] = [[.3,0],[0,.4]]
        g1 = MTDCalculator(e - d).gf()
        d = d.add_lead(0)
        g2 = MTDCalculator(e - d).gf()
        testing.assert_allclose(g1,g2[:g1.shape[0],:g1.shape[1]])
        d = d.add_lead(1)
        g3 = MTDCalculator(e - d).gf()
        testing.assert_allclose(g1,g3[:g1.shape[0],:g1.shape[1]])
        
    def test_gf3(self):
        tb = self.__model_tb__([[0]],[[1]])
        d = tb.periodic_device()
        d.center.diagonal[0,0] += 0.3
        d2 = d.add_lead(0).add_lead(1)
        d3 = d2.add_lead(0)
        d5 = d3.add_lead(1)
        for e in numpy.linspace(0,3,30)+1e-3j:
            gf3 = MTDCalculator(e - d3).gf()
            gf5 = MTDCalculator(e - d5).gf()
                
            testing.assert_allclose(gf3,gf5[:gf3.shape[0],:gf3.shape[1]])
        
    def __test_1__(self, a, b, size = 1):
        m = self.__model__([[a]],[[b]], size = size)
        
        for calc in (GreensFunctionCalculator, BlochCalculator):
            
            for e in numpy.linspace(0.01,3,30)+1e-9j:
                
                mtdc = MTDCalculator(e - m, calculator = calc)
                actual = -numpy.imag(numpy.trace(mtdc.gf()))
                expected = 4*abs(b)**2 - (e-a)**2
                expected = 1/expected**.5 if expected>0 else 0
                expected *= m.center.shape[0]
                testing.assert_allclose(actual,expected,atol = 1e-6*size)
                
                transmission = mtdc.transmission(0, 1)
                
                if abs(e-a)<2*abs(b):
                    testing.assert_allclose(transmission,1,atol = 1e-7)
                    
                else:
                    testing.assert_allclose(transmission,0, atol = 1e-7)
            
    def test_1_simple(self):
        self.__test_1__(0,1)
        
    def test_1_diag(self):
        self.__test_1__(0.1,1)
        
    def test_1_arb(self):
        self.__test_1__(0.1,0.5*(1+1j))
        
    def test_1_2x(self):
        self.__test_1__(0.1,0.5*(1+1j), size = 2)
        
    def test_1_3x(self):
        self.__test_1__(0.1,0.5*(1+1j), size = 3)
        
    def test_2_random(self):
        random.seed(0)
        tb = self.__model_tb__(random.rand(2,2) + 1j*random.rand(2,2) - 0.5 - 0.5j, random.rand(2,2) + 1j*random.rand(2,2) - 0.5 - 0.5j)
        tb = self.__model_tb__([[.1+.2j, .3-.4j], [-.5-.6j, .7+.8j]],[[.09+.1j, -.11-.12j], [-.13+.14j, .15+.16j]])
        tb = tb + tb.hc()
        assert tb == tb.hc()
        d1 = tb.periodic_device()
        assert d1.leads[0] == d1.leads[0].hc()
        assert d1.leads[1] == d1.leads[1].hc()
        d2 = tb.periodic_device(size = 2)
        
        for e in numpy.linspace(0.1,3,30)+1e-8j:
            
            gf1 = MTDCalculator(e - d1).gf()
            gf2 = MTDCalculator(e - d2).gf()
            testing.assert_allclose(gf1,gf2[:gf1.shape[0],:gf1.shape[1]])

    def test_2_random_ovlp(self):
        random.seed(0)
        d1 = tb.periodic_device()
        d1_o = tb_o.periodic_device()
        d2 = tb.periodic_device(size = 2)
        d2_o = tb_o.periodic_device(size = 2)
        d3 = tb.periodic_device(size = 3)
        d3_o = tb_o.periodic_device(size = 3)
        
        numpy.set_printoptions(linewidth = 200, precision = 6)
        
        for e in numpy.linspace(-3,3,30)+1e-2j:
            
            gf1 = MTDCalculator(d1_o*e - d1).gf()
            gf2 = MTDCalculator(d2_o*e - d2).gf()
            gf3 = MTDCalculator(d3_o*e - d3).gf()
            testing.assert_allclose(gf1,gf2[:gf1.shape[0],:gf1.shape[1]])
            testing.assert_allclose(gf2,gf3[:gf2.shape[0],:gf2.shape[1]])
            
    def test_transmission_reordered(self):
        random.seed(0)
        N = 5
        tb = self.__model_tb__(random.rand(N,N) + 1j*random.rand(N,N) - 0.5 - 0.5j, random.rand(N,N) + 1j*random.rand(N,N) - 0.5 - 0.5j)
        tb = tb + tb.hc()
        
        d1 = tb.periodic_device()
        
        d1_r = tb.periodic_device()
        
        new_order = random.permutation(d1_r.center.shape[0])
        d1_r.center = d1_r.center.foreach(lambda k,v: (k, v[new_order,:][:,new_order]))
        for i in range(len(d1_r.connections)):
            d1_r.connections[i] = d1_r.connections[i][:,new_order]
        
        for e in numpy.linspace(0.01,3,30)+1e-8j:
            t1 = MTDCalculator(e - d1).transmission(0,1)
            t2 = MTDCalculator(e - d1_r).transmission(0,1)
            testing.assert_allclose(abs(t1),round(abs(t1)),atol = 1e-6)
            testing.assert_allclose(t1,t2)

    def test_transmission_supercell(self):
        tb = self.__model_tb__([[0]],[[1]])
        d = tb.periodic_device()
        d.center.diagonal[0,0] += 0.3
        d1 = d.add_lead(0).add_lead(1)
        d2 = d1.add_lead(1)
        for e in numpy.linspace(0.01,3,30)+1e-8j:
            t1 = MTDCalculator(e - d1).transmission(0,1)
            t2 = MTDCalculator(e - d2).transmission(0,1)
            testing.assert_allclose(t1,t2, atol = 1e-8)
            
    def test_transmission_matrix_triv(self):
        tb = self.__model_tb__([[0]],[[1]])
        d = tb.periodic_device()
        for e in numpy.linspace(0.01,3,30)+1e-8j:
            t1 = MTDCalculator(e - d).transmission(0,1)
            t2 = MTDCalculator(e - d, calculator = BlochCalculator).tmatrix(0,1)
            testing.assert_allclose(t1,numpy.trace(t2.dot(t2.conj().T)), atol = 1e-14)
      
    def test_transmission_matrix_defect(self):
        tb = self.__model_tb__([[0]],[[1]])
        d = tb.periodic_device()
        d.center.diagonal[0,0] += 0.3
        d1 = d.add_lead(0).add_lead(1)
        for e in numpy.linspace(0.01,3,30)+1e-8j:
            t1 = MTDCalculator(e - d1).transmission(0,1)
            t2 = MTDCalculator(e - d1, calculator = BlochCalculator).tmatrix(0,1)
            testing.assert_allclose(t1,numpy.trace(t2.dot(t2.conj().T)), atol = 1e-14)
            
if __name__ == '__main__':
    unittest.main()
