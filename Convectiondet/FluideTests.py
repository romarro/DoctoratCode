import unittest
from __init__ import um, Q_
import pint
from Fluide import *

class Test_FluideTests(unittest.TestCase):
    
    def assertAlmostEquals(self, first, second, places = None, msg = None, delta = None):
        if not isinstance(first,type(second)):
            return self.fail("first and the second doen't have the same type")
        if isinstance(first, np.ndarray):
            if len(first)!=len(second):
                return self.fail("fist length is different than second's length")
            for f,s in zip(first,second):
                self.assertAlmostEquals(f,s,places,msg,delta)    
            return
        if isinstance(first,um.Quantity):
            if first.units!=second.units:
                try:
                    first=first.to(second.units)
                except pint.DimensionalityError:
                    return self.fail("first and  second have incompatible units")

            return self.assertAlmostEquals(first.magnitude,second.magnitude,places,msg,delta)

        if isinstance(first,un.AffineScalarFunc):
            diff=abs(first-second)
            if places:
                return self.assertTrue(round(diff.n,places)<round(diff.s,places),"the nominal value of the difference is greater than the error. first: %s, second:%s"%(first,second))
            elif delta:
                return self.assertTrue(diff<=delta,"the diff is greater than delta first: %s, second:%s"%(first,second))     
            else:
                return self.assertTrue(diff.n<diff.s,"diff nominal value greater than error first: %s, second:%s"%(first,second))    
        
        return super(Test_FluideTests, self).assertAlmostEquals(first, second, places, msg, delta)

    def test_assertAlmosEquals(self):
        self.assertAlmostEquals(2.,2.001,2)
        self.assertAlmostEquals(Q_(20.2,um.degC),Q_(20.3,um.degC),delta=0.5)
        self.assertAlmostEquals(un.ufloat(2.4,0.2),un.ufloat(2.38,0.4),1)
        self.assertAlmostEquals(un.ufloat(2.4,0.2),un.ufloat(2.38,0.4),delta=0.3)

        sir1 = np.array([un.ufloat(999.7496259274476,0.0017598342895507812),un.ufloat(998.2523477831393,0.004131866455078125)])
        sir2 = np.array([un.ufloat(999.749,0.0017598342895507812),un.ufloat(998.252347,0.004131866455078125)])
        self.assertAlmostEquals(sir1,sir2,4)

    def test_Water_dens(self):
        
        water = Water(0)
        temps=np.array([un.ufloat(10,0.02),un.ufloat(20,0.02)])
        densst=np.array([un.ufloat(999.7496259274476,0.0017598342895507812),un.ufloat(998.2523477831393,0.004131866455078125)])*um.kg/um.m**3
        dens=water.densit(Q_(temps,um.degC))
        self.assertAlmostEquals(densst,dens,3)

    def test_Water_cp(self):
        water = Water()
        temps=Q_(np.array([un.ufloat(10,0.02),un.ufloat(20,0.02)]),um.degC)
        cp=np.array([un.ufloat(4196.46185279,0.0),  un.ufloat(4184.46877474,0.0)])*um.J/um.kg/um.K
        self.assertAlmostEquals(cp,water.cp(temps),delta=5.0)


    def test_Water_Pr(self):

        water =Water()
        temps=np.array([un.ufloat(10,0.02),un.ufloat(20,0.02)])
        pr=np.array([un.ufloat(9.463021047858366,0.006206084509439321),un.ufloat(7.006353702684516,0.0038692636351550565)])*um.dimensionless
        self.assertAlmostEquals(pr,water.Pr(Q_(temps,um.degC)),3)

if __name__ == '__main__':
    unittest.main()
