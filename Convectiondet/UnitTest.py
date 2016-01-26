# -*- coding: utf-8 -*-
"""
Created on Fri May 16 12:53:26 2014

@author: Vlad
"""
import unittest

class FooTest(unittest.TestCase):
    
    def testFoo(self):
        self.failUnless(False)

def main():
    unittest.main()

if __name__ == '__main__':
    main()
