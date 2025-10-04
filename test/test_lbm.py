import unittest
from lbm import *
import numpy as np

class TestStringMethods(unittest.TestCase):
    ## for template / exmample only ##
    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

class TestStream(unittest.TestCase):
    def __init__():
        self.test = 0 
        # TODO
       

if __name__ == '__main__':
    unittest.main()