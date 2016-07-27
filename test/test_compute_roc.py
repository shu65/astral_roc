'''
Created on 2016/03/05

@author: shu
'''
import unittest
from src.compute_roc import is_true_positive, is_false_positive, get_key


class TestROCComputer(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_is_true_positive_true(self):
        classification1 = {"class": "a", "fold":"1", "superfamily":"1", "family":"1"}
        classification2 = {"class": "a", "fold":"1", "superfamily":"1", "family":"1"}
        self.assertTrue(is_true_positive(classification1, classification2))
        
    def test_is_true_positive_false(self):
        classification1 = {"class": "a", "fold":"1", "superfamily":"1", "family":"1"}
        classification2 = {"class": "a", "fold":"1", "superfamily":"2", "family":"1"}
        self.assertFalse(is_true_positive(classification1, classification2))
        
    def test_is_false_positive_diff_class(self):
        classification1 = {"class": "a", "fold":"1", "superfamily":"1", "family":"1"}
        classification2 = {"class": "b", "fold":"1", "superfamily":"1", "family":"1"}
        self.assertTrue(is_false_positive(classification1, classification2))
        
    def test_is_false_positive_diff_fold(self):
        classification1 = {"class": "a", "fold":"1", "superfamily":"1", "family":"1"}
        classification2 = {"class": "a", "fold":"2", "superfamily":"1", "family":"1"}
        self.assertTrue(is_false_positive(classification1, classification2))
        
                
    def test_is_false_positive_diff_superfamily(self):
        classification1 = {"class": "a", "fold":"1", "superfamily":"1", "family":"1"}
        classification2 = {"class": "a", "fold":"1", "superfamily":"2", "family":"1"}
        self.assertFalse(is_false_positive(classification1, classification2))
        
        
    def test_get_key(self):
        classification = {"class": "a", "fold":"1", "superfamily":"1", "family":"1"}
        self.assertEqual(("a", "1", "1"), get_key(classification))
                

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()