import numpy as npy

def f1(x):
    """
    Test function 1
    """
    return x*x*x - npy.pi*x + npy.e/100

def f2(x):
    """
    Test function 2
    """
    return -1.13 + npy.tanh(x-2) + 4*npy.exp(-x)*npy.sin((1/8.)*x**3) \
           *x + .1*npy.exp((1/35.)*x)