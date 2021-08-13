## pool.py
from multiprocessing import Pool, cpu_count
import sys
import os

def square(x):
    """Function to return the square of the argument"""
    return x*x

if __name__ == "__main__":
    # print the number of cores
    print("Number of cores available equals %d" % cpu_count())

    # create a pool of workers
    pool = Pool()

    # create an array of 5000 integers, from 1 to 5000
    a = range(1,5001)

    result = pool.map( square, a )

    total = reduce( lambda x,y: x+y, result )

    print("The sum of the square of the first 5000 integers is %d" % total)
