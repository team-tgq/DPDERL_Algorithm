import numpy

a = numpy.zeros((2, 3))


def s(a):
    a[1, 1] = 1
    return a


a = s(a)
print(a)
