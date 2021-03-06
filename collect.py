#!/bin/env python
import collections

d = {}

with open('hh1.dat') as f:
    for line in f:
        xs = [float(x) for x in line.split(' ') if x]
        d[(round(xs[1]), round(xs[0]))] = (xs[2])

with open('hh2.dat') as f:
    for line in f:
        xs = [float(x) for x in line.split(' ') if x]
        d[(round(xs[1]), round(xs[0]))] = (d[(round(xs[1]), round(xs[0]))],  xs[2])

with open('f_1sq.dat') as f:
    for line in f:
        xs = [float(x) for x in line.split(' ') if x]
        d[(round(xs[1]), round(xs[0]))] = (d[(round(xs[1]), round(xs[0]))],  xs[2])

with open('f_2sq_cl.dat') as f:
    for line in f:
        xs = [float(x) for x in line.split(' ') if x]
        d[(round(xs[1]), round(xs[0]))] = (d[(round(xs[1]), round(xs[0]))],  xs[2])

with open('data.txt', 'w') as f:
    od = collections.OrderedDict(sorted(d.items()))
    for key in od:
        (y, x) = key
        (((h1, h2), f1), f2) = d[key]
        print(x, y, h1, h2, f1, f2, file=f)
