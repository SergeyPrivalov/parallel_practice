#!/bin/env python

d = {}


with open('hh1.dat') as f:
    for line in f:
        xs = [float(x) for x in line.split(' ') if x]
        d[(round(xs[0]), round(xs[1]))] = (xs[2])
        
with open('hh2.dat') as f:
    for line in f:
        xs = [float(x) for x in line.split(' ') if x]
        d[(round(xs[0]), round(xs[1]))] = (d[(round(xs[0]), round(xs[1]))],  xs[2])

with open('f_1sq.dat') as f:
    for line in f:
        xs = [float(x) for x in line.split(' ') if x]
        d[(round(xs[0]), round(xs[1]))] = (d[(round(xs[0]), round(xs[1]))],  xs[2])

with open('f_2sq_cl.dat') as f:
    for line in f:
        xs = [float(x) for x in line.split(' ') if x]
        d[(round(xs[0]), round(xs[1]))] = (d[(round(xs[0]), round(xs[1]))],  xs[2])
        
with open('data.txt', 'w') as f:
    for key, val in d.items():
        (x, y) = key
        (((h1, h2), f1), f2) = val
        print(x, y, h1, h2, f1, f2, file=f)
