import numpy as np 


def main():
    N = 51
    X = np.linspace(-50, 50, N)
    Y = np.linspace(-50, 50, N)
    
    H1 = {}
    H2 = {}
    SIGMA = {}
    F = {}
    A = {}

    with open('H1.dat', 'w') as f:
        with open('H2.dat', 'w') as f2:
            with open('sigma.dat', 'w') as f3:
                for x in X:
                    for y in Y:
                        H1[(x,y)] = 200 * np.e ** (-(np.abs(x) + np.abs(y))/75) + 1000
                        print('{0} {1} {2:.1f}'.format(x, y, H1[(x,y)]), file=f)

                        H2[(x,y)] = 100 * np.e ** (-(np.abs(x) + np.abs(y))/75) + 1500
                        print('{0} {1} {2:.1f}'.format(x, y, H2[(x,y)]), file=f2)

                        SIGMA[(x,y)] = 20 if x == 0 and y == 0 else 50.0 / np.sqrt(x*x + y*y) * np.sin((x + y)/10)
                        print('{0} {1} {2:.1f}'.format(x, y, SIGMA[(x,y)]), file=f3)


    f = 0.00667
    with open('F.dat', 'w') as f4:  
        for xk in X:
            for yl in Y:
                F = 0
                for xi in X:
                    for yj in Y:
                        F += (
                                    ((xk - xi) ** 2 + (yl - yj) ** 2 + H1[(xi, yj)] ** 2) ** -0.5 -
                                    ((xk - xi) ** 2 + (yl - yj) ** 2 + H2[(xi, yj)] ** 2) ** -0.5
                                 ) * SIGMA[(xi,yj)] * (X[1] - X[0]) ** 2
                print('{0} {1} {2:.10f}'.format(xk, yl, F), file=f4)



main()

