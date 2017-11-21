import numpy as np 


def main():
    N = 51
    X = np.linspace(-50, 50, N)
    Y = np.linspace(-50, 50, N)
    with open('H1.dat', 'w') as f:
        with open('H2.dat', 'w') as f2:
            for x in X:
                for y in Y:
                    print('{0} {1} {2:.1f}'.format(x, y, 200 * np.e ** (-(np.abs(x) + np.abs(y)) ** 0.5) + 1000), file=f)
                    print('{0} {1} {2:.1f}'.format(x, y, 100 * np.e ** (-(np.abs(x) + np.abs(y)) ** 0.45) + 1500), file=f2)



main()




