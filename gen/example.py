import numpy as np 


def main():
    N = 51
    X = np.linspace(-50, 50, N)
    Y = np.linspace(-50, 50, N)
    
    with open('H1.dat', 'w') as f:
        with open('H2.dat', 'w') as f2:
        	with open('sigma.dat', 'w') as f3:
	            for x in X:
	                for y in Y:
	                    print('{0} {1} {2:.1f}'.format(x, y, 200 * np.e ** (-(np.abs(x) + np.abs(y))/75) + 1000), file=f)
	                    print('{0} {1} {2:.1f}'.format(x, y, 100 * np.e ** (-(np.abs(x) + np.abs(y))/75) + 1500), file=f2)
	                    print('{0} {1} {2:.1f}'.format(x, y, 50 if x == 0 and y == 0 else 50.0 / np.sqrt(x*x + y*y) * np.sin(x + y)), file=f3)


main()




