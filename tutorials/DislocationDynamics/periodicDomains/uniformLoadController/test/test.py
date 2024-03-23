import numpy as np

originalSFEmean = 300
a = np.loadtxt('data.txt')+originalSFEmean
print(f'average = {np.mean(a):.2f}')
