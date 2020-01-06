import numpy as np
a = np.array([1,3,5,7])
b = np.array([1,2,3,4,5,6,7,8,9,10])
for i,j in zip(a, np.append(a[1:], None)):
    print(i)
    print(j)
    print(b[i:j])