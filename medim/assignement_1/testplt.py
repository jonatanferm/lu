import cv2
import numpy as np


A = np.array([[1, 3, 5, 7],[1, 3, 5, 7]])
B = np.array([[2, 4, 6, 8],[2, 4, 6, 8]])

C = np.vstack((A, B))

A[:,::2] = 0
print(A)
B[:,1::2] = 0
print(B)

print(A + B)


