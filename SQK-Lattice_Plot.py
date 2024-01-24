from locale import normalize
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
from pandas import pivot

L = 10          # REMEMBER TO ENTER THE VALUE OF L 
N = 6*L*L      # The number of sites

with open("square_kagome_positions.txt") as f :
    i = 0
    
    pos_coordinates = np.zeros((N, 2))
    for line in f:
        x = line.split()
        for j in range(2):
            pos_coordinates[i][j] = x[j]
        i = i + 1
f.close


x_pos = []
y_pos = []

for i in range(N):
    
    x_pos = np.append(x_pos,pos_coordinates[i][0])
    y_pos = np.append(y_pos,pos_coordinates[i][1])


fig, ax = plt.subplots()

ax.scatter(x_pos,y_pos, color = 'r')

plt.show()