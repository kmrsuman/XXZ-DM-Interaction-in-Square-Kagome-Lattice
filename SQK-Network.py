from turtle import width
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from pandas import pivot
from matplotlib.colors import Normalize

G=nx.Graph()


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


for i in range(N):
    G.add_node(i, pos = (x_pos[i], y_pos[i]))


#nx.read_edgelist("sqk_edges.edgelist")

G.add_edge(0,4)
G.add_edge(0,3)
G.add_edge(3,64)
G.add_edge(60,64)
G.add_edge(60,63)
G.add_edge(63,124)
G.add_edge(120,124)
G.add_edge(120,123)
G.add_edge(123,184)
G.add_edge(180,184)
G.add_edge(180,183)
G.add_edge(183,244)
G.add_edge(240,244)
G.add_edge(240,243)
G.add_edge(243,304)
G.add_edge(300,304)
G.add_edge(300,303)
G.add_edge(303,364)
G.add_edge(360,364)
G.add_edge(360,363)
G.add_edge(363,424)
G.add_edge(420,424)
G.add_edge(420,423)
G.add_edge(423,484)
G.add_edge(480,484)
G.add_edge(480,483)
G.add_edge(483,544)
G.add_edge(540,544)
G.add_edge(540,543)
G.add_edge(542,543)
G.add_edge(542,545)
G.add_edge(545,549)
G.add_edge(549,548)
G.add_edge(546,549)
G.add_edge(548,551)
G.add_edge(551,555)
G.add_edge(555,554)
G.add_edge(554,557)
G.add_edge(557,561)
G.add_edge(561,560)
G.add_edge(560,563)
G.add_edge(563,567)
G.add_edge(567,566)
G.add_edge(566,569)
G.add_edge(569,573)
G.add_edge(573,572)
G.add_edge(572,575)
G.add_edge(575,579)
G.add_edge(579,578)
G.add_edge(578,581)
G.add_edge(581,585)
G.add_edge(585,584)
G.add_edge(584,587)
G.add_edge(587,591)
G.add_edge(591,590)
G.add_edge(590,593)
G.add_edge(593,597)
G.add_edge(597,596)
G.add_edge(596,599)






for i in range(N):
    if i!=60 and i !=120 and i!=180 and i!= 240 and i!= 300 and i!=360 and i!= 420 and i!=480 and i!= 540 and i != 0:

        if i%6 == 0:
            a = (i+1)%N
            b = (i-1+N)%N
            c = (i+3)%N
            d = (i+4)%N
            G.add_edge(i,a)
            G.add_edge(i,b)
            G.add_edge(i,c)
            G.add_edge(i,d)
            
        elif i%6 == 1:
            a = (i+1)%N
            b = (i+3)%N
            c = (i+4)%N
            d = (i-1+N)%N
            G.add_edge(i,a)
            G.add_edge(i,b)
            G.add_edge(i,c)
            G.add_edge(i,d)

    if i!=542 and i!=548 and i!=554 and i!=560 and i!=566 and i!=572 and i!=578 and i!=584 and i!=590 and i!=596:

        if i%6 == 2:
            a = (i+1)%N
            b = (i+3)%N
            c = (i-1+N)%N
            d = (i + 6*L +2)%N
            G.add_edge(i,a)
            G.add_edge(i,b)
            G.add_edge(i,c)
            G.add_edge(i,d)

    if i!=3 and i!=63 and i!=123 and i!=183 and i!= 243 and i!=303 and i!=363 and i!=423 and i!=483 and i!= 543 :
        if i!=549 and i!=555 and i!=561 and i!=567 and i!= 573 and i!=579 and i!=585 and i!=591 and i!=597 :

            if i%6 == 3:
                a = (i - 1 + N)%N
                b = (i - 3 + N)%N
                c = (i - 4 + N)%N
                d = (i + 6*L +1)%N
                G.add_edge(i,a)
                G.add_edge(i,b)
                G.add_edge(i,c)
                G.add_edge(i,d)
        

plt.rcParams["figure.figsize"] = [8.00, 8.00]
plt.rcParams["figure.autolayout"] = True

'''
options = {
    "font_size":10,
    "node_size": 3,
    "node_color": "black",
    "edgecolors": "black",
    "linewidths": 5,
    "width": 2,
}
'''
pos=nx.get_node_attributes(G,'pos')

fig, ax = plt.subplots(figsize=(10,10))
nx.draw(G,pos, node_size = 8, node_color = "k", style='--', linewidths=1) #, **options
#ax.margins(0.10)


########################################################################

with open("square_kagome_anglesD2.5Del2.5.txt") as g :
    i = 0
    
    angles = np.zeros((N, 2))
    for line in g:
        x = line.split()
        for j in range(2):
            angles[i][j] = x[j]
        i = i + 1
g.close

Sx = []
Sy = []
Sz = []

for i in range(N):
    
    Sx = np.append(Sx, np.cos(angles[i][1]))
    Sy = np.append(Sy, np.sin(angles[i][1]))
    Sz = np.append(Sz, np.cos(angles[i][0]))

color = Sz

NORMZ = Normalize(vmin=-1.0, vmax=1.0)


diagram =ax.quiver(x_pos,y_pos, Sx, Sy, color,norm=NORMZ,cmap='brg',pivot ='mid',scale = 30,headwidth=5,headlength=5)
fig.colorbar(diagram)
ax.scatter(x_pos,y_pos, color = 'k', s = 2)
plt.title("SQK Spins delta = 2.5 D/J = -2.5",fontsize = 30)


#####################################################################


plt.axis("off")
plt.savefig("SQ_kagome_spins.jpg", dpi = 500)
plt.show()