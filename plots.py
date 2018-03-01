#%%
%matplotlib inline

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

statsA = {'A1': 979515, 'A2': 1458954, 'A3': 2356877, 'A4': 1192324, 'A5': 864906, 'A6': 1862990, 'A7': 833566, 'A8': 1336569, 'A9': 5188707, 'A10': 869596, 'A11': 1885650, 'A12': 2640814, 'B1': 1461334, 'B2': 1430574, 'B3': 1461151, 'B4': 1387233, 'B5': 787276, 'B6': 2516497, 'B7': 3051455, 'B8': 2534009, 'B9': 2196957, 'B10': 1672517, 'B11': 2238618, 'B12': 2239936}

x_labels, y = zip(*statsA.items())

avg = sum([count for count in y])/len(statsA)
avg_y = [avg for _ in range(len(statsA))]

plt.figure(figsize=(12, 6))
plt.plot(
    y, 'g',
    avg_y, 'b',
)

plt.xticks(range(len(y)), x_labels)
plt.ylim(ymin=0)
plt.xlabel('Cells')
plt.ylabel('Reads')
plt.title('Demultiplexing Cell Sequences Pool A')
plt.show()

#%%
statsB = {'C1': 1512722, 'C2': 1097656, 'C3': 1199330, 'C4': 1345924, 'C5': 721761, 'C6': 1427646, 'C7': 855286, 'C8': 1639615, 'C9': 2516984, 'C10': 803838, 'C11': 1201528, 'C12': 1826806, 'D1': 827358, 'D2': 862868, 'D3': 759212, 'D4': 1961655, 'D5': 453363, 'D6': 1516236, 'D7': 2341031, 'D8': 1589464, 'D9': 1831497, 'D10': 1376601, 'D11': 2287288, 'D12': 4153978}
x_labels, y = zip(*statsB.items())

avg = sum([count for count in y])/len(statsB)
avg_y = [avg for _ in range(len(statsB))]

plt.figure(figsize=(12, 6))
plt.plot(
    y, 'g',
    avg_y, 'b',
)

plt.xticks(range(len(y)), x_labels)
plt.ylim(ymin=0)
plt.xlabel('Cells')
plt.ylabel('Reads')
plt.title('Demultiplexing Cell Sequences Pool B')
plt.show()


#%%
matches = {0: 784599, 1: 117785, 2: 25794, 3: 25239, 4: 19309, 5: 6086, 6: 703}
labels, sizes = zip(*matches.items())
colors = [
    '#1a9850',
    '#91cf60',
    '#d9ef8b',
    '#ffffbf',
    '#fee08b',
    '#fc8d59',
    '#d73027',
]
plt.figure(figsize=(12, 6))
plt.pie(
    sizes,
    labels=labels,
    colors=colors,
    autopct='%1.1f%%',
    startangle=180
)
plt.title('Cell A1')
plt.axis('equal')
plt.show()
