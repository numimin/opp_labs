#!/usr/bin/env python

import matplotlib.pyplot as plt
from hurry.filesize import size
import sys

#read data:
count = int(sys.stdin.readline())

data = [sys.stdin.readline() for _ in range(count * 2)]

arguments = [data[i * 2] for i in range(count)]
fig, ax = plt.subplots()

ax.plot(arguments, [float(data[i * 2 + 1]) for i in range(count)])

plt.ylabel('Time (seconds)')
plt.xlabel('Cores [rows x cols]')

plt.setp(ax.xaxis.get_minorticklabels(), rotation=90)
plt.xticks(rotation=90)

ax.tick_params(axis='both', which='major', labelsize=5)
ax.tick_params(axis='both', which='minor', labelsize=4)

fig.savefig("./lab2_plot.png", dpi=150)
