import pandas as pd
import matplotlib.pyplot as plt
import sys

if (len(sys.argv) != 3):
    print(f"Usage: python plot-rdf.py <rdf.dat> <xlabel>")
    exit()

df = pd.read_csv(sys.argv[1], sep='\s+', header=None, skiprows=2)
xlabel = sys.argv[2]

plt.plot(df.iloc[:,0], df.iloc[:,1])
plt.xlabel(F'r({xlabel})', fontsize=16)
plt.ylabel('g(r)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=14)

plt.savefig(f'rdf-{xlabel}.png')
