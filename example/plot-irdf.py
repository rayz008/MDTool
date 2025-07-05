import pandas as pd
import matplotlib.pyplot as plt
import sys

if (len(sys.argv) != 3):
    print(f"Usage: python plot-irdf.py <irdf.dat> <xlabel>")
    exit()

filename = sys.argv[1]
xlabel = sys.argv[2]

with open(filename, 'r') as f:
    lines = f.readlines()

# parsing irdf.dat file
irdf_data = {}
current_irdf = None
data_lines = []

for line in lines:
    line = line.strip()
    if line.startswith("iRDF:"):
        if current_irdf is not None and data_lines:
            df = pd.read_csv(pd.io.common.StringIO('\n'.join(data_lines)), 
                           sep='\t', header=None)
            irdf_data[current_irdf] = df
        
        current_irdf = int(line.split(":")[1].strip())
        data_lines = []
    elif line and not line.startswith("distance:") and not line.startswith("100"):
        data_lines.append(line)

if current_irdf is not None and data_lines:
    df = pd.read_csv(pd.io.common.StringIO('\n'.join(data_lines)), 
                   sep='\t', header=None)
    irdf_data[current_irdf] = df

# plot all iRDFs in one 
plt.figure(figsize=(8, 6))
for irdf_idx, df in irdf_data.items():
    plt.plot(df.iloc[:,0], df.iloc[:,1], label=f'iRDF {irdf_idx}', linewidth=2)

plt.xlabel(f'r({xlabel})', fontsize=16)
plt.ylabel('g(r)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend()

plt.savefig(f'irdf-{xlabel}.png')

