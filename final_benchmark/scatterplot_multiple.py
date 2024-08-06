import matplotlib.pyplot as plt
import sys
import numpy as np
from sklearn.metrics import r2_score 

files = sys.argv[1:]

plt.figure(figsize=(12, 8))

colors = ['red', 'blue', 'green', 'orange']

for i, file in enumerate(files):
    with open(file) as f:
        proteins = {}
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            id = line[0]
            size = line[1]
            contacts = line[2]
            time = line[3]
            
            proteins[id] = []
            proteins[id].append(size)
            proteins[id].append(time)
    
    ids = list(proteins.keys())
    sizes = [float(entry[0]) for entry in proteins.values()]
    times = [float(entry[1]) for entry in proteins.values()]
    #contacts = [float(entry[1]) for entry in proteins.values()]

    sizes = np.array(sizes)
    times = np.array(times)

    model_linear = np.poly1d(np.polyfit(sizes, times, 1)) 
    model_quadratic = np.poly1d(np.polyfit(sizes, times, 2))
    
    polyline = np.linspace(min(sizes), max(sizes), 100) 

    linear_polyline = model_linear(polyline)
    quadratic_polyline = model_quadratic(polyline) 

    print(f"Model {i+1}:")
    print(file, colors[i])
    r2_linear = r2_score(times, model_linear(sizes))
    print(f"R2 linear: {r2_linear}")

    r2_quadratic = r2_score(times, model_quadratic(sizes))
    print(f"R2 quadratic: {r2_quadratic}")
    print()

    plt.scatter(sizes, times, color=colors[i], edgecolor='k', alpha=0.7)
    plt.plot(polyline, linear_polyline, label=f'Model{i+1} ({colors[i]})\n    Linear Fit\n    $R^2 = {r2_linear:.2f}$') 
    plt.plot(polyline, quadratic_polyline, label=f'    Quadratic Fit\n    $R^2 = {r2_quadratic:.2f}$\n')

plt.xlabel('Protein Size', fontsize=14)
plt.ylabel('Time (s)', fontsize=14)
plt.title('Time vs Protein Size', fontsize=16)

# Annotate each point with the protein ID
for contact, time, id in zip(sizes, times, ids):
    plt.annotate(id, (contact, time), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=9)

plt.tight_layout()
plt.legend()
plt.show()

# models
# print(model_linear)
# print(model_quadratic) 