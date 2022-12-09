from repairanalysis import medrasrepair
import io
from contextlib import redirect_stdout
import numpy as np

sddpath = "SDD_Data"

print('\n\nMisrepair spectrum analysis')
#medrasrepair.repairSimulation(sddpath,'Spectrum')

# Options for misrepair spectrum output
medrasrepair.doPlot = False
medrasrepair.allFragments = False
medrasrepair.simulationLimit = 24 # Hours, time at which to simulate misrepair
medrasrepair.listAcentrics = False

# Get output from misrepair spectrum analysis
f = io.StringIO()
with redirect_stdout(f):
    medrasrepair.repairSimulation(sddpath,'Spectrum')
out = f.getvalue()

print(out)

lines = out.split("\n")

fileIndex = [i for i, e in enumerate(lines) if e.startswith("Data for:")]

files = []

headers_complete = lines[fileIndex[0]+1].split("\t")
headers = [e for e in headers_complete if e not in ('Index', 'Single-Junction Chromosomes', 'Multi-Junction Chromosomes', 'Normal Chromosomes', 'Initial DNA Fragmentation', 'Potential DNA loss')]
headers_indices = [i for i in range(len(headers_complete)) if headers_complete[i] in headers]

if len(fileIndex) == 1:
    files.append(lines[fileIndex[0]+2:-1]) # Skip header in first file
else:
    for f in range(len(fileIndex)):
        if f==0:    # First 
            files.append(lines[fileIndex[f]+2:fileIndex[f+1]-1]) # Skip header
        elif f==len(fileIndex)-1: # Last
            files.append(lines[fileIndex[-1]+1:-1])   # Last file
        else:   # Middle
            files.append(lines[fileIndex[f]+1:fileIndex[f+1]-1])

for file in range(len(files)):
    for line in range(len(files[file])):
        files[file][line] = files[file][line].split("\t")

for file in files:
    summaries = np.zeros(len(headers)+1) # Add one for single hit repair
    for line in file:
        if line[1]!="No misrepair!":
            for n in range(len(headers)):
                summaries[n]+=int(line[headers_indices[n]])
            summaries[-1] += int(line[1])-int(line[2])-int(line[3]) # Total breaks - Residual - Misrepaired
    for n, field in enumerate(headers):
        print("total "+field+" "+str(summaries[n]))
    print("Single Hit Repairs "+str(summaries[-1]))
    print("\n")
