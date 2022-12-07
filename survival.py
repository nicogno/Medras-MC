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

headers = lines[fileIndex[0]+1].split("\t")
headers = headers[1:-2]

for f in range(len(fileIndex)-1):
    if f==0:
        files.append(lines[fileIndex[f]+2:fileIndex[f+1]-1]) # Skip header
    else:
        files.append(lines[fileIndex[f]+1:fileIndex[f+1]-1])
files.append(lines[fileIndex[-1]+1:-1])   # Last file

for file in range(len(files)):
    for line in range(len(files[file])):
        files[file][line] = files[file][line].split("\t")


for file in files:
    summaries = np.zeros(len(headers))
    for line in file:
        if line[1]!="No misrepair!":
            for n, field in enumerate(line[1:-2]):
                summaries[n] += int(field)
    for n, field in enumerate(headers):
        print("total "+field+" "+str(summaries[n]))
    print("\n")
