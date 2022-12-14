from repairanalysis import medrasrepair
from repairanalysis import analyzeAberrations
import io
from contextlib import redirect_stdout
import numpy as np
import matplotlib.pyplot as plt
import os

sddpath = "SDD_Data"

apoptotic_base_rate = 0.011909
nhej_misrep_rate = 0.01463
single_hit_misrep_lethality = 0.189
single_hit_misrep_rate = nhej_misrep_rate*single_hit_misrep_lethality

#print('\n\nMisrepair spectrum analysis')
# medrasrepair.repairSimulation(sddpath,'Spectrum')

# Options for misrepair spectrum output
medrasrepair.doPlot = False
medrasrepair.allFragments = False
medrasrepair.simulationLimit = 24*14  # Hours, time at which to simulate misrepair
medrasrepair.listAcentrics = False

out = []
nbOfReps = 300

# Get output from misrepair spectrum analysis
f = io.StringIO()
with redirect_stdout(f):
    for i in range(nbOfReps):
        analyzeAberrations.checkHeader.__defaults__ = None, [False] # Prevent mutable arguments in checkHeader from printing headers only in the first iteration
        medrasrepair.repairSimulation(sddpath, 'Spectrum')
out = f.getvalue()

survivals = []
lines = out.split("\n")
senescent = []
#senescent_1 = []

#print(out)

repIndex = [i for i, e in enumerate(lines) if e.startswith("From folder:")]
repIndex.append(len(lines)-1) # Last line
reps = []

for index in range(len(repIndex)-1):
    reps.append(lines[repIndex[index]:repIndex[index+1]+1])

for rep in reps:
    fileIndex = [i for i, e in enumerate(rep) if e.startswith("Data for:")]

    files = []

    # Filter out unnecessary fields. Leave only: 
    # total Breaks, Residual, Misrepairs, 
    # Large Misrepairs, Inter-Chromosome Misrepairs, Acentric Linear, 
    # Large Acentric Fragment, Multi-Centric, Centric Ring, 
    # Acentric Ring, Multi-Centric Ring, Large Ring Fragment, 
    # Single Hit Repairs

    headers_complete = rep[fileIndex[0]+1].split("\t")

    headers = [e for e in headers_complete if e not in (
        'Index', 'Single-Junction Chromosomes', 'Multi-Junction Chromosomes', 'Normal Chromosomes', 'Initial DNA Fragmentation', 'Potential DNA loss')]
    headers_indices = [i for i in range(
        len(headers_complete)) if headers_complete[i] in headers]

    if len(fileIndex) == 1:
        files.append(rep[fileIndex[0]+2:-1])  # Skip header in first file
    else:
        for f in range(len(fileIndex)):
            if f == 0:    # First
                files.append(rep[fileIndex[f]+2:fileIndex[f+1]-1])  # Skip header
            elif f == len(fileIndex)-1:  # Last
                files.append(rep[fileIndex[-1]+1:-1])   # Last file
            else:   # Middle
                files.append(rep[fileIndex[f]+1:fileIndex[f+1]-1])

    for file in range(len(files)):
        for line in range(len(files[file])):
            files[file][line] = files[file][line].split("\t")

    for file in files:
        summaries = np.zeros(len(headers)+1)  # Add one for single hit repair
        survival = 0
        for line in file:
            if line[1] != "No misrepair!":
                for n in range(len(headers)):
                    summaries[n] += int(line[headers_indices[n]])
        #for n, field in enumerate(headers):
            #print("total "+field+" "+str(summaries[n]))
        #print("\n")
        single_hit_misrep_survival = np.exp(-single_hit_misrep_rate*summaries[13])
        #single_hit_arrest_survival = np.exp(-apoptotic_base_rate*summaries[1])
        single_hit_arrest_survival = np.exp(-apoptotic_base_rate*summaries[0])
        aberration_survival=np.exp(-(summaries[6]+summaries[7]+summaries[8]+summaries[10]+summaries[11]))
        #aberration_survival=np.exp(-sum(summaries[6:11])) # Test with acentric rings
        survival=single_hit_misrep_survival*single_hit_arrest_survival*aberration_survival
        survivals.append(survival)
        senescent.append(1-(single_hit_misrep_survival*single_hit_arrest_survival))
        #senescent_1.append(1-(single_hit_arrest_survival))

print(survivals)
print("\nAvg survival "+str(np.average(survivals))+"\n")
print("\nAvg senescent "+str(np.average(senescent))+"\n")
#print("\nAvg senescent single hit arrest"+str(np.average(senescent_1))+"\n")
plt.hist(survivals, 100, range=[0., 1.], label="Survival, n = "+str(nbOfReps))
plt.hist(senescent, 100, range=[0., 1.], label="Senescent, n = "+str(nbOfReps))
#plt.hist(senescent_1, 100, range=[0., 1.], label="Senescent single hit arrest, n = "+str(nbOfReps))
plt.legend()
plt.show() 
