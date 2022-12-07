# Add repairanalysis to path as preferred
from repairanalysis import medrasrepair

sddpath = "SDD_Data"

#print('Fidelity analysis')
#medrasrepair.repairSimulation(sddpath,'Fidelity')

#print('\n\nSeparation analysis')
#medrasrepair.repairSimulation(sddpath,'Separation')

print('\n\nMisrepair spectrum analysis')
medrasrepair.repairSimulation(sddpath,'Spectrum')
