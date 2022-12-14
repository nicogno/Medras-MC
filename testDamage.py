# Add damagegenerator to path as preferred
from damagegenerator import damageModel

#damageModel.basicXandIon()  # Generate an illustrative dataset
#damageModel.generateExposure(10.0, 4.58, 1, 1, 10) # Generate event for 10 MeV proton
damageModel.generateExposure(0.001, 0, 15, 0, 1) # Generate event for 15 Gy photons, 1 KeV
