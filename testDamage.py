# Add damagegenerator to path as preferred
from damagegenerator import damageModel

damageModel.basicXandIon()  # Generate an illustrative dataset
damageModel.generateExposure(10.0, 4.58, 1, 1, 10) # Generate event for 10 MeV proton
