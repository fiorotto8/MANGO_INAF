
def GEMsGen():

    ElementGems=["U238","Th232","U235","K40","Co60","Cs137"]
    ZElementGems=[92,90,92,19,27,55]
    NElements=[238,232,235,40,60,137]
    RadEl="GEMs"
    
    for i in range(len(ElementGems)):
        with open(RadEl+"_"+ElementGems[i]+".mac", 'w') as f:
            f.write("/process/eLoss/StepFunction 0.01 0.1 mm"+"\n")
            f.write("/run/verbose 0"+"\n")
            f.write("/event/verbose 0"+"\n")
            f.write("/random/setSeeds 575223199 8381234975"+"\n")
            f.write("/output/OutFile "+RadEl+"_"+ElementGems[i]+"\n")
            f.write("#Elements: Cathodes GEMs Rings Vessel Lens Sensors"+"\n")
            f.write("/detector/RadElement "+RadEl+"\n")
            f.write("/isotope/AtomicNumber "+str(ZElementGems[i])+"\n")
            f.write("/isotope/MassNumber "+str(NElements[i])+"\n")
            f.write("/run/beamOn 10000000"+"\n")

def VesselGen():

    ElementGems=["U238","Th232","K40"]
    ZElementGems=[92,90,19]
    NElements=[238,232,40]
    RadEl="Vessel"
    
    for i in range(len(ElementGems)):
        with open(RadEl+'_'+ElementGems[i]+".mac", 'w') as f:
            f.write("/process/eLoss/StepFunction 0.01 0.1 mm"+"\n")
            f.write("/run/verbose 0"+"\n")
            f.write("/event/verbose 0"+"\n")
            f.write("/random/setSeeds 575223199 8381234975"+"\n")
            f.write("/output/OutFile "+RadEl+"_"+ElementGems[i]+"\n")
            f.write("#Elements: Cathodes GEMs Rings Vessel Lens Sensors"+"\n")
            f.write("/detector/RadElement "+RadEl+"\n")
            f.write("/isotope/AtomicNumber "+str(ZElementGems[i])+"\n")
            f.write("/isotope/MassNumber "+str(NElements[i])+"\n")
            f.write("/run/beamOn 10000000"+"\n")
        

def FieldCageGen():

    ElementGems=["U238","Th232","K40","Co60","Cs137"]
    ZElementGems=[92,90,19,27,55]
    NElements=[238,232,40,60,137]
    RadEl="Rings"
    
    for i in range(len(ElementGems)):
        with open(RadEl+'_'+ElementGems[i]+".mac", 'w') as f:
            f.write("/process/eLoss/StepFunction 0.01 0.1 mm"+"\n")
            f.write("/run/verbose 0"+"\n")
            f.write("/event/verbose 0"+"\n")
            f.write("/random/setSeeds 575223199 8381234975"+"\n")
            f.write("/output/OutFile "+RadEl+"_"+ElementGems[i]+"\n")
            f.write("#Elements: Cathodes GEMs Rings Vessel Lens Sensors"+"\n")
            f.write("/detector/RadElement "+RadEl+"\n")
            f.write("/isotope/AtomicNumber "+str(ZElementGems[i])+"\n")
            f.write("/isotope/MassNumber "+str(NElements[i])+"\n")
            f.write("/run/beamOn 10000000"+"\n")

def CathodesGen():

    ElementGems=["U238","U234"]
    ZElementGems=[92,92]
    NElements=[238,234]
    RadEl="Cathodes"
    
    for i in range(len(ElementGems)):
        with open(RadEl+'_'+ElementGems[i]+".mac", 'w') as f:
            f.write("/process/eLoss/StepFunction 0.01 0.1 mm"+"\n")
            f.write("/run/verbose 0"+"\n")
            f.write("/event/verbose 0"+"\n")
            f.write("/random/setSeeds 575223199 8381234975"+"\n")
            f.write("/output/OutFile "+RadEl+"_"+ElementGems[i]+"\n")
            f.write("#Elements: Cathodes GEMs Rings Vessel Lens Sensors"+"\n")
            f.write("/detector/RadElement "+RadEl+"\n")
            f.write("/isotope/AtomicNumber "+str(ZElementGems[i])+"\n")
            f.write("/isotope/MassNumber "+str(NElements[i])+"\n")
            f.write("/run/beamOn 10000000"+"\n")

def LensGen():

    ElementGems=["U238","Th232"]
    ZElementGems=[92,90]
    NElements=[238,232]
    RadEl="Lens"
    
    for i in range(len(ElementGems)):
        with open(RadEl+'_'+ElementGems[i]+".mac", 'w') as f:
            f.write("/process/eLoss/StepFunction 0.01 0.1 mm"+"\n")
            f.write("/run/verbose 0"+"\n")
            f.write("/event/verbose 0"+"\n")
            f.write("/random/setSeeds 575223199 8381234975"+"\n")
            f.write("/output/OutFile "+RadEl+"_"+ElementGems[i]+"\n")
            f.write("#Elements: Cathodes GEMs Rings Vessel Lens Sensors"+"\n")
            f.write("/detector/RadElement "+RadEl+"\n")
            f.write("/isotope/AtomicNumber "+str(ZElementGems[i])+"\n")
            f.write("/isotope/MassNumber "+str(NElements[i])+"\n")
            f.write("/run/beamOn 10000000"+"\n")

def SensorGen():

    ElementGems=["U238","Th232","K40","Co60","Cs137"]
    ZElementGems=[92,90,19,27,55]
    NElements=[238,232,40,60,137]
    RadEl="Sensors"
    
    for i in range(len(ElementGems)):
        with open(RadEl+'_'+ElementGems[i]+".mac", 'w') as f:
            f.write("/process/eLoss/StepFunction 0.01 0.1 mm"+"\n")
            f.write("/run/verbose 0"+"\n")
            f.write("/event/verbose 0"+"\n")
            f.write("/random/setSeeds 575223199 8381234975"+"\n")
            f.write("/output/OutFile "+RadEl+"_"+ElementGems[i]+"\n")
            f.write("#Elements: Cathodes GEMs Rings Vessel Lens Sensors"+"\n")
            f.write("/detector/RadElement "+RadEl+"\n")
            f.write("/isotope/AtomicNumber "+str(ZElementGems[i])+"\n")
            f.write("/isotope/MassNumber "+str(NElements[i])+"\n")
            f.write("/run/beamOn 10000000"+"\n")


if __name__ == "__main__":

    GEMsGen()
    VesselGen()
    FieldCageGen()
    CathodesGen()
    LensGen()
    SensorGen()
