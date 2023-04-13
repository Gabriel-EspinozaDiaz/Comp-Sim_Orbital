from Solar import Solar
from Solar import SolarSatelliteScan
from Solar import SolarNullPlanets

def main():

    
    s = Solar(satellites=False)
    sScan = SolarSatelliteScan(satellites=True)
    sNull = SolarNullPlanets(satellites=False)

    s.runWoAnimation()
    

main()