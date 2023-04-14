from Solar import Solar
from Solar import SolarSatelliteScan
from Solar import SolarNullPlanets

def main():

    
    s = SolarNullPlanets(satellite=False)

    s.run()
    

main()