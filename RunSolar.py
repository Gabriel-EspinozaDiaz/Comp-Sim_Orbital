from Solar import Solar
#from Solar import SolarSatelliteScan
from Solar import SolarNullPlanets

def main():

    
    s = Solar(satellite=True)

    s.run()
    

main()