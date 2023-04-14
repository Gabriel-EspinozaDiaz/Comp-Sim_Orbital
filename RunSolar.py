from Solar import Solar
from Solar import SolarNullPlanets

def main():

    s = Solar(satellite=True)
    sNull = SolarNullPlanets(satellite=False)

    s.run()
    
main()