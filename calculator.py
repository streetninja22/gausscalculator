import math

#the maximum test value for p, the conic parameter or semi-latus rectum
G_CONIC_PARAMETER_TEST_VALUE = 340649140000000000

#solves the Gauss problem via the p-iteration technique
def solveGauss(gravitationalParameter, radiusInitial, radiusTarget, velocityInitial, velocityTarget, deltaTrueAnomaly, timeOfFlight):
    
    #absolute values for radii
    rIAbsolute = abs(radiusInitial)
    rTAbsolute = abs(radiusTarget)

    #generically named constants. I'm not sure exactly what these values represent
    kConstant = rIAbsolute * rTAbsolute * (1 - math.cos(deltaTrueAnomaly))
    lConstant = rIAbsolute + rTAbsolute
    mConstant = rIAbsolute * rTAbsolute * (1 + math.cos(deltaTrueAnomaly))
    
    
    #minimum and maximum values of the semi-latus rectum depending on value of the change in true anomaly
    parameterMin = 0
    parameterMax = 0
    
    if deltaTrueAnomaly < math.pi:
        parameterMin = kConstant / (lConstant + math.sqrt(2 * mConstant))
        parameterMax = G_CONIC_PARAMETER_TEST_VALUE
    else:
        parameterMin = 0
        parameterMax = kConstant / (lConstant - math.sqrt(2 * mConstant))
    
    while True:
        parameter = (parameterMin + parameterMax) / 2
        print("Testing for p value of " + str(parameter) + " (minimum of " + str(parameterMin) + ", maximum of " + str(parameterMax) + ")")
        
        semiMajorAxis = (mConstant * kConstant * parameter) / ((2 * mConstant - lConstant**2) * parameter**2 + 2 * kConstant * lConstant * parameter - kConstant**2)
        
        fFunction = 1 - (radiusTarget / parameter) * (1 - math.cos(deltaTrueAnomaly))
        gFunction = (radiusInitial * radiusTarget * math.sin(deltaTrueAnomaly)) / math.sqrt(gravitationalParameter * parameter)
        
        #test value for time-of-flight
        timeOfFlightTest = 0
        
        if semiMajorAxis >= 0:
            print("Positive a")
            eccentricAnomaly = math.acos(1 - (radiusInitial / semiMajorAxis) * (1 - fFunction))
            
            timeOfFlightTest = gFunction + math.sqrt((semiMajorAxis**3) / gravitationalParameter) * (eccentricAnomaly - math.sin(eccentricAnomaly))
        else:
            print("a: " + str(semiMajorAxis))
            eccentricAnomaly = math.acosh(1 - (radiusInitial / semiMajorAxis) * (1 - fFunction))
            print("delta F: " + str(eccentricAnomaly))
            
            timeOfFlightTest = gFunction + math.sqrt(((-semiMajorAxis)**3) / gravitationalParameter) * (math.sinh(eccentricAnomaly) - eccentricAnomaly)
        
        if math.isclose(timeOfFlight, timeOfFlightTest, rel_tol=0.05):
            return parameter
        elif timeOfFlightTest > timeOfFlight:
            parameterMin = parameter
        else:
            parameterMax = parameter
            
        print("t: " + str(timeOfFlightTest))
        


print("Standard gravitational parameter of parent body:")
gravitationalParameter = float(input())

print("Intitial radius:")
radiusInitial = float(input())
print("Radius of target body:")
radiusTarget = float(input())

print("Initial velocity:")
velocityInitial = float
(input())
print("Velocity of target body:")
velocityTarget = float(input())

print("Enter the change in true anomaly (radians)")
deltaTrueAnomaly = float(input())

print("Input transfer time (seconds)")
timeOfFlight = input()

#test values
gravitationalParameter = 1.1723328e18
radiusInitial =  13599840256
radiusTarget =  19669121365
velocityInitial = 9285
velocityTarget = 7147
deltaTrueAnomaly = 2.82743
timeOfFlight = 5.53e+6


print("P-Value: " + str(solveGauss(gravitationalParameter, radiusInitial, radiusTarget, velocityInitial, velocityTarget, deltaTrueAnomaly, timeOfFlight)))