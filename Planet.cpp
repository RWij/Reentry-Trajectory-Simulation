#include "Planet.h"

/* ============================================= Earth ============================================= */
Earth::Earth() : Planet(9.8062, 6378000.0, 7250.0) 
{
    this->setAtmsDensity_kgpm3(1.226);
    this->setAtmsPressure_Pa(101325.0);
    this->setAtmsTemperature_K(288.19);
}

/* ============================================= Mars ============================================= */
Mars::Mars() : Planet(3.71, 3380000.0, 11100.0) 
{
    this->setAtmsDensity_kgpm3(0.057);
    this->setAtmsPressure_Pa(600.0);
    this->setAtmsTemperature_K(242.15);
}