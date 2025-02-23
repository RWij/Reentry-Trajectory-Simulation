#ifndef __PLANET__
#define __PLANET__

#include <cmath>

class Planet 
{
protected:
    // intrinsic data
    double gravity_mps2;
    double radius_m;
    double geoheight_m;

    // atmospheric data
    double atms_temperature_k;    
    double atms_density_kgpm3;
    double atms_pressure_Pa;

    // heat transfer data
    double qrad_wpcm2;
    double qconv_wpcm2;

    /* Constructors */
    Planet(double gravity_mps2, double radius_m,  double geoheight_m) :
        gravity_mps2(gravity_mps2), radius_m(radius_m), geoheight_m(geoheight_m) {}

public:
    /* Getters and Setters*/
    double getGravity_mps2() const { return gravity_mps2; }
    double getRadius_m()     const { return radius_m;}
    double getGeoHeight_m()     const { return geoheight_m;}

    double getAtmsDensity_kgpm3()                       { return atms_density_kgpm3;}
    void setAtmsDensity_kgpm3(double rho_kgpm3)         { atms_density_kgpm3 = rho_kgpm3;}

    double getAtmsPressure_Pa()                         { return atms_pressure_Pa;}
    void setAtmsPressure_Pa(double pressure_pa)         { atms_pressure_Pa = pressure_pa;}

    double getAtmsPressure_Atms()                       { return atms_pressure_Pa / Atms2Pa;}
    void setAtmsPressure_Atms(double pressure_atms)     { atms_pressure_Pa = pressure_atms * Atms2Pa;}

    double getAtmsTemperature_K()                       { return atms_temperature_k;}
    void setAtmsTemperature_K(double temperature_k)     { atms_temperature_k = temperature_k;}

    double getQrad_wpcm2()                              { return qrad_wpcm2;}
    void setQrad_wpcm2(double new_qrad_wpcm2)           { qrad_wpcm2 = new_qrad_wpcm2;}

    double getQconv_wpcm2()                             { return qconv_wpcm2;}
    void setQconv_wpcm2(double new_qconv_wpcm2)         { qconv_wpcm2 = new_qconv_wpcm2;}

    /* Override */
    virtual void calculate_atmospheric_properties(double altitude_m) = 0;
    virtual void calculate_radiation_heat(double velocity_mps, double rho_kgpm3, double rnose_m) = 0;
    virtual void calculate_convective_heat_transfer(double velocity_mps, double rnose_m) = 0;
};

/* ============================================= Earth ============================================= */
class Earth : public Planet
{
private: 
    const double vedat[19] = 
    {
        9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750, 11000,
        11500, 12000, 12500, 13000, 13500, 14000, 14500, 15000, 15500, 16000
    };
    const double fevdat[19] = 
    {
        1.5, 4.3, 9.7, 19.5, 35., 55., 81., 115., 151., 238., 359.,
        495., 660., 850., 1065., 1313., 1550., 1780., 2040.
    };
public:
    Earth();

    void calculate_atmospheric_properties(double altitude_m) override
    {
        double temperature_k = 0.0;
        double pressure_pa = 0.0;
        double pressure_atms = 0.0;
        double rho_kgpm3 = 0.0;

        if (altitude_m <= 11000) 
        {
            temperature_k = 288.19 - 0.00649 * altitude_m;
            pressure_pa = 101290 * pow(temperature_k / 288.08, 5.256);
            if (temperature_k < 0) 
            {
                temperature_k = 2.73;
            }
            rho_kgpm3 = fabs((pressure_pa / 1000) / (.2869 * (temperature_k)));
        }
        else if (altitude_m > 11000 && altitude_m < 25000) 
        {
            temperature_k = 216.69;
            pressure_pa = 0.02265 * exp(1.73 - ((0.000157) * altitude_m));
            if (temperature_k < 0) 
            {
                temperature_k = 2.73;
            }
            rho_kgpm3 = fabs((pressure_pa / 1000) / (.2869 * (temperature_k)));
        }
        else if (altitude_m >= 25000) 
        {
            temperature_k = 141.94 + 0.00229 * altitude_m;
            pressure_pa = 2488 * pow(temperature_k / 216.6, -11.388);

            if (temperature_k < 0) 
            {
                temperature_k = 2.73;
            }
            rho_kgpm3 = fabs((pressure_pa / 1000) / (0.2869 * temperature_k));
        }

        if (temperature_k < 0) 
        {
            temperature_k = 2.73;
        }
    
        if (pressure_pa < 0) 
        {
            pressure_pa = 0.0;
        }

        pressure_atms = pressure_pa / 101325.;
    
        this->setAtmsTemperature_K(temperature_k);
        this->setAtmsPressure_Atms(pressure_atms);
        this->setAtmsDensity_kgpm3(rho_kgpm3);
    }

    void calculate_radiation_heat(double velocity_mps, double rho_kgpm3, double rnose_m) override
    {
        double a, b, c, fv, dfdv, qrad;

        a = 1.072e+06 * pow(velocity_mps, -1.88) * pow(rho_kgpm3, -0.325);
        a = std::min(a, 1.0);
        if (rnose_m > 1.0 && rnose_m < 2.0) 
        {
            a = std::min(a, 0.6);
        }
        if (rnose_m > 2.0 && rnose_m < 3.0) 
        {
            a = std::min(a, 0.5);
        }
        b = 1.22;
        c = 4.736e+04;
                    
        if (velocity_mps < vedat[0]) 
        {
            fv = fevdat[0];
        }
        else 
        {
            fv = fevdat[18];    // if v > any value listed in vedat
            for (size_t i = 1; i < 18; i++) 
            {
                if (velocity_mps < vedat[i]) 
                {
                    dfdv = (fevdat[i] - fevdat[i - 1]) / (vedat[i] - vedat[i - 1]);
                    fv = fevdat[i - 1] + (velocity_mps - vedat[i - 1]) * dfdv;
                }
            }
        }
        if (velocity_mps < 9000) 
        {
            c = 0.;
        }
        qrad = c * pow(rnose_m, a) * pow(rho_kgpm3, b) * fv;

        this->setQrad_wpcm2(qrad);
    }

    void calculate_convective_heat_transfer(double velocity_mps, double rnose_m) override 
    {
        double qconv;
        qconv = (1.74153e-4 * pow(this->getAtmsDensity_kgpm3() / rnose_m, 0.5) * pow(velocity_mps, 3)) / 10000.0;
        this->setQconv_wpcm2(qconv);
    }
};

/* ============================================= Mars ============================================= */
class Mars : public Planet
{
private:
    double vmdat[17] = 
    {
        6000, 6150, 6300, 6500, 6700, 6900, 7000, 7200, 7400, 7600,
        7800, 8000, 8200, 8400, 8600, 8800, 9000
    };
    double fmvdat[17] = 
    {
        0.2, 1.0, 1.95, 3.42, 5.1, 7.1, 8.1, 10.2, 12.5, 14.8, 17.1,
        19.2, 21.4, 24.1, 26.0, 28.9, 32.8
    };
public:
    Mars();

    void calculate_atmospheric_properties(double altitude_m) override
    {
        double temperature_k = 0.0;
        double pressure_pa = 0.0;
        double pressure_atms = 0.0;
        double rho_kgpm3 = 0.0;

        if (altitude_m < 7000) 
        {
            temperature_k = 242.15 - 0.00222 * altitude_m;
            if (temperature_k < 0) 
            {
                temperature_k = 2.73;
            }
            pressure_pa = 699 * exp(-0.00009 * altitude_m);
            rho_kgpm3 = fabs((pressure_pa / 1000) / (.1921 * (temperature_k)));
        }
        else if (altitude_m >= 7000) 
        {
            temperature_k = 249.75 - 0.00222 * altitude_m;
            if (temperature_k < 0) 
            {
                temperature_k = 2.73;
            }
            pressure_pa = 669 * exp(-0.00009 * altitude_m);
            rho_kgpm3 = fabs((pressure_pa / 1000) / (.1921 * (temperature_k)));
        }

        if (temperature_k < 0) 
        {
            temperature_k = 2.73;
        }
    
        if (pressure_pa < 0) 
        {
            pressure_pa = 0.0;
        }

        pressure_atms = pressure_pa / 101325.;
    
        this->setAtmsTemperature_K(temperature_k);
        this->setAtmsPressure_Atms(pressure_atms);
        this->setAtmsDensity_kgpm3(rho_kgpm3);
    }

    void calculate_radiation_heat(double velocity_mps, double rho_kgpm3, double rnose_m) override
    {
        double a, b, c, fv, dfdv, qrad;

        a = 0.526;
        b = 1.19;
        c = 2.35e+04;

        if (velocity_mps < vmdat[0]) 
        {
            fv = fmvdat[0];
        }
        else 
        {
            fv = fmvdat[16];        // if v > any value listed in vmdat
            for (size_t i = 1; i < 16; i++) 
            {
                if (velocity_mps < vmdat[i]) 
                {
                    dfdv = (fmvdat[i] - fmvdat[i - 1]) / (vmdat[i] - vmdat[i - 1]);
                    fv = fmvdat[i - 1] + (velocity_mps - vmdat[i - 1]) * dfdv;
                }
            }
        }
        if (velocity_mps < 5500) 
        {
            c = 0.;
        }
        qrad = c * pow(rnose_m, a) * pow(rho_kgpm3, b) * fv;
        this->setQrad_wpcm2(qrad);
    }

    void calculate_convective_heat_transfer(double velocity_mps, double rnose_m) override 
    {
        double qconv;
        qconv = (1.9027e-4 * pow(this->getAtmsDensity_kgpm3() / rnose_m, 0.5) * pow(velocity_mps, 3)) / 10000.0;
        this->setQconv_wpcm2(qconv);
    }
};

#endif // __PLANET__