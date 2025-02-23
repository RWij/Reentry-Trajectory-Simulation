#include <iostream>
#include <iomanip>
#include <cmath>

using errorCode = int;

static const long double PI      = acos(-1);
static const long double deg2rad = 1.0; //PI / 180.0;
static const long double rad2deg = 1.0; //180.0 / PI;
static const long double Atms2Pa = 101325.;

constexpr int PRECISION        = 5;
constexpr int WIDTH            = 11;

enum Planet_enum { Earth, Mars };

struct Variables 
{
    double flight_path_angle_deg;      // flight path angle, deg
    double velocity_mps;               // entry/mid-flight velocity, km/s
    double ballistic_coefficient_nd;   // ballistic coefficient, non-dimensional
    double nose_radius_m;              // effective nose radius, m
    double lift_over_drag_nd;          // lift over drag ratio, non-dimensionl
    double bank_angle_deg;             // bank angle, deg
    double altitude_m;                 // entry/mid-flight interface altitude, km
    
public:
    /** initialize input variables **/
    Variables() : flight_path_angle_deg(0), velocity_mps(0), 
        ballistic_coefficient_nd(0), nose_radius_m(0), lift_over_drag_nd(0),
        bank_angle_deg(0), altitude_m(0) {}

    Variables(double flight_path_angle_deg, double velocity_mps,
        double ballistic_coefficient_nd, double nose_radius_m,
        double lift_over_drag_nd, double bank_angle_deg, double altitude_m) :
        flight_path_angle_deg(0), velocity_mps(0), ballistic_coefficient_nd(0), 
        nose_radius_m(0), lift_over_drag_nd(0), bank_angle_deg(0), altitude_m(0) {}
};

struct Planet 
{
    double gravity_mps2;
    double atms_density_kgpm3;
    double atms_pressure_Pa;
    double radius_m;
    double geoheight_m;
public:
    Planet() : gravity_mps2(0), atms_density_kgpm3(0), 
    atms_pressure_Pa(0), radius_m(0), geoheight_m(0) {}

    Planet(double gravity_mps2, double atms_density_kgpm3, 
        double atms_pressure_Pa, double radius_m, double geoheight_m) : 
        gravity_mps2(gravity_mps2), atms_density_kgpm3(atms_density_kgpm3),
        atms_pressure_Pa(atms_pressure_Pa), radius_m(radius_m), geoheight_m(geoheight_m) {}
};

struct Atmosphere 
{
    long double temperature_k;
    long double pressure_atms;
    long double rho_kgpm3;
public:
    /** initialize input variables **/
    Atmosphere() : temperature_k(0), pressure_atms(0), rho_kgpm3(0) {}

    Atmosphere(double temperature_k, double pressure_atms, double rho_kgpm3) : 
    temperature_k(temperature_k), pressure_atms(pressure_atms), rho_kgpm3(rho_kgpm3) {}
};

struct RungeKuttaVariables 
{
    double velocity_mps;
    double flight_path_deg;
    double altitude_m;
    double deltaS;
    double delta_velocity_mps;
    double Vnmax;
    double hnmax;
    double nmax;

public:
    RungeKuttaVariables() : velocity_mps(0), flight_path_deg(0), 
    altitude_m(0), deltaS(0), delta_velocity_mps(0), 
    Vnmax(0), hnmax(0), nmax(0) {}

    RungeKuttaVariables(double velocity_mps, double flight_path_deg,
        double altitude_m, double deltaS, double delta_velocity_mps, 
        double Vnmax, double hnmax, double nmax) : velocity_mps(velocity_mps), flight_path_deg(flight_path_deg), 
        altitude_m(altitude_m), deltaS(deltaS), delta_velocity_mps(delta_velocity_mps), 
        Vnmax(Vnmax), hnmax(hnmax), nmax(nmax) {}
};

Planet initialize_atmospheric_conditions(Planet_enum planet_enum) 
{

    if (planet_enum == Earth)
     {
        Planet planet = 
        {
            9.8062,         // gravity, mps2
            1.226,          // atmospheric density, kgpm3
            101325.0,       // atmospheric pressure, Pa
            6378000.0,      // planet radius, m
            7250.0          // geopotential height, m
        };
        return planet;
       
    }
    else if (planet_enum == Mars) 
    {
        Planet planet = 
        {
            3.71,           // gravity, mps2
            0.057,          // atmospheric density, kgpm3
            600.0,          // atmospheric pressure, Pa
            3380000.0,      // planet radius, m
            11100.0         // geopotential height, m
        };
        return planet;
    }
}

Atmosphere calculate_atmospheric_properties(double altitude_m, Planet_enum planet_enum) 
{
    Atmosphere atmosphere;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    double pressure_atms = 0.0;
    long double rho_kgpm3 = 0.0;

    if (planet_enum == Earth) 
    {
        if (altitude_m <= 11000) 
        {
            temperature_k = 288.19 - 0.00649 * altitude_m;
            pressure_pa = 101290 * pow(temperature_k / 288.08, 5.256);
            if (temperature_k < 0) 
            {
                temperature_k = 2.73;
            }
            rho_kgpm3 = ((pressure_pa / 1000) / (.2869 * (temperature_k)));
        }
        else if (altitude_m > 11000 && altitude_m < 25000) 
        {
            temperature_k = 216.69;
            pressure_pa = 0.02265 * exp(1.73 - ((0.000157) * altitude_m));
            if (temperature_k < 0) 
            {
                temperature_k = 2.73;
            }
            rho_kgpm3 = ((pressure_pa / 1000) / (.2869 * (temperature_k)));
        }
        else if (altitude_m >= 25000) 
        {
            temperature_k = 141.94 + 0.00229 * altitude_m;
            pressure_pa = 2488 * pow(temperature_k / 216.6, -11.388);

            if (temperature_k < 0) 
            {
                temperature_k = 2.73;
            }
            rho_kgpm3 = (pressure_pa / 1000) / (0.2869 * temperature_k);
        }

        pressure_atms = pressure_pa / 101325.;
    }
    else if (planet_enum == Mars) 
    {
        if (altitude_m < 7000) 
        {
            temperature_k = 242.15 - 0.00222 * altitude_m;
            pressure_pa = 699 * exp(-0.00009 * altitude_m);
            if (temperature_k < 0) 
            {
                temperature_k = 2.73;
            }
            rho_kgpm3 = abs((pressure_pa / 1000) / (.1921 * (temperature_k)));
        }
        else if (altitude_m >= 7000) 
        {
            temperature_k = 249.75 - 0.00222 * altitude_m;
            pressure_pa = 669 * exp(-0.00009 * altitude_m);
            rho_kgpm3 = abs((pressure_pa / 1000) / (.1921 * (temperature_k)));
        }
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

    atmosphere.temperature_k = temperature_k;
    atmosphere.pressure_atms = pressure_atms;
    atmosphere.rho_kgpm3 = rho_kgpm3;

    return atmosphere;
}

RungeKuttaVariables calculate_next_runge_kutta_step(RungeKuttaVariables rkv, Variables var, Atmosphere& atms, Planet& planet, double g, double dt) 
{
    RungeKuttaVariables new_rkv;
    /** Runge-Kutta Numerical Scheme **/
    double y1k1, y2k1, y3k1, y1k2, y2k2, y3k2, y1k3,
        y2k3, y3k3, y1k4, y2k4, y3k4, y1k5, y2k5, y3k5, y1k6, y2k6, y3k6;
    double k11, k21, k31, k12, k22, k32, k13,
        k23, k33, k14, k24, k34, k15, k25, k35, k16, k26, k36;

    double rhot = atms.rho_kgpm3;
    double Re = planet.radius_m;
    double B = var.ballistic_coefficient_nd;
    double LoD = var.lift_over_drag_nd;
    double bank = var.bank_angle_deg * deg2rad;

    double Vt, gammat, alt_t, dS, dvdt; 

    y1k1 = rkv.velocity_mps;
    y2k1 = rkv.flight_path_deg * deg2rad;
    y3k1 = rkv.altitude_m;
    k11 = (-rhot * pow(y1k1,2)) / (2.0 * B) + g * sin(y2k1);
    k21 = (-y1k1 * cos(y2k1)) / (Re + y3k1)-((rhot * y1k1) / (2.0 * B)) * (LoD)*cos(bank)+(g / y1k1) * cos(y2k1);
    k31 = -y1k1 * sin(y2k1);

    y1k2 = y1k1 + 0.25 * k11 * dt;
    y2k2 = y2k1 + 0.25 * k21 * dt;
    y3k2 = y3k1 + 0.25 * k31 * dt;
    k12 = (-rhot * pow(y1k2,2)) / (2.0 * B) + g * sin(y2k2);
    k22 = (-y1k2 * cos(y2k2)) / (Re + y3k2) - ((rhot * y1k2) / (2.0 * B)) * (LoD)*cos(bank) + (g / y1k2) * cos(y2k2);
    k32 = -y1k2 * sin(y2k2);

    y1k3 = y1k1 + 0.125 * k11 * dt + 0.125 * k12 * dt;
    y2k3 = y2k1 + 0.125 * k21 * dt + 0.125 * k22 * dt;
    y3k3 = y3k1 + 0.125 * k31 * dt + 0.125 * k32 * dt;
    k13 = (-rhot * pow(y1k3, 2)) / (2.0 * B) + g * sin(y2k3);
    k23 = (-y1k3 * cos(y2k3)) / (Re + y3k3) - ((rhot * y1k3) / (2.0 * B)) * (LoD)*cos(bank) + (g / y1k3) * cos(y2k3);
    k33 = -y1k3 * sin(y2k3);

    y1k4 = y1k1 - 0.5 * k12 * dt + k13 * dt;
    y2k4 = y2k1 - 0.5 * k22 * dt + k23 * dt;
    y3k4 = y3k1 - 0.5 * k32 * dt + k33 * dt;
    k14 = (-rhot * pow(y1k4, 2)) / (2.0 * B) + g * sin(y2k4);
    k24 = (-y1k4 * cos(y2k4)) / (Re + y3k4) - ((rhot * y1k4) / (2.0 * B)) * (LoD)*cos(bank) + (g / y1k4) * cos(y2k4);
    k34 = -y1k4 * sin(y2k4);

    y1k5 = y1k1 + 0.1875 * k11 * dt + 0.5625 * k14 * dt;
    y2k5 = y2k1 + 0.1875 * k21 * dt + 0.5625 * k24 * dt;
    y3k5 = y3k1 + 0.1875 * k31 * dt + 0.5625 * k34 * dt;
    k15 = (-rhot * pow(y1k5, 2)) / (2.0 * B) + g * sin(y2k5);
    k25 = (-y1k5 * cos(y2k5)) / (Re + y3k5)-((rhot * y1k5) / (2.0 * B)) * (LoD)*cos(bank)+(g / y1k5) * cos(y2k5);
    k35 = -y1k5 * sin(y2k5);

    y1k6 = y1k1 - 0.428571429 * k11 * dt + 0.285714286 * k12 * dt + 1.71428571429 * k13 * dt - 1.71428571429 * k14 * dt 
        + 1.14285714286 * k15 * dt;
    y2k6 = y2k1 - 0.428571429 * k21 * dt + 0.285714286 * k22 * dt + 1.71428571429 * k23 * dt - 1.71428571429 * k24 * dt 
        + 1.14285714286 * k25 * dt;
    y3k6 = y3k1 - 0.428571429 * k31 * dt + 0.285714286 * k32 * dt + 1.71428571429 * k33 * dt - 1.71428571429 * k34 * dt 
        + 1.14285714286 * k35 * dt;
    k16 = (-rhot * pow(y1k6, 2)) / (2.0 * B) + g * sin(y2k6);
    k26 = (-y1k6 * cos(y2k6)) / (Re + y3k6) - ((rhot * y1k6) / (2.0 * B)) * (LoD)*cos(bank) + (g / y1k6) * cos(y2k6);
    k36 = -y1k6 * sin(y2k6);

    Vt      = y1k1 + 0.011111111 * (7.0 * k11 + 32. * k13 + 12. * k14 + 32. * k15 + 7. * k16) * dt;
    gammat  = y2k1 + 0.011111111 * (7.0 * k21 + 32. * k23 + 12. * k24 + 32. * k25 + 7. * k26) * dt;
    alt_t   = y3k1 + 0.011111111 * (7.0 * k31 + 32. * k33 + 12. * k34 + 32. * k35 + 7. * k36) * dt;
    dS      = -1. / tan(y2k1) * (0.011111111 * (7.0 * k31 + 32. * k33 + 12. * k34 + 32. * k35 + 7. * k36) * dt);

    dvdt = -((rhot * pow(Vt,2)) / (2.0 * B)) + g * sin(gammat);
    new_rkv.nmax = std::min(new_rkv.nmax, dvdt / planet.gravity_mps2);

    if (new_rkv.nmax == dvdt / planet.gravity_mps2) 
    {
        new_rkv.Vnmax = y1k1;
        new_rkv.hnmax = y3k1;
    }

    new_rkv.deltaS              = dS;
    new_rkv.velocity_mps        = Vt;
    new_rkv.altitude_m          = alt_t;
    new_rkv.flight_path_deg     = gammat * rad2deg;
    new_rkv.delta_velocity_mps  = dvdt;
    
    return new_rkv;
}

double calculate_radiation_heat(double v, double rho, double rnose, Planet_enum planet_selection) 
{
    double a, b, c, fv, dfdv, qrad;

    double vedat[19] = 
    {
        9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750, 11000,
        11500, 12000, 12500, 13000, 13500, 14000, 14500, 15000, 15500, 16000
    };
    double fevdat[19] = 
    {
        1.5, 4.3, 9.7, 19.5, 35., 55., 81., 115., 151., 238., 359.,
        495., 660., 850., 1065., 1313., 1550., 1780., 2040.
    };
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
        
    if (planet_selection == Earth) 
    {
        a = 1.072e+06 * pow(v, -1.88) * pow(rho, -0.325);
        a = std::min(a, 1.0);
        if (rnose > 1.0 && rnose < 2.0) 
        {
            a = std::min(a, 0.6);
        }
        if (rnose > 2.0 && rnose < 3.0) 
        {
            a = std::min(a, 0.5);
        }
        b = 1.22;
        c = 4.736e+04;
                    
        if (v < vedat[0]) 
        {
            fv = fevdat[0];
        }
        else 
        {
            fv = fevdat[19];    // if v > any value listed in vedat
            for (size_t i = 1; i < 18; i++) 
            {
                if (v < vedat[i]) 
                {
                    dfdv = (fevdat[i] - fevdat[i - 1]) / (vedat[i] - vedat[i - 1]);
                    fv = fevdat[i - 1] + (v - vedat[i - 1]) * dfdv;
                }
            }
        }

        if (v < 9000) 
        {
            c = 0.;
        }
    }
    else if (planet_selection == Mars) 
    {
        a = 0.526;
        b = 1.19;
        c = 2.35e+04;

        if (v < vmdat[0]) 
        {
            fv = fmvdat[0];
        }
        else 
        {
            fv = fmvdat[16];        // if v > any value listed in vmdat
            for (size_t i = 1; i < 16; i++) 
            {
                if (v < vmdat[i]) 
                {
                    dfdv = (fmvdat[i] - fmvdat[i - 1]) / (vmdat[i] - vmdat[i - 1]);
                    fv = fmvdat[i - 1] + (v - vmdat[i - 1]) * dfdv;
                }
            }
        }
        if (v < 5500) 
        {
            c = 0.;
        }
    }         
    
    qrad = c * pow(rnose, a) * pow(rho, b) * fv;
    return qrad;
}

errorCode run_reentry(Planet_enum planet_selection, Variables& var, double dtime_sec, double dtime_out_sec) 
{
    Atmosphere atmosphere;
    RungeKuttaVariables rkv;
    Planet planet = initialize_atmospheric_conditions(planet_selection);

    errorCode returnCode = 0; // 0 = returned successfully, -1 = error encountered
    
    bool isRunning  = true;

    /** Initialize Loop Variables **/
    double time     = 0.0;
    double Cp       = 1.27;
    double prntint  = 1.0;
    double tprnt    = 1.0;
    double nmax     = 100.0;
    double n        = 0.0;
    double qtotal   = 0.0;
    double qradbc   = 0.0;
    double s        = 0.0;
    double Hrec     = (pow(var.velocity_mps, 2.) / 2.) / 2325.854324;
    double ruCh     = 1.0e-6;
    double lam      = 4;
    double g        = 0.0;
    double dynprs   = 0.0;
    double dynprs0  = 0.0;
    double surfP    = 0.0; 
    double acceleration                 = 0.0;

    double prev_altitude_m              = var.altitude_m;
    double curr_altitude_m              = var.altitude_m;
    double curr_velocity_mps            = var.velocity_mps;
    double prev_velocity_mps            = var.velocity_mps;
    double prev_flight_path_angle_deg   = var.flight_path_angle_deg;
    double curr_flight_path_angle_deg   = var.flight_path_angle_deg;
    double curr_rho                     = planet.atms_density_kgpm3;
    double prev_qconv                   = 0.0;  
    double curr_qconv                   = 0.0;
    double prev_qrad                    = 0.0;
    double curr_qrad                    = 0.0;
    double prev_qtotal                  = 0.0;
    double curr_qtotal                  = 0.0;
    double prev_deltaS                  = 0.0;
    double curr_deltaS                  = 0.0;
    double prev_delta_velocity_mps      = 0.0;
    double curr_delta_velocity_mps      = 0.0;

    std::cout << std::setw(WIDTH) << "Time(s)"          ;
    std::cout << std::setw(WIDTH) << "Alt(m)"           ;
    std::cout << std::setw(WIDTH) << "Vel(m/s)"         ;
    std::cout << std::setw(WIDTH) << "Patm(N/m2)"       ;
    std::cout << std::setw(WIDTH) << "Tatm(K)"          ;
    std::cout << std::setw(WIDTH) << "rho(kg/m3)"       ;
    std::cout << std::setw(WIDTH) << "Gam(deg)"         ;
    std::cout << std::setw(WIDTH) << "Decel(g)"         ;
    std::cout << std::setw(WIDTH) << "Qc(W/cm2)"        ;
    std::cout << std::setw(WIDTH) << "Qr(W/Cm2)"        ;
    std::cout << std::setw(WIDTH) << "Qt(W/cm2)"        ;
    std::cout << std::setw(WIDTH) << "Qbc(W/cm2)"       ;
    std::cout << std::setw(WIDTH) << "S(m)"             ;
    std::cout << std::setw(WIDTH) << "Pdyn(N/m2)"       ;
    std::cout << std::setw(WIDTH) << "Hrec(J/K)"        ;
    std::cout << std::setw(WIDTH) << "ruCh"             ;
    std::cout << std::setw(WIDTH) << "Psurf(N/m2)"      ;
    std::cout << std::endl;
    std::cout << "-----------------------------------------------";
    std::cout << "-----------------------------------------------";
    std::cout << "-----------------------------------------------";
    std::cout << "-----------------------------------------------" << std::endl;

    // initialize atmosphere and for printing initial values
    atmosphere = calculate_atmospheric_properties(curr_altitude_m, planet_selection);

    surfP = atmosphere.pressure_atms + dynprs0 * Cp;
    
    std::cout << std::setprecision(PRECISION);
    std::cout << std::setw(WIDTH) << std::round(time);
    std::cout << std::setw(WIDTH) << curr_altitude_m;
    std::cout << std::setw(WIDTH) << curr_velocity_mps;
    std::cout << std::setw(WIDTH) << atmosphere.pressure_atms * Atms2Pa;
    std::cout << std::setw(WIDTH) << atmosphere.temperature_k;
    std::cout << std::setw(WIDTH) << atmosphere.rho_kgpm3;
    std::cout << std::setw(WIDTH) << curr_flight_path_angle_deg;
    std::cout << std::setw(WIDTH) << n;
    std::cout << std::setw(WIDTH) << curr_qconv;
    std::cout << std::setw(WIDTH) << curr_qrad;
    std::cout << std::setw(WIDTH) << curr_qtotal;
    std::cout << std::setw(WIDTH) << qradbc;
    std::cout << std::setw(WIDTH) << s;
    std::cout << std::setw(WIDTH) << dynprs;
    std::cout << std::setw(WIDTH) << Hrec;
    std::cout << std::setw(WIDTH) << ruCh;
    std::cout << std::setw(WIDTH) << surfP;
    std::cout << std::endl;

    while (curr_altitude_m > 0.0 && isRunning) 
    {
        /** setup starting conditions for the next interval **/
        rkv.altitude_m          = curr_altitude_m;
        rkv.velocity_mps        = curr_velocity_mps;
        rkv.flight_path_deg     = curr_flight_path_angle_deg;
        rkv.delta_velocity_mps  = curr_delta_velocity_mps;

        // print the starting conditions
        if (time >= tprnt) {
            std::cout << std::setprecision(PRECISION);
            std::cout << std::setw(WIDTH) << std::round(time);
            std::cout << std::setw(WIDTH) << curr_altitude_m;
            std::cout << std::setw(WIDTH) << curr_velocity_mps;
            std::cout << std::setw(WIDTH) << atmosphere.pressure_atms * Atms2Pa;
            std::cout << std::setw(WIDTH) << atmosphere.temperature_k;
            std::cout << std::setw(WIDTH) << atmosphere.rho_kgpm3;
            std::cout << std::setw(WIDTH) << curr_flight_path_angle_deg;
            std::cout << std::setw(WIDTH) << n;
            std::cout << std::setw(WIDTH) << curr_qconv;
            std::cout << std::setw(WIDTH) << curr_qrad;
            std::cout << std::setw(WIDTH) << curr_qtotal;
            std::cout << std::setw(WIDTH) << qradbc;
            std::cout << std::setw(WIDTH) << s;
            std::cout << std::setw(WIDTH) << dynprs;
            std::cout << std::setw(WIDTH) << Hrec;
            std::cout << std::setw(WIDTH) << ruCh;
            std::cout << std::setw(WIDTH) << surfP;
            std::cout << std::endl;
            tprnt = tprnt + prntint;
        }

        g = planet.gravity_mps2 / pow(1 + (curr_altitude_m / planet.radius_m), 2); 

        rkv = calculate_next_runge_kutta_step(rkv, var, atmosphere, planet, g, dtime_sec);

        s += rkv.deltaS;

        dynprs = .5 * atmosphere.rho_kgpm3 * pow(rkv.velocity_mps, 2);

        curr_rho = atmosphere.rho_kgpm3;

        curr_qrad = calculate_radiation_heat(curr_velocity_mps, atmosphere.rho_kgpm3, var.nose_radius_m, planet_selection);
        qradbc = curr_qrad * 0.881;

        if (planet_selection == Earth) {
            curr_qconv = (1.74153e-4 * pow(curr_rho / var.nose_radius_m, 0.5) * pow(curr_velocity_mps, 3)) / 10000.0;
        }
        else if (planet_selection == Mars) {
            curr_qconv = (1.9027e-4 * pow(curr_rho / var.nose_radius_m, 0.5) * pow(curr_velocity_mps, 3)) / 10000.0;
        }

        curr_qtotal = prev_qtotal + ((prev_qconv + curr_qconv) / 2.0) * dtime_sec + ((prev_qrad + curr_qrad) / 2.0) * dtime_sec;
        
        Hrec = pow(curr_velocity_mps, 2.) / 2.;

        ruCh = (curr_qconv * 10000) / Hrec;

        Hrec = Hrec / 2325.854324;

        ruCh = ruCh / 4.882;

        surfP = atmosphere.pressure_atms + dynprs0 * Cp;

        rkv.delta_velocity_mps = -(curr_rho * pow(curr_velocity_mps, 2)) / (2.0 * var.ballistic_coefficient_nd) + g * sin(curr_flight_path_angle_deg * deg2rad);
        
        n = rkv.delta_velocity_mps / g;

        // update loop variables
        prev_altitude_m             = curr_altitude_m;
        curr_altitude_m             = rkv.altitude_m;
        prev_velocity_mps           = curr_velocity_mps;
        curr_velocity_mps           = rkv.velocity_mps;
        prev_flight_path_angle_deg  = curr_flight_path_angle_deg;
        curr_flight_path_angle_deg  = rkv.flight_path_deg;
        prev_deltaS                 = curr_deltaS;
        curr_deltaS                 = rkv.deltaS;
        prev_delta_velocity_mps     = curr_delta_velocity_mps;
        curr_delta_velocity_mps     = rkv.delta_velocity_mps;
        prev_qconv                  = curr_qconv;
        prev_qrad                   = curr_qrad;
        prev_qtotal                 = curr_qtotal;

        acceleration = (curr_velocity_mps - prev_velocity_mps) / dtime_sec / g;
        if (curr_altitude_m > prev_altitude_m + 100) {
            std::cout << "Flight path angle too shallow for the given entry velocity." << std::endl; 
            std::cout << "Vehicle will skip out of the atmosphere, make flight path angle larger" << std::endl;
            isRunning = false;
            returnCode = -1;
            break;
        }

        time += dtime_sec;

        // update atmosphere data
        atmosphere = calculate_atmospheric_properties(curr_altitude_m, planet_selection);
    }

    return returnCode;
}

int main()
{
    Variables vars;
    vars.velocity_mps               = 5000.0;
    vars.altitude_m                 = 60000.0;
    vars.flight_path_angle_deg      = 50.0;
    vars.ballistic_coefficient_nd   = 3000.0;
    vars.nose_radius_m              = 3.0;
    vars.lift_over_drag_nd          = 2.3;
    vars.bank_angle_deg             = 0.0;    

    /** select planet **/
    Planet_enum planet_selection = Earth; // Either 'Earth' or 'Mars'

    /** time intervals **/
    double dtime_sec           = 0.05; // Time steps in calculations
    double dtime_out_sec       = 1.00; // Time steps for print out

    errorCode simErrorCode;
    simErrorCode = run_reentry(planet_selection, vars, dtime_sec, dtime_out_sec);
    std::cout << "Returned " << simErrorCode << "\n" << std::endl;
    return simErrorCode;
}
