
#ifndef CONSTANTS_H
#define CONSTANTS_H

/* Namespace containing constants used in several files. */
namespace cnst {
    const int num_dims = 3;                       // Use dimensions x, y and z
    const double G = 6.67408E-11*2.97449433E-19;  // [m^3*kg^(-1)*s^(-2)] --> [AU^3*kg^(-1)*yr^(-2)]
    const double mass_sun = 2.0E+30;
    const double pi = acos(-1.0);
    const double four_pi_sq = 4*acos(-1.0)*acos(-1.0);
    const double c = 173.0*365.0;               // Speed of light [AU/yr]
    const double mercury_year = 90.0/365.0;     // Time to get around sun for mercury [yr]
}

#endif // CONSTANTS_H
