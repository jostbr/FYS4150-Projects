
#ifndef PLANET_H
#define PLANET_H

# include <iostream>

class Planet{
    public:
        Planet();
        Planet(std::string planet_name, double x, double y, double v_x, double v_y);
        ~Planet();
    private:
        std::string planet_name;
        double r[2], v[2], a[2];
};

#endif // PLANET_H
