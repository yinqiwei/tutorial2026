#include <vector>

class Mass {
    public:
        std::vector<double> mass = {
            0.0000000, 1.007825, 4.002602, 7.01600, 9.012182, 11.00931, 12.00000,
            14.00307, 15.99491, 18.99840, 19.99244 
        };

    double printmass (int a);
}; 

//Mass::mass[i]