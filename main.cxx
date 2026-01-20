// main.cpp
#include <iostream>
#include "shape_utils.h"
#include <iomanip> // For std::fixed and std::setprecision

int main() {
    double radius = 5.0;
    double length = 4.0;
    double width = 6.0;
    
    std::cout << "hi Qiwei" << std::endl;
    // Set output to fixed point with 2 decimal places
    std::cout << std::fixed << std::setprecision(2);

    std::cout << "Circle calculations for radius " << radius << ":" << std::endl;
    std::cout << "  Area: " << calculate_circle_area(radius) << std::endl;
    std::cout << "  Perimeter: " << calculate_circle_perimeter(radius) << std::endl;

    std::cout << "\nRectangle calculations for length " << length << " and width " << width << ":" << std::endl;
    std::cout << "  Area: " << calculate_rectangle_area(length, width) << std::endl;
    std::cout << "  Perimeter: " << calculate_rectangle_perimeter(length, width) << std::endl;

    return 0;
}

