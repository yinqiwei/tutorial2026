// shape_utils.cpp
#include "shape_utils.h"
#include <cmath> // For M_PI

// Function definitions
double calculate_circle_area(double radius) {
    return M_PI * radius * radius;
}

double calculate_circle_perimeter(double radius) {
    return 2 * M_PI * radius;
}

double calculate_rectangle_area(double length, double width) {
    return length * width;
}

double calculate_rectangle_perimeter(double length, double width) {
    return 2 * (length + width);
}

