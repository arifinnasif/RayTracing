#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <cmath>

#include <iostream>

using namespace std;

#ifndef UTILS
#define UTILS

#define EPSILON 0.0000001

double deg2rad(double deg) {
    return deg/180*M_PI;
}

double rad2deg(double rad) {
    return rad*180/M_PI;
}

struct point {
    GLdouble x, y, z;
};

typedef struct {
    GLdouble r, g, b;
} color;

typedef struct point vec;
typedef struct point point;


vec cross_product(vec a, vec b) {
    vec ret;
    ret.x = a.y*b.z - a.z*b.y;
    ret.y = a.z*b.x - a.x*b.z;
    ret.z = a.x*b.y - a.y*b.x;
    return ret;
}

double dot_product(vec a, vec b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

double mod(vec arg) {
    return sqrt(arg.x * arg.x + arg.y * arg.y + arg.z * arg.z);
}

vec scale_to_r(vec arg, GLdouble r){
    vec ret;
    double m = mod(arg);
    ret.x = arg.x * r / m;
    ret.y = arg.y * r / m;
    ret.z = arg.z * r / m;
    return ret;
}

void print_vec(vec arg) {
    // print only 6 decimal places
    printf("(%.6lf, %.6lf, %.6lf)\n", arg.x, arg.y, arg.z);
}

// function that solves simultaneous equations with 3 variables
// returns a result
bool solve_simultaneous_equations(double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double a3, double b3, double c3, double d3, double *x, double *y, double *z) {
    double det = a1*b2*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1 - a2*b1*c3 - a1*b3*c2;
    if (det == 0) return false;
    *x = (d1*b2*c3 + d2*b3*c1 + d3*b1*c2 - d3*b2*c1 - d2*b1*c3 - d1*b3*c2)/det;
    *y = (a1*d2*c3 + a2*d3*c1 + a3*d1*c2 - a3*d2*c1 - a2*d1*c3 - a1*d3*c2)/det;
    *z = (a1*b2*d3 + a2*b3*d1 + a3*b1*d2 - a3*b2*d1 - a2*b1*d3 - a1*b3*d2)/det;
    return true;
}

// operator overloading to add two vectors
vec operator+(vec a, vec b) {
    vec ret;
    ret.x = a.x + b.x;
    ret.y = a.y + b.y;
    ret.z = a.z + b.z;
    return ret;
}

// operator overloading to subtract two vectors
vec operator-(vec a, vec b) {
    vec ret;
    ret.x = a.x - b.x;
    ret.y = a.y - b.y;
    ret.z = a.z - b.z;
    return ret;
}

// operator overloading to negate a vector
vec operator-(vec a) {
    vec ret;
    ret.x = -a.x;
    ret.y = -a.y;
    ret.z = -a.z;
    return ret;
}


vec get_reflection_direction(vec incident, vec normal) {
    return scale_to_r(incident - scale_to_r(normal, 2*dot_product(incident, normal)), 1);
}

#endif