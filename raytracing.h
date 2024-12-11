#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <tuple>

using namespace std;

class Vector {
public:
    tuple<double, double, double> v;

    Vector();
    Vector(double vx, double vy, double vz);
    friend Vector operator+(const Vector& v1, const Vector& other);
    friend Vector operator-(const Vector& v1, const Vector& other);
    friend Vector operator*(const Vector& v1, double scalar);
    friend double operator*(const Vector& v1, const Vector& other);
    Vector clamp(double min_val, double max_val) const;
    Vector normalize() const;
    double magnitude() const;
    friend ostream& operator<<(ostream& os, const Vector& vector);
};

class Sphere {
public:
    Vector center;
    double r;
    Vector color;
    double specular;
    double reflective;
    double transparency; // New
    double refractiveIndex; // New

    Sphere(Vector vcenter, double vr, Vector vcolor, double vspecular, double vreflective, double vtransparency = 0.0, double vrefractiveIndex = 1.0);
};

class Light {
public:
    string lightType;
    double intensity;
    Vector position;
    Vector direction;

    Light(const string& type, double intens);
    Light(const string& type, double intens, const Vector& pos);
    Light(const string& type, double intens, const Vector& dir, bool isDirectional);
};

class Triangle {
public:
    Vector a, b, c;
    Vector color;
    double specular;
    double reflective;
    double transparency; // New
    double refractiveIndex; // New

    Triangle(Vector va, Vector vb, Vector vc, Vector vcolor, double vspecular, double vreflective, double vtransparency = 0.0, double vrefractiveIndex = 1.0);
    Vector getNormal() const;
};

extern const int Cw, Ch;
extern const double Vw, Vh, d;
extern const Vector O, BackgroundColor;

extern Sphere scene[5];
extern int amountSphere;

extern Light sceneLights[4];
extern int amountLight;

extern Triangle sceneTriangles[2];
extern int amountTriangle;

Vector CanvasToViewport(int x, int y);
double dot(Vector x, Vector y);
Vector cross(const Vector& v1, const Vector& v2);
Vector ReflectRay(Vector R, Vector N);
Vector RefractRay(Vector I, Vector N, double n1, double n2);
pair<double, double> IntersectRaySphere(Vector O, Vector D, const Sphere& sphere);
pair<double, double> IntersectRayTriangle(const Vector& O, const Vector& D, const Triangle& triangle);
pair<Sphere*, double> ClosestIntersection(const Vector& O, const Vector& D, double t_min, double t_max);
pair<Triangle*, double> ClosestTriangleIntersection(const Vector& O, const Vector& D, double t_min, double t_max);
double ComputeLighting(const Vector& P, const Vector& N, const Vector& V, double s);
Vector TraceRay(const Vector& O, const Vector& D, double t_min, double t_max, int depth);

#endif // RAYTRACING_H
