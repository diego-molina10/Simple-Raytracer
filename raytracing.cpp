#include "raytracing.h"

const int Cw = 700, Ch = 700;
const double Vw = 2.0, Vh = 2.0, d = 1.0;
const Vector O = {0, 0, 0};
const Vector BackgroundColor = {0, 0, 0};

Sphere scene[5] = {
    Sphere({0, -1, 3}, 1, {255, 0, 0}, 500, 0.2, 0.5, 1.5), // Semi-transparent
    Sphere({2, 0, 4}, 1, {0, 0, 255}, 500, 0.3, 1.0, 1.0), 
    Sphere({-2, 0, 4}, 1, {0, 255, 0}, 10, 0.4, 0.3, 2.42),
    Sphere({0, -5001, 0}, 5000, {255, 255, 0}, 1000, 0.5),
    Sphere({3, 4, 8}, 2, {255, 0, 255}, 1000, 1.0, 0.0, 1.0)
};

int amountSphere = sizeof(scene) / sizeof(scene[0]);

Light sceneLights[4] = {
    Light("ambient", 0.2),
    Light("point", 0.6, {2, 1, 0}),
    Light("directional", 0.2, {1, 4, 4}, true),
    Light("point", 0.4, {17, 6, 8})
};

int amountLight = sizeof(sceneLights) / sizeof(sceneLights[0]);

Triangle sceneTriangles[2] = {
    Triangle({0,0,2}, {1,2,2}, {-1,2,2}, {255,0,255}, 750, 0.6, 0.4, 1.0), // Semi-transparent
    Triangle({0, -1, 3}, {2, 0, 4}, {-2, 0, 4}, {255,255,255}, 1000, 0.2)
};

int amountTriangle = sizeof(sceneTriangles) / sizeof(sceneTriangles[0]);

Vector::Vector() : v(0.0, 0.0, 0.0) {}

Vector::Vector(double vx, double vy, double vz) : v(vx, vy, vz) {}

Vector operator+(const Vector& v1, const Vector& other) {
        return Vector(
            get<0>(v1.v) + get<0>(other.v),
            get<1>(v1.v) + get<1>(other.v),
            get<2>(v1.v) + get<2>(other.v)
        );
    }

Vector operator-(const Vector& v1, const Vector& other) {
        return Vector(
            get<0>(v1.v) - get<0>(other.v),
            get<1>(v1.v) - get<1>(other.v),
            get<2>(v1.v) - get<2>(other.v)
        );
    }

Vector operator*(const Vector& v1, double scalar) {
        return Vector(
            get<0>(v1.v) * scalar,
            get<1>(v1.v) * scalar,
            get<2>(v1.v) * scalar
        );
    }

double operator*(const Vector& v1, const Vector& other) {
        return get<0>(v1.v) * get<0>(other.v) +
               get<1>(v1.v) * get<1>(other.v) +
               get<2>(v1.v) * get<2>(other.v);
    }

Vector Vector::clamp(double min_val, double max_val) const {
        return Vector(
            min(max(get<0>(v), min_val), max_val),
            min(max(get<1>(v), min_val), max_val),
            min(max(get<2>(v), min_val), max_val)
        );
    }

Vector Vector::normalize() const {
        double mag = magnitude();
        return Vector(
            get<0>(v) / mag,
            get<1>(v) / mag,
            get<2>(v) / mag
        );
    }

double Vector::magnitude() const {
        return sqrt(
            get<0>(v) * get<0>(v) +
            get<1>(v) * get<1>(v) +
            get<2>(v) * get<2>(v)
        );
    }

ostream& operator<<(ostream& os, const Vector& vector) {
        os << "(" << get<0>(vector.v) << ", "
           << get<1>(vector.v) << ", "
           << get<2>(vector.v) << ")";
        return os;
    }

Sphere::Sphere(Vector vcenter, double vr, Vector vcolor, double vspecular, double vreflective, double vtransparency, double vrefractiveIndex)
    : center(vcenter), r(vr), color(vcolor), specular(vspecular), reflective(vreflective), transparency(vtransparency), refractiveIndex(vrefractiveIndex) {}

Light::Light(const string& type, double intens) 
    : lightType(type), intensity(intens), position(0, 0, 0), direction(0, 0, 0)
{}

Light::Light(const string& type, double intens, const Vector& pos)
    : lightType(type), intensity(intens), position(pos), direction(0, 0, 0)
{}

Light::Light(const string& type, double intens, const Vector& dir, bool isDirectional)
    : lightType(type), intensity(intens), position(0, 0, 0), direction(dir)
{}

Triangle::Triangle(Vector va, Vector vb, Vector vc, Vector vcolor, double vspecular, double vreflective, double vtransparency, double vrefractiveIndex)
    : a(va), b(vb), c(vc), color(vcolor), specular(vspecular), reflective(vreflective), transparency(vtransparency), refractiveIndex(vrefractiveIndex) {}

Vector Triangle::getNormal() const {
        Vector ab = b - a;
        Vector ac = c - a;
        
        Vector normal = Vector(
            get<1>(ab.v) * get<2>(ac.v) - get<2>(ab.v) * get<1>(ac.v),
            get<2>(ab.v) * get<0>(ac.v) - get<0>(ab.v) * get<2>(ac.v),
            get<0>(ab.v) * get<1>(ac.v) - get<1>(ab.v) * get<0>(ac.v)
            );
        
        return normal.normalize();
    }
    
Vector CanvasToViewport(int x, int y) {
    return {x * Vw / Cw, -y * Vh / Ch, d};
}

double dot(Vector x, Vector y) {
    return get<0>(x.v) * get<0>(y.v) + get<1>(x.v) * get<1>(y.v) + get<2>(x.v) * get<2>(y.v);
}

Vector cross(const Vector& v1, const Vector& v2) {
        return Vector(
            get<1>(v1.v) * get<2>(v2.v) - get<2>(v1.v) * get<1>(v2.v),
            get<2>(v1.v) * get<0>(v2.v) - get<0>(v1.v) * get<2>(v2.v),
            get<0>(v1.v) * get<1>(v2.v) - get<1>(v1.v) * get<0>(v2.v)
        );
    }

Vector ReflectRay(Vector R, Vector N) {
    return N * 2 * dot(N, R) - R;
}

Vector RefractRay(Vector I, Vector N, double n1, double n2) { 
    double n = n1 / n2; // Relation between refractive indexes 
    double cosI = -dot(N, I); // Is negative because we need it for the incident ray towards the surface
    double sinT2 = n * n * (1 - cosI * cosI); // Snell's law and sinI2 = (1 - cosI * cosI)
    
    if (sinT2 > 1.0) { // Total internal reflection
        return ReflectRay(I, N);
    }
    
    double cosT = sqrt(1.0 - sinT2); 
    return I * n + N * (n * cosI - cosT); 
    // I * n adjusts the incident ray by the proportion of the refraction indexes
    // N * (n * cosI - cosT) adjusts the direction of the ray to align with the new refraction angle 
}

pair<double, double> IntersectRaySphere(Vector O, Vector D, const Sphere& sphere) {
    double r = sphere.r;
    auto CO = Vector(
        get<0>(O.v) - get<0>(sphere.center.v),
        get<1>(O.v) - get<1>(sphere.center.v),
        get<2>(O.v) - get<2>(sphere.center.v)
    );

    double a = dot(D, D);
    double b = 2 * dot(CO, D);
    double c = dot(CO, CO) - r * r;

    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0)
        return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};

    double t1 = (-b + sqrt(discriminant)) / (2 * a);
    double t2 = (-b - sqrt(discriminant)) / (2 * a);

    return {t1, t2};
}

pair<double, double> IntersectRayTriangle(const Vector& O, const Vector& D, const Triangle& triangle) {
    const double EPSILON = 1e-8;
    
    Vector edge1 = triangle.b - triangle.a;
    Vector edge2 = triangle.c - triangle.a;
    
    Vector h = cross(D, edge2); // cross product
    double a = edge1 * h; // dot product
    
    if (fabs(a) < EPSILON) {
        return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()}; 
    }
    
    double f = 1.0 / a;
    Vector s = O - triangle.a;
    double u = f * (s * h); 
    
    if (u < 0.0 || u > 1.0) {
        return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()}; 
    }
    
    Vector q = cross(s, edge1); 
    double v = f * (D * q); 
    
    if (v < 0.0 || u + v > 1.0) {
        return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()}; 
    }
    
    double t = f * (edge2 * q);
    
    if (t > EPSILON) {
        return {t, numeric_limits<double>::infinity()}; 
    } else {
        return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()}; 
    }
}

pair<Sphere*, double> ClosestIntersection(const Vector& O, const Vector& D, double t_min, double t_max) {
    double closest_T = numeric_limits<double>::infinity();
    const Sphere* closestSphere = nullptr;

    for (int i = 0; i < amountSphere; ++i) {
        auto [t1, t2] = IntersectRaySphere(O, D, scene[i]);
        if (t_min <= t1 && t1 <= t_max && t1 < closest_T) {
            closest_T = t1;
            closestSphere = &scene[i];
        }
        if (t_min <= t2 && t2 <= t_max && t2 < closest_T) {
            closest_T = t2;
            closestSphere = &scene[i];
        }
    }

    return {const_cast<Sphere*>(closestSphere), closest_T};
}

pair<Triangle*, double> ClosestTriangleIntersection(const Vector& O, const Vector& D, double t_min, double t_max) {
    double closest_T = numeric_limits<double>::infinity();
    const Triangle* closestTriangle = nullptr;

    for (int i = 0; i < amountTriangle; ++i) {
        auto [t, _] = IntersectRayTriangle(O, D, sceneTriangles[i]);
        if (t_min <= t && t <= t_max && t < closest_T) {
            closest_T = t;
            closestTriangle = &sceneTriangles[i];
        }
    }

    return {const_cast<Triangle*>(closestTriangle), closest_T};
}

double ComputeLighting(const Vector& P, const Vector& N, const Vector& V, double s) {
    double intensity = 0.0;

    for (int i = 0; i < amountLight; i++) {
        Light light = sceneLights[i];

        Vector L;
        double t_max;

        if (light.lightType == "ambient") {
            intensity += light.intensity;
        } else {
            if (light.lightType == "point") {
                L = light.position - P;
                t_max = 1.0;
            } else {
                L = light.direction;
                t_max = numeric_limits<double>::infinity();
            }

            auto [shadow_sphere, shadow_t] = ClosestIntersection(P, L, 0.001, t_max);
            if (shadow_sphere != nullptr) {
                continue;
            }

            double n_dot_l = dot(N, L);
            if (n_dot_l > 0) {
                double length_N = N.magnitude();
                double length_L = L.magnitude();
                intensity += light.intensity * (n_dot_l / (length_N * length_L));
            }

            if (s != -1) {
                Vector R = N * 2 * dot(N, L) - L;
                double r_dot_v = dot(R, V);
                if (r_dot_v > 0) {
                    double length_R = R.magnitude();
                    double length_V = V.magnitude();
                    intensity += light.intensity * pow(r_dot_v / (length_R * length_V), s);
                }
            }
        }
    }

    return intensity;
}

Vector TraceRay(const Vector& O, const Vector& D, double t_min, double t_max, int depth) {
    if (depth <= 0) {
        return BackgroundColor;
    }

    auto [closest_sphere, closest_t_sphere] = ClosestIntersection(O, D, t_min, t_max);
    auto [closest_triangle, closest_t_triangle] = ClosestTriangleIntersection(O, D, t_min, t_max);

    if (closest_sphere == nullptr && closest_triangle == nullptr) {
        return BackgroundColor;
    }

    Vector P, N;
    Vector local_color;
    double specular, reflective, transparency, refractiveIndex;

    if (closest_sphere != nullptr && closest_t_sphere < closest_t_triangle) {
        P = O + D * closest_t_sphere;
        N = (P - closest_sphere->center).normalize();
        local_color = closest_sphere->color;
        specular = closest_sphere->specular;
        reflective = closest_sphere->reflective;
        transparency = closest_sphere->transparency;
        refractiveIndex = closest_sphere->refractiveIndex;
    } else if (closest_triangle != nullptr) {
        P = O + D * closest_t_triangle;
        N = closest_triangle->getNormal();
        local_color = closest_triangle->color;
        specular = closest_triangle->specular;
        reflective = closest_triangle->reflective;
        transparency = closest_triangle->transparency;
        refractiveIndex = closest_triangle->refractiveIndex;
    }

    Vector V = D * (-1);
    double lighting = ComputeLighting(P, N, V, specular);
    local_color = local_color * lighting;

    Vector reflected_color = BackgroundColor;
    if (reflective > 0) {
        Vector R = ReflectRay(V, N);
        reflected_color = TraceRay(P, R, 0.001, numeric_limits<double>::infinity(), depth - 1);
    }

    Vector refracted_color = BackgroundColor;
    if (transparency > 0) {
        double currentRefractiveIndex = 1.0; // Assumes air as the initial medium
        Vector R = RefractRay(V, N, currentRefractiveIndex, refractiveIndex); // Resulting Refracted Ray 
        refracted_color = TraceRay(P, R, 0.001, numeric_limits<double>::infinity(), depth - 1);
    }

    return local_color * (1 - reflective - transparency) + reflected_color * reflective + refracted_color * transparency;
}
