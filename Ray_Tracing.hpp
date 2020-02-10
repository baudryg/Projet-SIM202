#ifndef RAY_TRACING_H
#define RAY_TRACING_H

#include "Image.hpp"
#include "Point.hpp"
#include "Maillage.hpp"

class RayTracing
{
protected:
    Maillage maillageRT_;
    Image imgRT_;
    Point srcLumRT_; // Unique source de lumiere pour l'instant
    Point obsRT_;
    double Ia, Is;
    double ka, kd, kr;
    int n;

public:
    RayTracing(const Maillage & maillageRT, const Image & imgRT, const Point & srcLumRT, const Point & obsRT)
    : maillageRT_(maillageRT), imgRT_(imgRT), srcLumRT_(srcLumRT), obsRT_(obsRT),
      Ia(1.0), Is(1.0), ka(0.3), kd(0.5), kr(0.0), n(6) {}

    double intensiteFromAngle(double, double);
    double intensiteFromAngle();
    double intensiteGouraud(int, const Point &, const Point &);
    void colorImagePlat();
    void colorImageGouraud();
    void afficherImage();
};




#endif // RAY_TRACING_H
