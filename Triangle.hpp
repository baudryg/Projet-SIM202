#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <vector>
#include "Point.hpp"

using namespace std;

class Triangle
{
public:
    Point A_;
    Point B_;
    Point C_;
    Point bary_;
    Point projBary_; // uniquement pour le cas de l'observateur a l'infini

public:
    Triangle();
    Triangle(const Point &, const Point &, const Point &);
    Triangle & operator=(const Triangle &);
    void print() const;

    double maxz();
    bool contain(const Point &);
    bool isIntersected(const Point &); // uniquement pour le cas de l'observateur a l'infini
    bool isIntersected(const Point &, const Point &);
    bool isIntersectedObs(const Point &, const Point &);
    Point ptIntersection(const Point &); // uniquement pour le cas de l'observateur a l'infini
    Point ptIntersection(const Point &, const Point &);
    Point ptIntersectionObs(const Point &, const Point &);
    Point vectNormal();
    double aire();
    std::vector<double> coordBarycentric(const Point &);
    vector<double> angleOmbragePlat(const Point &);

};

#endif // M
