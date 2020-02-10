#ifndef MAILLAGE_H
#define MAILLAGE_H

#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include "Point.hpp"
#include "Triangle.hpp"

using namespace std;

class Maillage {
protected:
    std::vector<Point> vPoint_;
    std::vector<Triangle> vTriangle_;

    std::vector<std::vector<int> > vPointAroundTriangle_;
    std::vector<std::vector<int> > vTriangleAroundPoint_;

public:
    Maillage() {vPoint_.resize(0); vTriangle_.resize(0);}
    Maillage(vector<Point> vPoint, vector<Triangle> vTriangle) : vPoint_(vPoint), vTriangle_(vTriangle) {}
    int gNbPt() {return vPoint_.size();}
    int gNbTri() {return vTriangle_.size();}

    Point operator[](int) const;
    Triangle operator()(int) const;

    void loadGMSH(const std::string &);

    bool crossTriangle(const Point &); // uniquement pour le cas de l'observateur a l'infini
    bool crossTriangle(const Point &, const Point &);
    bool crossTriangleObs(const Point &, const Point &);
    bool crossOtherTriangle(int, const Point &, const Point &);
    Triangle closestIntersectedTriangle(const Point &); // uniquement pour le cas de l'observateur a l'infini
    Triangle closestIntersectedTriangle(const Point &, const Point &);
    int closestIntersectedTriangleId(const Point & P); // uniquement pour le cas de l'observateur a l'infini
    int closestIntersectedTriangleId(const Point &, const Point &);
    int closestIntersectedTriangleIdObs(const Point &, const Point &);

    std::vector<int> getPointAroundTriangle(int);
    std::vector<int> getTriangleAroundPoint(int);

    Point vectNormPoint(int);
    vector<double> angleOmbrage(int, const Point &, const Point &);
};

template <typename Out>
static inline void split(const std::string & s, char delim, Out result)
{
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

static inline std::vector<std::string> split(const std::string & s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

#endif // MAILLAGE_H
