#include "Point.hpp"
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

double Point::operator()(int i) const
{
    switch (i)
      {
         case 1:
            return(x_);
         case 2:
            return(y_);
         case 3:
            return(z_);
         default:
            cout << "Indice different de 1, 2 ou 3." << endl;
            exit(1);
      }
}

double & Point::operator()(int i)
{
    switch (i)
      {
         case 1:
            return(x_);
         case 2:
            return(y_);
         case 3:
            return(z_);
         default:
            cout << "Indice different de 1, 2 ou 3." << endl;
            exit(1);
      }
}

double Point::operator[](int i) const
{
    switch (i)
      {
         case 0:
            return(x_);
         case 1:
            return(y_);
         case 2:
            return(z_);
         default:
            cout << "Indice different de 0, 1 ou 2." << endl;
            exit(1);
      }
}

double & Point::operator[](int i)
{
    switch (i)
      {
         case 0:
            return(x_);
         case 1:
            return(y_);
         case 2:
            return(z_);
         default:
            cout << "Indice different de 0, 1 ou 2." << endl;
            exit(1);
      }
}

Point & Point::operator=(const Point & P)
{
    x_ = P.x_;
    y_ = P.y_;
    z_ = P.z_;
    return(*this);
}

Point Point::operator+=(const Point & P)
{
    x_ += P.x_;
    y_ += P.y_;
    z_ += P.z_;
    return(*this);
}

Point Point::operator-=(const Point & P)
{
    x_ -= P.x_;
    y_ -= P.y_;
    z_ -= P.z_;
    return(*this);
}

Point Point::operator*=(double a)
{
    x_ *= a;
    y_ *= a;
    z_ *= a;
    return(*this);
}

Point Point::operator/=(double a)
{
    if (a == 0.0) {
        cout << "Division par zero." << endl;
        exit(1);
    }
    x_ /= a;
    y_ /= a;
    z_ /= a;
    return(*this);
}

Point Point::operator^=(const Point & P)
{
    double x(y_*P.z_ - z_*P.y_);
    double y(z_*P.x_ - x_*P.z_);
    double z(x_*P.y_ - y_*P.x_);
    x_ = x; y_ = y; z_ = z;
    return(*this);
}

void Point::print() const
{
    cout << "(" << x_ << ", " << y_ << ", " << z_ << ")" << endl;
}

Point operator+(const Point & P, const Point & Q)
{
    Point R(P);
    R += Q;
    return(R);
}

Point operator-(const Point & P, const Point & Q)
{
    Point R(P);
    R -= Q;
    return(R);
}

Point operator*(double a, const Point & P)
{
    Point R(P);
    R *= a;
    return(R);
}

Point operator*(const Point & P, double a)
{
    Point R(P);
    R *= a;
    return(R);
}

Point operator/(const Point & P, double a)
{
    Point R(P);
    R /= a;
    return(R);
}

double operator|(const Point & P, const Point & Q)
{
    double ps;
    ps = P[0] * Q[0] + P[1] * Q[1] + P[2] * Q[2];
    return(ps);
}

Point operator^(const Point & P, const Point & Q)
{
    Point R(P);
    R ^= Q;
    return(R);
}

double norm_euc(const Point & P)
{
    return(sqrt(P | P));
}

double SignedVolume(const Point & A, const Point & B, const Point & C, const Point & D)
{
    return((1./6.) * ((B-A) ^ (C-A)) | (D-A));
}


