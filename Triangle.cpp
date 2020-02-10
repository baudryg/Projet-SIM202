#include "Triangle.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

Triangle::Triangle()
{
    Point A, B, C;
    A_ = A;
    B_ = B;
    C_ = C;
    bary_ = (A_ + B_ + C_)/3;
    projBary_ = bary_;
    projBary_[2] = 0.0;
}

Triangle::Triangle(const Point & A, const Point & B, const Point & C)
: A_(A), B_(B), C_(C)
{
    bary_ = (A_ + B_ + C_)/3;
    projBary_ = bary_;
    projBary_[2] = 0.0;
}

Triangle & Triangle::operator=(const Triangle & T)
{
    A_ = T.A_;
    B_ = T.B_;
    C_ = T.C_;
    bary_ = (A_ + B_ + C_)/3;
    projBary_ = bary_;
    projBary_[2] = 0.0;
    return(*this);
}

void Triangle::print() const
{
    cout << "A : ";
    A_.print();
    cout << "B : ";
    B_.print();
    cout << "C : ";
    C_.print();
}

double Triangle::maxz()
{
    return(max(max(A_[2], B_[2]), C_[2]));
}

bool Triangle::contain(const Point & P)
{
    return(true); // A REMPLIR !
}

bool Triangle::isIntersected(const Point & P1)
{
    Point P2(P1);
    P2[2] = 2*maxz();
    double vol_1, vol_2, vol_3, vol_4, vol_5;
    vol_1 = SignedVolume(P1, A_, B_, C_);
    vol_2 = SignedVolume(P2, A_, B_, C_);
    vol_3 = SignedVolume(P1, P2, A_, B_);
    vol_4 = SignedVolume(P1, P2, B_, C_);
    vol_5 = SignedVolume(P1, P2, C_, A_);
    return((vol_1 * vol_2 < 0) && (vol_3 * vol_4 > 0) && (vol_4 * vol_5 > 0));
}

bool Triangle::isIntersected(const Point & P1, const Point & P2)
{
    double vol_1, vol_2, vol_3, vol_4, vol_5;
    vol_1 = SignedVolume(P1, A_, B_, C_);
    vol_2 = SignedVolume(P2, A_, B_, C_);
    vol_3 = SignedVolume(P1, P2, A_, B_);
    vol_4 = SignedVolume(P1, P2, B_, C_);
    vol_5 = SignedVolume(P1, P2, C_, A_);
    return((vol_1 * vol_2 < 0) && (vol_3 * vol_4 > 0) && (vol_4 * vol_5 > 0));
}

bool Triangle::isIntersectedObs(const Point & P_obs, const Point & P_pixel)
{
    Point P_max((P_pixel - P_obs) * (1 + 2 * maxz()/(P_pixel[2] - P_obs[2])) + P_obs);
    double vol_1, vol_2, vol_3, vol_4, vol_5;
    vol_1 = SignedVolume(P_obs, A_, B_, C_);
    vol_2 = SignedVolume(P_max, A_, B_, C_);
    vol_3 = SignedVolume(P_obs, P_max, A_, B_);
    vol_4 = SignedVolume(P_obs, P_max, B_, C_);
    vol_5 = SignedVolume(P_obs, P_max, C_, A_);
    return((vol_1 * vol_2 < 0) && (vol_3 * vol_4 > 0) && (vol_4 * vol_5 > 0));
}

Point Triangle::ptIntersection(const Point & P1)
{
    Point P2(P1);
    P2[2] = 2*maxz();
    Point N((B_-A_) ^ (C_-A_));
    double t = - ((P1-A_) | N) / ((P2-P1) | N);
    return(P1 + t*(P2-P1));
}

Point Triangle::ptIntersection(const Point & P1, const Point & P2)
{
    Point N((B_-A_) ^ (C_-A_));
    double t = - ((P1-A_) | N) / ((P2-P1) | N);
    return(P1 + t*(P2-P1));
}

Point Triangle::ptIntersectionObs(const Point & P_obs, const Point & P_pixel)
{
    Point P_max((P_pixel - P_obs) * (1 + 2 * maxz()/(P_pixel[2] - P_obs[2])) + P_obs);
    Point N((B_-A_) ^ (C_-A_));
    double t = - ((P_obs-A_) | N) / ((P_max-P_obs) | N);
    return(P_obs + t*(P_max-P_obs));
}

Point Triangle::vectNormal()
{
    Point norm((A_-B_) ^ (B_-C_));
    if (norm_euc(norm) == 0.0) {
        cout << "Triangle plat -> vecteur normal nul" << endl;
        exit(1);
    }
    norm /= norm_euc(norm);
    return(norm);
}

double Triangle::aire()
{
    return(0.5 * norm_euc((B_-A_) ^ (C_-A_)));
}

std::vector<double> Triangle::coordBarycentric(const Point & P_intersection)
{
    Triangle T1(P_intersection, B_, C_);
    Triangle T2(A_, P_intersection, C_);
    Triangle T3(A_, B_, P_intersection);
    double aireT = aire();
    std::vector<double> vCoordBary(3);
    vCoordBary[0] = T1.aire()/aireT;
    vCoordBary[1] = T2.aire()/aireT;
    vCoordBary[2] = T3.aire()/aireT;
    return(vCoordBary);
}

vector<double> Triangle::angleOmbragePlat(const Point & P_lum)
{
    Point norm = vectNormal();
    vector<double> vAng(2);
    vAng[0] = acos((norm | (P_lum - bary_)) / norm_euc(P_lum - bary_)); //theta
    vAng[1] = acos((norm | (projBary_-bary_))/norm_euc(projBary_-bary_)) - vAng[0]; //alpha
    return(vAng);
}


