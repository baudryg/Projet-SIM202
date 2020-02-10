#include "Ray_Tracing.hpp"
#include <cmath>
#include "Triangle.hpp"
#include <iostream>


double RayTracing::intensiteFromAngle(double theta, double alpha)
{
    return(Ia*ka + Is * (kd * cos(theta) + kr * pow(cos(alpha), n)));
}

double RayTracing::intensiteFromAngle()
{
    return(Ia*ka);
}

double RayTracing::intensiteGouraud(int id_T, const Point & P_intersection, const Point & P_pixel)
{
    double intensite;
    vector<double> coordBary = maillageRT_(id_T).coordBarycentric(P_intersection);
    vector<double> vAng;
    std::vector<int> vPointAroundTriangle = maillageRT_.getPointAroundTriangle(id_T);
    for (int i=0; i<3; ++i) {
        vector<double> vAng = maillageRT_.angleOmbrage(vPointAroundTriangle[i], srcLumRT_, P_pixel);
        intensite += coordBary[i] * intensiteFromAngle(vAng[0], vAng[1]);
    }
    return(intensite);
}

void RayTracing::colorImagePlat()
{
    for (int i=0; i<imgRT_.gnbx(); ++i) {
        if (i%2 == 0)
        {
            cout << i/2 << "%" << endl;
        }
        for (int j=0; j<imgRT_.gnby(); ++j) {
            Point Pij(imgRT_.glx()/double(imgRT_.gnbx())*(i+0.5), imgRT_.gly()/double(imgRT_.gnby())*(imgRT_.gnby() - 1 - j + 0.5), 0.0);

            if (maillageRT_.crossTriangleObs(obsRT_, Pij)) {
                int id_Tij = maillageRT_.closestIntersectedTriangleIdObs(obsRT_, Pij);
                Triangle Tij = maillageRT_(id_Tij);
                Point P_intersection = Tij.ptIntersectionObs(obsRT_, Pij);
                if (maillageRT_.crossOtherTriangle(id_Tij, P_intersection, srcLumRT_)) {
                    for (int k=0; k<3; ++k) {
                        (imgRT_.vPixel_[i][j])[k] = uint8_t(intensiteFromAngle()*255);
                    }

                    // Affichage
                    /*cout << "i = " << i << " j = " << j << endl;
                    Tij.print();
                    cout << "intensite : " << intensiteFromAngle() << " ";
                    Pij.print();*/
                }
                else {
                    vector<double> Angij = Tij.angleOmbragePlat(srcLumRT_);
                    for (int k=0; k<3; ++k) {
                        (imgRT_.vPixel_[i][j])[k] = uint8_t(intensiteFromAngle(Angij[0], Angij[1])*255);
                    }

                    //Affichage
                    /*cout << "i = " << i << " j = " << j << endl;
                    Tij.print();
                    cout << "theta = " << Angij[0] << "alpha = " << Angij[1] << endl;
                    cout << "intensite : " << intensiteFromAngle(Angij[0], Angij[1]) << " ";
                    Pij.print();*/
                }
            }
            else {
                for (int k=0; k<3; ++k) {
                    (imgRT_.vPixel_[i][j])[k] = uint8_t(0);
                }
                /*cout <<"pas de croisement : ";
                Pij.print();*/
            }
        }
    }
}

void RayTracing::colorImageGouraud()
{
    for (int i=0; i<imgRT_.gnbx(); ++i) {
        for (int j=0; j<imgRT_.gnby(); ++j) {

            Point Pij(imgRT_.glx()/double(imgRT_.gnbx())*(i+0.5), imgRT_.gly()/double(imgRT_.gnby())*(imgRT_.gnby() - 1 - j + 0.5), 0.0);

            if (maillageRT_.crossTriangleObs(obsRT_, Pij)) {
                int id_Tij = maillageRT_.closestIntersectedTriangleIdObs(obsRT_, Pij);
                Triangle Tij = maillageRT_(id_Tij);
                Point P_intersection = Tij.ptIntersectionObs(obsRT_, Pij);

                if (maillageRT_.crossOtherTriangle(id_Tij, P_intersection, srcLumRT_)) {

                }
                else {
                    PointAroundTriangle = maillageRT_.getPointAroundTriangle(id_Tij);

                }
            }
            else {
                for (int k=0; k<3; ++k) {
                    (imgRT_.vPixel_[i][j])[k] = uint8_t(0);
                }
                /*cout <<"pas de croisement : ";
                Pij.print();*/
            }
        }
    }
}

void RayTracing::afficherImage()
{
    colorImageGouraud();
    imgRT_.displayImageBMP();
}


