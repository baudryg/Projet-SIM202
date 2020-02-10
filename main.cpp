#include <iostream>
#include "Maillage.hpp"
#include "Image.hpp"
#include "bitmap_image.hpp"
#include"Ray_Tracing.hpp"

#include "Triangle.hpp"

using namespace std;

int main()
{
    /*Point P1(1.0, 1.0, 1.0);
    Point P2(1.0, 2.0, 2.0);
    Point P3(0.0, 1.0, 2.0);
    Point P4(2.0, 1.0, 2.0);
    Point P5(1.0, 0.0, 2.0);

    Triangle T1(P1, P2, P3);
    Triangle T2(P1, P4, P2);
    Triangle T3(P1, P5, P4);
    Triangle T4(P1, P3, P5);

    Point P_lum(0.0, 3.0, 0.0);
    Point P_obs(1.0, 1.0, -2.0);

    cout << T1.angleOmbragePlat(P_lum)[0] << " " << T1.angleOmbragePlat(P_lum)[1] << endl;
    cout << T2.angleOmbragePlat(P_lum)[0] << " " << T2.angleOmbragePlat(P_lum)[1] << endl;
    cout << T3.angleOmbragePlat(P_lum)[0] << " " << T3.angleOmbragePlat(P_lum)[1] << endl;
    cout << T4.angleOmbragePlat(P_lum)[0] << " " << T4.angleOmbragePlat(P_lum)[1] << endl;

    vector<Point> vPts(5);
    vPts[0] = P1; vPts[1] = P2; vPts[2] = P3; vPts[3] = P4; vPts[4] = P5;

    vector<Triangle> vTris(4);
    vTris[0] = T1; vTris[1] = T2; vTris[2] = T3; vTris[3] = T4;

    Maillage M(vPts, vTris);

    Point P_pixel(0.9, 1.1, 0.0);
    Point P_test(0.7, 1.3, 4.0);

    Point Pi(0.9, 1.1, 0.0);
    Point Pj(-0.1, 0.1, -1.0);
    Point Pk(0.1, 0.9, 2.0);
    Point Pl = Pj ^ Pk;
    Pl.print();
    double ps = Pi | Pl;
    cout << "(Pi | Pl) = " << ps << endl;
    cout << "Is T1 intersected by P_obs --- P_pixel ? " << T1.isIntersectedObs(P_obs, P_pixel) << endl;
    cout << "Is T1 intersected by Pi ? " << T1.isIntersected(P_obs, P_test) << endl;
    Point P_intersection = T1.ptIntersectionObs(P_obs, P_pixel);
    cout << "Intersection worked ? " << T1.isIntersectedObs(P_intersection, Pk) << endl;
    T1.vectNormal().print();*/

    Image img(2.0, 2.0, 300, 300);
    Point P_src(0.0, 3.0, 0.0);
    Point P_observ(1.0, 1.0, -2.0);

    Maillage M;
    M.loadGMSH("Fichiers_GMSH/test_sphere.msh");

    M.vectNormPoint(0).print();

    RayTracing RT(M, img, P_src, P_observ);
    RT.afficherImage();

    /*cout << "Nb de points : " << M.gNbPt() << endl;
    cout << "Nb de triangles : " << M.gNbTri() << endl;
    for (int i=0; i<40; ++i)
    {
        cout << "Num : " << i <<endl;
        M[i].print();
        cout << endl;
    }*/

    return 0;
}
