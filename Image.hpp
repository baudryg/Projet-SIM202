#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include "Pixel.hpp"

using namespace std;

class Image
{
public:
    double lx_, ly_;
    int nbx_, nby_;
    vector<vector<Pixel> > vPixel_; // ATTENTION : tout est en public



public:
    Image(double, double, int, int);

    int gnbx() {return(nbx_);}
    int gnby() {return(nby_);}
    int glx() {return(lx_);}
    int gly() {return(ly_);}

    void displayImageBMP();
};

#endif // IMAGE_H
