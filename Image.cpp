#include "Image.hpp"
#include <vector>
#include "Triangle.hpp"
#include "Pixel.hpp"
#include "bitmap_image.hpp"

using namespace std;

Image::Image(double lx, double ly, int nbx, int nby) : lx_(lx), ly_(ly), nbx_(nbx), nby_(nby)
{
    vPixel_.resize(nbx);
    for (int i=0; i<nbx; ++i) {
        vPixel_[i].resize(nby);
    }
}

void Image::displayImageBMP()
{
    bitmap_image imgBMP(nbx_, nby_);
    imgBMP.clear();

    for (unsigned int y = 0; y < imgBMP.height(); ++y)
    {
       for (unsigned int x = 0; x < imgBMP.width(); ++x)
       {
            imgBMP.set_pixel(x, y, (vPixel_[x][y])[0], (vPixel_[x][y])[1], (vPixel_[x][y])[2]);
       }
    }
    imgBMP.set_pixel(0, 0, 255, 0, 0);

    imgBMP.save_image("test.bmp");
}
