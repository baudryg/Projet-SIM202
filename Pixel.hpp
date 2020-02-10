#ifndef PIXEL_H
#define PIXEL_H

#include "Maillage.hpp"
#include "Triangle.hpp"
#include "Pixel.hpp"

class Pixel
{
public:
    uint8_t red_, green_, blue_;


public:
    Pixel() : red_(0), green_(0), blue_(0) {}
    Pixel(uint8_t red, uint8_t green, uint8_t blue) : red_(red), green_(green), blue_(blue) {}

    uint8_t operator[](int i) const;
    uint8_t & operator[](int);
};

#endif // M
