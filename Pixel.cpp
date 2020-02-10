#include "Pixel.hpp"
#include <iostream>
#include <cstdlib>

using namespace std;

uint8_t Pixel::operator[](int i) const
{
    switch (i)
      {
         case 0:
            return(red_);
         case 1:
            return(green_);
         case 2:
            return(blue_);
         default:
            cout << "Indice different de 1, 2 ou 3." << endl;
            exit(1);
      }
}

uint8_t & Pixel::operator[](int i)
{
    switch (i)
      {
         case 0:
            return(red_);
         case 1:
            return(green_);
         case 2:
            return(blue_);
         default:
            cout << "Indice different de 1, 2 ou 3." << endl;
            exit(1);
      }
}
