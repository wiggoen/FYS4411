#ifndef HERMITE_H
#define HERMITE_H

#include <vector>
#include <cmath>
#include <array>


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
 template<typename T> static T H0(const T& x) {return 1;}
 template<typename T> static T H1(const T& x) {return 2*x;}
 template<typename T> static T H2(const T& x) {return 4*pow(x, 2) - 2;}
 template<typename T> static T H3(const T& x) {return 8*pow(x, 3) - 12*x;}
 template<typename T> static T H4(const T& x) {return 16*pow(x, 4) - 48*pow(x, 2) + 12;}
#pragma GCC diagnostic pop
template<typename T> static T H(const T& x, const int& n) {
   if (n > 4) {
       return -1;
   }
   switch(n) {
       case 0: return H0(x);
       case 1: return H1(x);
       case 2: return H2(x);
       case 3: return H3(x);
       case 4: return H4(x);
       default: return 0;
   }
}
#endif /* HERMITE_H */
