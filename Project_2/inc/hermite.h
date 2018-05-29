#ifndef HERMITE
#define HERMITE

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
 template<typename T> static T H5(const T& x) {return 32*pow(x, 5) - 160*pow(x, 3) + 120*x;}
 template<typename T> static T H6(const T& x) {return 64*pow(x, 6) - 480*pow(x, 4) + 720*pow(x, 2) - 120;}
 template<typename T> static T H7(const T& x) {return 128*pow(x, 7) - 1344*pow(x, 5) + 3360*pow(x, 3) - 1680*x;}
 template<typename T> static T H8(const T& x) {return 256*pow(x, 8) - 3584*pow(x, 6) + 13440*pow(x, 4) - 13440*pow(x, 2) + 1680;}
 template<typename T> static T H9(const T& x) {return 512*pow(x, 9) - 9216*pow(x, 7) + 48384*pow(x, 5) - 80640*pow(x, 3) + 30240*x;}
 template<typename T> static T H10(const T& x) {return 1024*pow(x, 10) - 23040*pow(x, 8) + 161280*pow(x, 6) - 403200*pow(x, 4) + 302400*pow(x, 2) - 30240;}
 template<typename T> static T H11(const T& x) {return 2048*pow(x, 11) - 56320*pow(x, 9) + 506880*pow(x, 7) - 1774080*pow(x, 5) + 2217600*pow(x, 3) - 665280*x;}
 template<typename T> static T H12(const T& x) {return 4096*pow(x, 12) - 135168*pow(x, 10) + 1520640*pow(x, 8) - 7096320*pow(x, 6) + 13305600*pow(x, 4) - 7983360*pow(x, 2) + 665280;}
 template<typename T> static T H13(const T& x) {return 8192*pow(x, 13) - 319488*pow(x, 11) + 4392960*pow(x, 9) - 26357760*pow(x, 7) + 69189120*pow(x, 5) - 69189120*pow(x, 3) + 17297280*x;}
 template<typename T> static T H14(const T& x) {return 16384*pow(x, 14) - 745472*pow(x, 12) + 12300288*pow(x, 10) - 92252160*pow(x, 8) + 322882560*pow(x, 6) - 484323840*pow(x, 4) + 242161920*pow(x, 2) - 17297280;}
 template<typename T> static T H15(const T& x) {return 32768*pow(x, 15) - 1720320*pow(x, 13) + 33546240*pow(x, 11) - 307507200*pow(x, 9) + 1383782400*pow(x, 7) - 2905943040*pow(x, 5) + 2421619200*pow(x, 3) - 518918400*x;}
 template<typename T> static T H16(const T& x) {return 65536*pow(x, 16) - 3932160*pow(x, 14) + 89456640*pow(x, 12) - 984023040*pow(x, 10) + 5535129600*pow(x, 8) - 15498362880*pow(x, 6) + 19372953600*pow(x, 4) - 8302694400*pow(x, 2) + 518918400;}
 template<typename T> static T H17(const T& x) {return 131072*pow(x, 17) - 8912896*pow(x, 15) + 233963520*pow(x, 13) - 3041525760*pow(x, 11) + 20910489600*pow(x, 9) - 75277762560*pow(x, 7) + 131736084480*pow(x, 5) - 94097203200*pow(x, 3) + 17643225600*x;}
#pragma GCC diagnostic pop
template<typename T> static T H(const T& x, const int& n) {
   if (n > 17) {
       return -1;
   }
   switch(n) {
       case 0: return H0(x);
       case 1: return H1(x);
       case 2: return H2(x);
       case 3: return H3(x);
       case 4: return H4(x);
       case 5: return H5(x);
       case 6: return H6(x);
       case 7: return H7(x);
       case 8: return H8(x);
       case 9: return H9(x);
       case 10: return H10(x);
       case 11: return H11(x);
       case 12: return H12(x);
       case 13: return H13(x);
       case 14: return H14(x);
       case 15: return H15(x);
       case 16: return H16(x);
       case 17: return H17(x);
       default: return 0;
   }
}
#endif /* HERMITE */
