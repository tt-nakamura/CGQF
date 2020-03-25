// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/vec_RR.h>
using namespace NTL;

static RR y,z;

void RGCD(RR& x, const RR& a, const RR& b, double eps)
// output: x = GCD of real numbers a and b
//         such that a/x and b/x are both integers
// input: eps = lower bound of x
{
    abs(x,a);
    abs(y,b);
    while(y>=eps) {
        div(z,x,y);
        floor(z,z);
        z *= y;
        sub(z,x,z);
        x = y;
        y = z;
    }
}

void RGCD(RR& x, const vec_RR& v, double eps)
// output: x = GCD of array of real numbers v
// input: eps = lower bound of x
{
    if(v.length()==0) { clear(x); return; }
    long i;
    abs(x, v[0]);
    for(i=1; i<v.length(); i++)
        RGCD(x, x, v[i], eps);
}
