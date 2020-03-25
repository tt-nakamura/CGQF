// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/ZZ.h>
using namespace NTL;

long Jacobi(long a, long b)
// Jacobi symbol (a/b), 0<=a<b, b odd
{
    long c,t(0);
    while(a) {
        for(c=0; (a&1)==0; c++) a>>=1;
        if(c&1) {
            c = b&7;
            if(c==3 || c==5) t ^= 1;
        }
        if(a&b&2) t ^= 1;
        c = b%a;
        b = a;
        a = c;
    }
    if(b!=1) return 0;
    else if(t) return -1;
    else return 1;
}

long Kronecker(const ZZ& a, long b)
// Kronecker symbol (a/b), b>0
{
    if(b&1) return Jacobi(a%b, b);
    if(!IsOdd(a)) return 0;
    long c,j;
    for(c=0; (b&1)==0; c++) b>>=1;
    j = Jacobi(a%b, b);
    b = trunc_long(a,3);
    if(c&1 && (b==3 || b==5)) return -j;
    else return j;
}

long SqrRootMod(long a, long p)
// input:
//   a = integer, 0 <= a < p
//   p = odd prime, p < 2^{30}
// return:
//   x such that x^2 = a (mod p)
//   by Cipolla method
{
    if(a==0) return 0;
    long d, s((p+1)>>1);
    long b,c0(1),c1(0),x0,x1(1);
    do {
        b = rand()%p;
        d = (b*b-a)%p;
        if(d<0) d+=p;
    } while(Jacobi(d,p) >= 0);
    for(x0=b;;) {
        if(s&1) {
            b = (c0*x0 + c1*x1%p*d)%p;
            c1 = (c0*x1 + c1*x0)%p;
            c0 = b;
        }
        s>>=1;
        if(s==0) return b;
        b = (x0*x0 + x1*x1%p*d)%p;
        x1 = (x0*x1<<1)%p;
        x0 = b;
    }
}
