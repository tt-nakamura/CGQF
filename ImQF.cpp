// uses NTL
//   http://www.shoup.net/ntl

#include "ImQF.h"
#include "HashMap.h"

ZZ ImQF::D;// discriminant
ZZ ImQF::h;// class number

long ImQF::init(const ZZ& d)
// set discriminant d<0, d=0,1(mod 4)
{
    if(sign(d)>=0) return 0;
    long k;
    k = trunc_long(d,2);
    if(k!=0 && k!=3) return 0;
    D = d;
    clear(h);
    return 1;
}

long ImQF::init(long d)
// set discriminant d<0, d=0,1(mod 4)
{
    if(d>=0 || d&2) return 0;
    D = d;
    clear(h);
    return 1;
}

void set_c(ImQF& f)
// c = (b^2 - D)/(4a)
{
    sqr(f.c, f.b);
    f.c -= ImQF::D;
    f.c /= f.a;
    f.c >>= 2;
}

void set(ImQF& f)
// f = principal form (a=1)
{
    set(f.a);
    negate(f.c, ImQF::D);
    if(IsOdd(f.c)) { set(f.b); f.c++; }
    else clear(f.b);
    f.c >>= 2;
}

void reduce(ImQF& f)
// f = reduced form of f
{
    long i;
    ZZ t,q,r;
    for(i=0;; i=1) {
        if(f.a > f.c) {
            negate(f.b, f.b);
            swap(f.a, f.c);
        }
        else if(i) break;
        negate(t, f.a);
        if(f.b <= f.a && f.b > t) break;
        LeftShift(t, f.a, 1);
        DivRem(q, r, f.b, t);
        if(r > f.a) { r-=t; q++; }
        add(t, f.b, r);
        t >>= 1;
        MulSubFrom(f.c, t, q);
        f.b = r;
    }
    if(f.a==f.c && sign(f.b)<0) negate(f.b, f.b);
}

void inv(ImQF& f, const ImQF& g)
// f = inverse of g
{
    if(&f!=&g) { f.a=g.a; f.c=g.c; }
    negate(f.b, g.b);
}

void SetPrime(ImQF& f, long p)
// f = prime form above p, Kronecker(D,p)>=0
{
    f.a = p;
    if(p==2) {
        if(IsOdd(ImQF::D)) set(f.b);// D=1 mod 8
        else if(NumTwos(ImQF::D) >= 3) clear(f.b);// D=0 mod 8
        else f.b=2;// D=4 mod 8
    }
    else {
        f.b = SqrRootMod(ImQF::D%p, p);
        if(IsOdd(f.b)^IsOdd(ImQF::D)) sub(f.b, p, f.b);
    }
    set_c(f);
}

void mul(ImQF& f, const ImQF& g, const ImQF& h)
// f = composition of g and h, in reduced form
// &f==&g is allowed
{
    if(IsOne(g)) { f=h; return; }
    if(IsOne(h)) { f=g; return; }
    ZZ s,n,d,u,v;
    XGCD(d,u,v, g.a, h.a);
    add(s, g.b, h.b);
    s>>=1;
    sub(n, g.b, s);
    u *= n;
    XGCD(d,s,v,s,d);
    div(n, h.a, d);
    u %= n;
    u *= v;
    MulAddTo(u, s, g.c);
    div(v, g.a, d);
    u %= n;
    u *= v;
    u <<= 1;
    sub(f.b, g.b, u);
    mul(f.a, n, v);
    set_c(f);
    reduce(f);
}

void sqr(ImQF& f, const ImQF& g)
// f = composition of g and g, in reduced form
{
    if(divide(g.b, g.a)) { set(f); return; }
    ZZ x,y,d;
    XGCD(d,x,y, g.b, g.a);
    x *= g.c;
    div(f.a, g.a, d);
    x %= f.a;
    x *= f.a;
    x <<= 1;
    sub(f.b, g.b, x);
    sqr(f.a, f.a);
    set_c(f);
    reduce(f);
}

void power(ImQF& f, const ImQF& g, long n)
// f = g^n, in reduced form
{
    if(n==0) { set(f); return; }
    if(n<0) { inv(f,g); power(f,f,-n); return; }
    if(&f==&g) { ImQF h(g); power(f,h,n); return; }
    long m(1L<<(NumBits(n)-1));
    f=g;
    for(m>>=1; m; m>>=1) {
        sqr(f,f);;
        if(n&m) mul(f,f,g);
    }
}

void power(ImQF& f, const ImQF& g, const ZZ& n)
// f = g^n, in reduced form
{
    if(IsZero(n)) { set(f); return; }
    if(sign(n)<0) { inv(f,g); power(f,f,-n); return; }
    if(&f==&g) { ImQF h(g); power(f,h,n); return; }
    long i;
    f=g;
    for(i=NumBits(n)-2; i>=0; i--) {
        sqr(f,f);;
        if(bit(n,i)) mul(f,f,g);
    }
}

std::ostream& operator<<(std::ostream& s, const ImQF& f) {
    s << '[' << f.a << ' ' << f.b;
    s << ' ' << f.c << ']';
    return s;
}

const ZZ& ImQF::ClassNum(long P)
// return class number h
// If h was computed previously, h is not recomputed.
// Otherwise, h is computed only APPROXIMATELY
//   using Dirichlet formula with upper bound P.
//   (default P = 2^{18})
{
    if(!IsZero(h)) return h;
    if(D>-15) { set(h); return h; }
    static double PI(4*atan(1));
    long p,k;
    double x;
    conv(x,D);
    x = sqrt(-x)/PI;
    PrimeSeq ps;
    while((p = ps.next()) <= P) {
        k = Kronecker(D,p);
        if     (k>0) x /= (1 - 1./p);
        else if(k<0) x /= (1 + 1./p);
    }
    conv(h, round(x));
    return h;
}

void order(ZZ& n, const ImQF& f)
// n = order of f using Terr method
{
    long j,k;
    HashMap<ImQF, long> B(NumBits(ImQF::D)>>1);
    ImQF b(f), g(f), e;
    set(e);
    B.install(e,0);
    set(n);
    for(j=1; !B.lookup(g,k); j++, n+=j) {
        B.install(b,j);
        b *= f;
        g *= b;
    }
    n -= k;
}
