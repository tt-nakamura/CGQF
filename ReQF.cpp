// uses NTL
//   http://www.shoup.net/ntl

#include "ReQF.h"
#include "HashTab.h"

ZZ ReQF::D;// discriminant
ZZ ReQF::h;// class number
ZZ ReQF::S;// floor[sqrt(D)]
ZZ ReQF::q;// b = 2aq + r
RR ReQF::R;// regulator = log(fundamental unit)
RR ReQF::SD;// sqrt(D)
long ReQF::LD;// ceil[log_2(D)]
long ReQF::E;// norm of fundamental unit

static ZZ t,r;

long ReQF::init(const ZZ& d)
// set discriminant, d>0, d=0,1(mod 4)
{
    if(d<=4 || bit(d,1)) return 0;
    SqrRoot(S,d);
    sqr(D,S);
    if(D==d) return 0;
    D=d;
    conv(SD,D);
    SqrRoot(SD,SD);
    LD = long(ceil(log(D)/log(2.)));
    if(LD&1) LD++;
    clear(h);
    clear(R);
    E=0;
    return 1;
}

long ReQF::init(long d)
// set discriminant, d>0, d=0,1(mod 4)
{ return init(t=d); }

void set(ReQF& f)
// f = principal form (a=1)
{
    set(f.a);
    f.b = ReQF::S;
    if(IsOdd(ReQF::D)^IsOdd(f.b)) f.b--;
    set_c(f);
}

void set_c(ReQF& f)
// c = (b^2 - D)/(4a)
{
    sqr(f.c, f.b);
    f.c -= ReQF::D;
    f.c >>= 2;
    f.c /= f.a;
}

long IsReduced(const ReQF& f) {
    if(sign(f.b) <= 0 || f.b > ReQF::S) return 0;
    LeftShift(t, f.a, 1);
    t -= ReQF::S;
    if(t == f.b) return 1;
    abs(t,t);
    return t < f.b;
}

void normalize(ReQF& f) {
    LeftShift(t, f.a, 1);
    sub(f.b, ReQF::S, f.b);
    sub(f.b, ReQF::S, f.b%=t);
    set_c(f);
}

void inv(ReQF& f, const ReQF& g)
// f = inverse of g
{
    if(&f!=&g) {
        f.a = g.a;
        f.c = g.c;
    }
    negate(f.b, g.b);
}

void rho_inv(ReQF& f)
// inverse continued fraction expansion of f
{
    LeftShift(t, f.a, 1);
    add(r, f.b, ReQF::S);
    DivRem(ReQF::q, r, r, t);
    if(sign(f.a)<0 && IsZero(r)) {
        ReQF::q--;
        r=t;
    }
    sub(r, ReQF::S, r);
    f.b -= r;
    f.b >>= 1;
    MulSubFrom(f.c, f.b, ReQF::q);
    swap(f.a, f.c);
    negate(f.a, f.a);
    negate(f.c, f.c);
    f.b = r;
}

void rho(ReQF& f)
// continued fraction expansion of f
{
    swap(f.a, f.c);
    negate(f.a, f.a);
    negate(f.c, f.c);
    LeftShift(t, f.a, 1);
    add(r, f.b, ReQF::S);
    DivRem(ReQF::q, r, r, t);
    if(sign(f.a)<0 && IsZero(r)) {
        ReQF::q--;
        r=t;
    }
    sub(r, ReQF::S, r);
    f.b -= r;
    f.b >>= 1;
    MulSubFrom(f.c, f.b, ReQF::q);
    f.b = r;
}

void SetPrime(ReQF& f, long p)
// f = prime form above p, Kronecker(D,p)>=0
{
    f.a = p;
    if(p==2) {
        if(IsOdd(ReQF::D)) set(f.b);// D=1 mod 8
        else if(NumTwos(ReQF::D)>=3) clear(f.b);// D=0 mod 8
        else f.b = 2;// D=4 mod 8
    }
    else {
        f.b = SqrRootMod(ReQF::D%p, p);
        if(IsOdd(f.b)^IsOdd(ReQF::D)) sub(f.b, p, f.b);
    }
    normalize(f);
}

void mul(ReQF& f, const ReQF& g, const ReQF& h)
// f = composition of g and h (not reduced)
// &f==&g is allowed
{
    if(IsOne(g)) { f=h; return; }
    if(IsOne(h)) { f=g; return; }
    ZZ s,n,d,u,v;
    XGCD(d,u,v,g.a,h.a);
    add(s, g.b, h.b); s>>=1;
    sub(n, g.b, s); u *= n;
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
}

void sqr(ReQF& f, const ReQF& g)
// f = composition of g and g (not reduced)
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
}

void div(ReQF& f, const ReQF& g, const ReQF& h)
// f = composition of g and h^{-1} (not reduced)
{
    ReQF& h1(const_cast<ReQF&>(h));
    negate(h1.b, h1.b);
    mul(f,g,h1);
    negate(h1.b, h1.b);
}

void PowerRed(ReQF& f, const ReQF& g, long n) {
// f = g^n, in reduced form
    if(n==0) { set(f); return; }
    if(n<0) { inv(f,g); PowerRed(f,f,-n); return; }
    if(&f==&g) { ReQF h(g); PowerRed(f,h,n); return; }
    long m(1L<<(NumBits(n)-1));
    f=g;
    for(m>>=1; m; m>>=1) {
        SqrRed(f,f);
        if(n&m) MulRed(f,f,g);
    }
}

void PowerRed(ReQF& f, const ReQF& g, const ZZ& n)
// f = g^n, in reduced form
{
    if(IsZero(n)) { set(f); return; }
    if(sign(n)<0) { inv(f,g); PowerRed(f,f,-n); return; }
    if(&f==&g) { ReQF h(g); PowerRed(f,h,n); return; }
    long i;
    f=g;
    for(i=NumBits(n)-2; i>=0; i--) {
        SqrRed(f,f);
        if(bit(n,i)) MulRed(f,f,g);
    }
}

std::ostream& operator<<(std::ostream& s, const ReQF& f) {
    s << '[' << f.a << ' ' << f.b;
    s << ' ' << f.c << ']';
    return s;
}

void close(ReQF& f)
// reduce f to closest form
// referece:
//   J. Buchmann and U. Vollmer
//    "Binary Quadratic Forms: An Algorithmic Approach"
//     Algorithm 10.1
{
    ZZ u,v;
    set(u);
    while(!IsReduced(f)) {
        rho_inv(f);
        t = v;
        v *= ReQF::q;
        v += u;
        u = t;
    }
    mul(t, u, f.a); t<<=1;
    MulAddTo(t, v, f.b);
    if(sign(t) + sign(v) == 0) do {
        rho_inv(f);
        t = v;
        v *= ReQF::q;
        v += u;
        u = t;
        mul(t, u, f.a); t<<=1;
        MulAddTo(t, v, f.b);
    } while(sign(t) + sign(v) == 0);
    else for(ReQF g(f); sign(t) + sign(v); f=g) {
        rho(g);
        t = u;
        u *= ReQF::q;
        u -= v;
        negate(v,t);
        mul(t, u, g.a); t<<=1;
        MulAddTo(t, v, g.b);
    }
}

long IsEquiv(const ReQF& a, const ReQF& b)
// test if a and b are equivalent or not
// reference:
//   J. Buchmann and U. Vollmer
//    "Binary Quadratic Forms: An Algorithmic Approach"
//     Algorithm 10.3
{
    long i,j,n;
    double R;
    ReQF f,g;
    DivRed(g,a,b);
    if(IsOne(g)) return 1;
    HashTable<ReQF> H(NumBits(ReQF::S));
    set(f);
    for(i=0; i<ReQF::LD; i++) {
        rho_inv(f);
        if(f==g) return 1;
    }
    conv(R, ReQF::regulator());
    n = long(ceil(sqrt(2*R/log(2.)))) + 1;
    for(i=0; i<n; i++) {
        div(g,g,f);
        close(g);
        if(H.lookup(g)) return 1;
        for(j=0; j<2; j++) {
            rho_inv(f);
            if(f==g) return 1;
            H.install(f);
        }
    }
    return 0;
}

const ZZ& ReQF::ClassNum()
// return class number h
// If h was computed previously, h is not recomputed.
// Otherwise, h is computed only APPROXIMATELY
//   using Dirichlet formula.
{
    if(!IsZero(h)) return h;
    if(D<40) { set(h); return h; }
    long n(500),m(6),i,j,k,p;
    double x,y,eps(1.e-1);
    PrimeSeq ps;
    conv(x, SD);
    conv(y, regulator());
    x /= 2.*y;
    for(j=0; j<m;) {
        for(i=0; i<n; i++) {
            p = ps.next();
            k = Kronecker(D,p);
            if     (k>0) x /= (1 - 1./p);
            else if(k<0) x /= (1 + 1./p);
        }
        y = round(x);
        if(fabs(x-y) <= eps) j++; else j=0;
    }
    conv(h,y);
    return h;
}
