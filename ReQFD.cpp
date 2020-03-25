// uses NTL
//   http://www.shoup.net/ntl

#include "ReQF.h"
#include "HashTab.h"

static RR x,y;

void conv(ReQFD& f, const ReQF& g)
// distance of f = 1
{
    f.a = g.a;
    f.b = g.b;
    f.c = g.c;
    set(f.d);
}

void set(ReQFD& f)
// f = principal form (a=1), with distance = 1
{
    set((ReQF&)f);
    set(f.d);
}

void dist(RR& d, const ReQF& f)
// d = distance of f
{
    conv(d, f.b); d += ReQF::SD;
    sqr(d,d);
    conv(x, f.a); d /= x;
    conv(x, f.c); d /= x;
    d /= 4;
}

void rho(ReQFD& f)
// continued fraction expansion of f
//   with distance updated
{
    dist(y,f);
    f.d /= y;
    rho((ReQF&)f);
}

void rho_inv(ReQFD& f)
// inverse continued fraction expansion of f
//   with distance updated
{
    rho_inv((ReQF&)f);
    dist(y,f);
    f.d *= y;
}

void mul(ReQFD& f, const ReQFD& g, const ReQFD& h)
// f = composition of g and h, with distance updated
// &f==&g is allowd
{
    mul((ReQF&)f,g,h);
    mul(f.d, g.d, h.d);
}

void div(ReQFD& f, const ReQFD& g, const ReQFD& h)
// f = composition of g and h^{-1}, with distance updated
// &f==&g is allowd
{
    div((ReQF&)f,g,h);
    div(f.d, g.d, h.d);
}

void close(ReQFD& f)
// reduce f to closest form, with distance updated
// referece:
//   J. Buchmann and U. Vollmer
//    "Binary Quadratic Forms: An Algorithmic Approach"
//     Algorithm 10.1
{
    ZZ u,v,t;
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
    else for(ReQFD g(f); sign(t) + sign(v); f=g) {
        rho(g);
        t = u;
        u *= ReQF::q;
        u -= v;
        negate(v,t);
        mul(t, u, g.a); t<<=1;
        MulAddTo(t, v, g.b);
    }
}

const RR& ReQF::regulator()
// return:
//   regulator R
//   If R is computed previously, R is not recomputed
//   Otherwise R is computed by Terr method
// referece:
//   J. Buchmann and U. Vollmer
//    "Binary Quadratic Forms: An Algorithmic Approach"
//     Algorithm 10.2
{
    if(!IsZero(R)) return R;
    long i;
    ReQFD f,g,*h;
    HashTable<ReQFD> H(NumBits(S));
    set(f);
    set(g);
    for(i=0; i<LD; i++) {
        rho_inv(f);
        if(IsOne(f)) break;
        H.install(f);
    }
    while(f!=g) {
        div(g,g,f);
        close(g);
        if(h = H.lookup(g)) {
            f.d = h->d;
            break;
        }
        rho_inv(f);
        H.install(f);
    }
    f.d /= g.d;
    E = sign(f.d);
    abs(R, f.d);
    log(R, R);
    R /= 2;
    return R;
}
