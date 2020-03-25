// uses NTL
//   http://www.shoup.net/ntl

#ifndef __ReQF_h__
#define __ReQF_h__

#include<NTL/ZZ.h>
#include<NTL/RR.h>
#include<NTL/pair.h>
using namespace NTL;

struct ReQF {// Real Quadratic Form
    ZZ a,b,c;// ax^2 + bxy + cy^2
    static ZZ D,S,h,q;
    static RR SD,R;
    static long LD,E;
    static long init(const ZZ&);
    static long init(long);
    static const ZZ& ClassNum();
    static const RR& regulator();
    inline long hash_val() const { return conv<long>(b) ^ conv<long>(c); }
};

struct ReQFD : ReQF { RR d; };// d = distance

void set(ReQF&);
void set_c(ReQF&);
void SetPrime(ReQF&, long);
void rho(ReQF&);
void rho_inv(ReQF&);
void normalize(ReQF&);
void inv(ReQF&, const ReQF&);
long IsReduced(const ReQF&);

void mul(ReQF&, const ReQF&, const ReQF&);
void sqr(ReQF&, const ReQF&);
void div(ReQF&, const ReQF&, const ReQF&);
void PowerRed(ReQF&, const ReQF&, const ZZ&);
void PowerRed(ReQF&, const ReQF&, long);
long IsEquiv(const ReQF&, const ReQF&);
void close(ReQF&);
std::ostream& operator<<(std::ostream&, const ReQF&);

void dist(RR&, const ReQF&);
void conv(ReQFD&, const ReQF&);
void set(ReQFD&);
void rho(ReQFD&);
void rho_inv(ReQFD&);
void mul(ReQFD&, const ReQFD&, const ReQFD&);
void div(ReQFD&, const ReQFD&, const ReQFD&);
void close(ReQFD&);

inline long operator==(const ReQF& f, const ReQF& g) { return f.a==g.a && f.b==g.b; }
inline long operator!=(const ReQF& f, const ReQF& g) { return f.a!=g.a || f.b!=g.b; }
inline long IsOne(const ReQF& f) { return IsOne(f.a); }
inline void reduce(ReQF& f)  { while(!IsReduced(f)) rho(f); }
inline void reduce(ReQFD& f) { while(!IsReduced(f)) rho(f); }

#define MulRed(f,g,h) { mul(f,g,h); reduce(f); }
#define DivRed(f,g,h) { div(f,g,h); reduce(f); }
#define SqrRed(f,g) { sqr(f,g); reduce(f); }

void ClassGroup(Vec<Pair<ReQF,ZZ> >&);
void ClassGroup2(Vec<Pair<ReQF,ZZ> >&, double=1);

long Kronecker(const ZZ&, long);
long SqrRootMod(long, long);

#endif // __ReQF_h__
