// uses NTL
//   http://www.shoup.net/ntl

#ifndef __ImQF_h__
#define __ImQF_h__

#include<NTL/ZZ.h>
#include<NTL/pair.h>
using namespace NTL;

struct ImQF {// Imaginary Quadratic Form
    ZZ a,b,c;// ax^2 + bxy + cy^2
    static ZZ D,h;
    static long init(const ZZ& D);
    static long init(long D);
    static const ZZ& ClassNum(long=1<<18);
    inline long hash_val() const { return conv<long>(b) ^ conv<long>(c); }
};

void set_c(ImQF&);
void set(ImQF&);
void SetPrime(ImQF&, long);
void reduce(ImQF&);
void inv(ImQF&, const ImQF&);

void mul(ImQF&, const ImQF&, const ImQF&);
void sqr(ImQF&, const ImQF&);
void power(ImQF&, const ImQF&, long);
void power(ImQF&, const ImQF&, const ZZ&);
void order(ZZ&, const ImQF&);
std::ostream& operator<<(std::ostream&, const ImQF&);

inline bool IsOne(const ImQF& f) { return IsOne(f.a); }
inline bool operator==(const ImQF& f, const ImQF& g) { return (f.a==g.a && f.b==g.b); }
inline void operator*=(ImQF& f, const ImQF& g) { mul(f,f,g); }

void ClassGroup(Vec<Pair<ImQF,ZZ> >&);
void ClassGroup2(Vec<Pair<ImQF,ZZ> >&, double=1);

long Kronecker(const ZZ&, long);
long SqrRootMod(long, long);

#endif // __ImQF_h__
