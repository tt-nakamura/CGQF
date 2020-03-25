// uses NTL
//   http://www.shoup.net/ntl

#include "ReQF.h"
#include<NTL/mat_ZZ.h>

void SmithNF(vec_ZZ&, mat_ZZ&, const mat_ZZ&, const ZZ&);

void HNF2SNF(Vec<Pair<ReQF,ZZ> >& G, const Vec<ReQF>& F, const mat_ZZ& U)
// convert Hermite normal form to Smith normal form
{
    long i,j,k,n(F.length());
    vec_ZZ d;
    mat_ZZ V;
    ReQF f;
    SmithNF(d,V,U, ReQF::h);
    for(k=d.length(); k; k--)
        if(!IsOne(d[k-1])) break;
    G.SetLength(k);
    for(i=0; i<k; i++) {
        set(G[i].a);
        for(j=0; j<n; j++) {
            PowerRed(f, F[j], V[i][j]);
            MulRed(G[i].a, G[i].a, f);
        }
        G[i].b = d[i];
    }
}    

void ClassGroup(Vec<Pair<ReQF,ZZ> >& G)
// output:
//   G = vector of ("generator", "order") pair
//       in decreasing order of "order"
//       by Terr-Shanks method
{
    long i,j,k,l,m,p,q;
    ZZ a,h1,h2;
    ZZ& h(ReQF::h);
    Vec<ReQF> F;
    vec_vec_ZZ U;
    mat_ZZ U1;
    ReQF f,g,f1,g1;
    Vec<Pair<ReQF, Vec<long> > > B,B1;

    G.SetLength(0);
    ReQF::ClassNum();
    if(IsOne(h)) return;
    RightShift(h2,h,2);
    sub(h1,h,h2);
    add(h2,h,h2);
    B.SetLength(1);
    set(B[0].a);
    PrimeSeq ps;
    set(h);
    for(l=0;;) {
        p = ps.next();
        if(Kronecker(ReQF::D, p) <= 0) continue;
        SetPrime(f1,p);
        reduce(f1);
        f = f1;
        g = f1;
        B1 = B;
        k = B.length();
        for(i=0; i<k; i++) B1[i].b.append(0);
        for(p=q=1;; q++, p+=q) {
            for(j=0; j<k; j++)
                if(IsEquiv(g, B1[j].a)) goto a;
            for(i=0; i<B.length(); i++, k++) {
                B1.SetLength(k+1);
                MulRed(B1[k].a, B[i].a, f);
                B1[k].b = B[i].b;
                B1[k].b.append(q);
            }
            MulRed(f,f,f1);
            MulRed(g,g,f);
        }
    a:  ;
        p -= B1[j].b[l];
        if(p==1) continue;
        mul(a,h,p);
        if(a>h2) continue;
        F.append(f1);
        U.SetLength(l+1);
        U[l].SetLength(l+1);
        for(k=0; k<l; k++) U[l][k] = B1[j].b[k];
        U[l][l] = p;
        if((h=a) >= h1) break;
        inv(f1,f1);
        f = f1;
        k = m = B.length();
        for(i=1; i<p; i++) {
            for(j=0; j<m; j++, k++) {
                B.SetLength(k+1);
                MulRed(B[k].a, B[j].a, f);
                B[k].b = B[j].b;
                B[k].b.append(i);
            }
            MulRed(f,f,f1);
        }
        for(i=0; i<m; i++) B[i].b.append(0);
        l++;
    }
    l = F.length();
    U1.SetDims(l,l);
    for(i=0; i<l; i++) VectorCopy(U1[i], U[i], l);
    HNF2SNF(G,F,U1);
}
