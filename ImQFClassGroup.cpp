// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/mat_ZZ.h>
#include "ImQF.h"
#include "HashMap.h"

void SmithNF(vec_ZZ&, mat_ZZ&, const mat_ZZ&, const ZZ&);

void HNF2SNF(Vec<Pair<ImQF,ZZ> >& G, const Vec<ImQF>& F, const mat_ZZ& U)
// convert Hermite norm form to Smith normal form
{
    long i,j,k,n(F.length());
    vec_ZZ d;
    mat_ZZ V;
    ImQF f;
    SmithNF(d,V,U, ImQF::h);
    for(k=d.length(); k; k--)
        if(!IsOne(d[k-1])) break;
    G.SetLength(k);
    for(i=0; i<k; i++) {
        set(G[i].a);
        for(j=0; j<n; j++) {
            power(f, F[j], V[i][j]);
            G[i].a *= f;
        }
        G[i].b = d[i];
    }
}

void ClassGroup(Vec<Pair<ImQF,ZZ> >& G)
// output:
//   G = vector of ("generator", "order") pair
//       in decreasing order of "order"
//       by Terr-Shanks method
// reference:
//   J. Buchmann and U. Vollmer
//    "Binary Quadratic Forms: An Algorithmic Approach"
//     Algorithm 9.3
{
    long i,j,k,l,m,p,q,c1,c2;
    double x;
    ZZ& h(ImQF::h);
    ZZ u,v,h1,h2;
    ImQF f,g,f1,g1;
    vec_long e;
    vec_vec_ZZ U;
    mat_ZZ U1;
    Vec<ImQF> F;
    Vec<Pair<ImQF, Vec<long> > > H1,H2,B;
    Vec<Pair<ImQF, Vec<ZZ> > > H;
    HashMap<ImQF, Vec<long> > B1(NumBits(ImQF::D)>>1);

    G.SetLength(0);
    ImQF::ClassNum();
    if(IsOne(h)) return;
    RightShift(u,h,2);
    sub(h1,h,u);
    add(h2,h,u);
    H.SetLength(1); set(H[0].a);
    B.SetLength(1); set(B[0].a);
    H1 = H2 = B;
    PrimeSeq ps;
    set(h);
    for(l=m=0, c1=c2=1;;) {
        p = ps.next();
        if(Kronecker(ImQF::D, p) <= 0) continue;
        SetPrime(f1,p);
        reduce(f1);
        f = f1;// baby element
        g = f1;// giant element
        B1.clear();
        for(i=0; i<B.length(); i++) {
            e = B[i].b;
            e.append(0);
            B1.install(B[i].a, e);
        }
        for(u=q=1;; q++, u+=q) {
            for(i=0; i<H.length(); i++) {
                mul(g1, g, H[i].a);
                if(B1.lookup(g1,e)) goto a;
            }
            for(i=0; i<B.length(); i++) {
                mul(g1, B[i].a, f);
                e = B[i].b;
                e.append(-q);
                B1.install(g1,e);
            }
            f *= f1;
            g *= f;
        }
    a:  ;
        u += e[l];
        if(IsOne(u)) continue;
        mul(v,h,u);
        if(v > h2) continue;
        F.append(f1);
        U.SetLength(l+1);
        U[l].SetLength(l+1);
        for(k=0; k<l; k++)
            add(U[l][k], H[i].b[k], e[k]);
        U[l][l] = u;
        h = v;
        if(h >= h1) break;
        for(j=0; j<H1.length(); j++) H1[j].b.append(0);
        for(j=0; j<H2.length(); j++) H2[j].b.append(0);
        conv(x,h);
        x = sqrt(x);
        if(l) {
            mul(v,c1,u);
            if(v < long(floor(x))) {
                conv(p,u);
                conv(c1,v);
                q = H1.length();
                inv(f1,f1);
                set(f);
                for(i=1, k=q; i<p; i++) {
                    f *= f1;
                    for(j=0; j<q; j++, k++) {
                        H1.SetLength(k+1);
                        mul(H1[k].a, H1[j].a, f);
                        H1[k].b = H1[j].b;
                        H1[k].b[l] += i;
                    }
                }
            }
            else {
                conv(p, U[m][m]);
                c2 *= p;
                q = H2.length();
                set(f);
                for(i=1, k=q; i<p; i++) {
                    f *= F[m];
                    for(j=0; j<q; j++, k++) {
                        H2.SetLength(k+1);
                        mul(H2[k].a, H2[j].a, f);
                        H2[k].b = H2[j].b;
                        H2[k].b[m] += i;
                    }
                }
                m = l;
            }
        }
        p = long(ceil(x/c1));
        q = long(ceil(x/c2));
        set(f);
        f1 = F[m];
        inv(f1,f1);
        k = B.length();
        for(i=0; i<k; i++) B[i].b.append(0);
        for(i=1; i<p; i++) {
            f *= f1;
            for(j=0; j<H1.length(); j++, k++) {
                B.SetLength(k+1);
                mul(B[k].a, H1[j].a, f);
                B[k].b = H1[j].b;
                B[k].b[m] += i;
            }
        }
        set(f);
        power(f1, F[m], p);
        k = H.length();
        for(i=0; i<k; i++) H[i].b.append(ZZ::zero());
        for(i=1, u=p; i<q; i++, u+=p) {
            f *= f1;
            for(j=0; j<H2.length(); j++, k++) {
                H.SetLength(k+1);
                mul(H[k].a, H2[j].a, f);
                conv(H[k].b, H2[j].b);
                H[k].b[m] += u;
            }
        }
        l++;
    }
    l = F.length();
    U1.SetDims(l,l);
    for(i=0; i<l; i++) VectorCopy(U1[i], U[i], l);
    HNF2SNF(G,F,U1);
}
