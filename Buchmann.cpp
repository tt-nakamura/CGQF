// uses NTL
//   http://www.shoup.net/ntl

#include "ReQF.h"
#include<NTL/HNF.h>
#include<NTL/LLL.h>
#include<NTL/mat_GF2.h>

#define BUC_EXPMAX 20
#define BUC_EXTRA  10

void RGCD(RR&, const vec_RR&, double=0.2);
void HNF2SNF(Vec<Pair<ReQF,ZZ> >&, const Vec<ReQF>&, const mat_ZZ&);

void ClassGroup2(Vec<Pair<ReQF,ZZ> >& G, double b)
// input:
//   b: controls size of factor base (default b=1)
// output:
//   G = vector of ("generator", "order") pair
//       in decreasing order of "order"
//       by Buchmann method
// reference:
//   H. Cohen, "A Course in Computational Algebraic Number Theory"
//     section 5.9
{
    const ZZ& D(ReQF::D);
    long i,j,k,l,p,N,K,M;
    double x,y,B;
    RR z;
    ZZ a,c;
    Vec<long> P,P0,e,v;
    vec_RR u;
    vec_GF2 s;
    mat_GF2 S;
    mat_ZZ A,H,V;
    mat_RR U;
    ReQF f,f0;
    ReQFD g;
    Vec<ReQF> F;
    Mat<ReQFD> FE;

    y = log(D);
    x = exp(sqrt(y*log(y)/8));
    y = b*y*y;
    if(x > NTL_MAX_LONG || y > NTL_MAX_LONG)
        Error("too large discriminant");
    M = long(x>y ? x:y);

    PrimeSeq ps;
    conv(y, ReQF::SD);
    B = y/2;
    for(i=0, x=1; (p = ps.next()) <= M;) {
        if(p==0) Error("too large discriminant");
        k = Kronecker(D,p);
        B /= 1 - k/double(p);
        if(k<0) continue;
        if(p==2) {
            k = trunc_long(D,4);
            if(k==0 || k==4) continue;
        }
        else if(divide(a,D,p) && divide(a,p)) continue;
        P.append(p);
        F.SetLength(i+1);
        SetPrime(F[i++], p);
        if(x>y || divide(D,p)) continue;
        P0.append(i-1);
        x *= p;
    }
    B *= sqrt(2);
    N = P.length();
    M = P0.length();
    K = N + BUC_EXTRA;

    FE.SetDims(M, BUC_EXPMAX);
    for(i=0; i<M; i++) {
        conv(FE[i][0], F[P0[i]]);
        reduce(FE[i][0]);
        for(j=1; j<BUC_EXPMAX; j++)
            MulRed(FE[i][j], FE[i][j-1], FE[i][0]);
    }
    e.SetLength(M);
    v.SetLength(N);
    u.SetLength(K);
    s.SetLength(K);
    A.SetDims(K,N);
    for(k=0;;) {
        for(; k<K; k++) {
            l = k%N;
        a:  ;
            f = F[l];
            reduce(f);
            for(i=0; i<M; i++) {
                e[i] = RandomBnd(BUC_EXPMAX);
                MulRed(f, f, FE[i][e[i]]);
                A[k][P0[i]] = e[i] + 1;
            }
            f0=f;
            for(j=0;; j++) {
                if(f.a.SinglePrecision()) {
                    conv(p, f.a);
                    for(i=0; i<N; i++) v[i] = 0;
                    for(i=0; i<N && p>=P[i]; i++)
                        while(p%P[i]==0) { p/=P[i]; v[i]++; }
                    if(p==1) break;
                }
                rho(f);
                if(((j&1)==0 && f.a==f0.a) ||
                   ((j&1)!=0 && f.b==f0.b)) goto a;
            }
            for(i=0; i<N; i++) {
                if(v[i]==0) continue;
                if(f.b%(P[i]<<1) <= P[i])
                     A[k][i] -= v[i];
                else A[k][i] += v[i];
            }
            A[k][l]++;
            if(k<N && IsZero(A[k][k])) {
                clear(A[k]);
                goto a;
            }
            if(!IsZero(ReQF::R)) continue;
            conv(g, F[l]);
            reduce(g);
            for(i=0; i<M; i++) MulRed(g, g, FE[i][e[i]]);
            while(j--) rho(g);
            s[k] = (sign(g.d)<0);
            abs(u[k], g.d);
            log(u[k], u[k]);
            u[k] /= 2;
        }
        if(IsZero(c)) {
            k = image(a,A,V);
            if(IsZero(ReQF::R)) {
                V.SetDims(K-k, K);
                LLL_XD(V);
                conv(U,V); mul(u,U,u);
                conv(S,V); mul(s,S,s);
                RGCD(ReQF::R, u);
                ReQF::E = (IsZero(s) ? 1:-1);
            }
            if(k<N) {
                clear(A);
                k=0;
                continue;
            }
            SqrRoot(c,a);
        }
        HNF(A,A,c);
        set(ReQF::h);
        for(i=0; i<N; i++) ReQF::h *= A[i][i];
        conv(x, ReQF::h);
        conv(y, ReQF::R); x*=y;
        if(x<B) break;
        A.SetDims(K,N);
        for(i=k=N; i<K; i++) clear(A[i]);
    }
    for(i=k=0; i<N; i++)
        if(!IsOne(A[i][i])) k++;
    H.SetDims(k,k);
    for(i=k=0; i<N; i++) {
        if(IsOne(A[i][i])) continue;
        for(j=l=0; j<=i; j++)
            if(!IsOne(A[j][j])) H[k][l++] = A[i][j];
        F[k++] = F[i];
    }
    F.SetLength(H.NumRows());
    HNF2SNF(G,F,H);
}
