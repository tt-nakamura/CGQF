// uses NTL
//   http://www.shoup.net/ntl

#include "ImQF.h"
#include<NTL/HNF.h>
#include<NTL/LLL.h>

#define MCC_EXPMAX 20
#define MCC_EXTRA  10

void HNF2SNF(Vec<Pair<ImQF,ZZ> >&, const Vec<ImQF>&, const mat_ZZ&);

void ClassGroup2(Vec<Pair<ImQF,ZZ> >& G, double b)
// input:
//   b: controls size of factor base (default b=1)
// output:
//   G = vector of ("generator", "order") pair
//       in decreasing order of "order"
//       by McCurley method
// reference:
//   H. Cohen, "A Course in Algorithmic Algebraic Number Theory"
//     section 5.5
{
    const ZZ& D(ImQF::D);
    long i,j,k,l,p,N,K,M;
    double x,y,B;
    ZZ a,c;
    Vec<long> P,P0,v;
    mat_ZZ A,H;
    ImQF f;
    Vec<ImQF> F;
    Mat<ImQF> FE;

    conv(x,D);
    y = log(-x);
    x = exp(sqrt(y*log(y)/8));
    y = b*y*y;
    if(x > NTL_MAX_LONG || y > NTL_MAX_LONG)
        Error("too large discriminant");
    M = long(x>y ? x:y);

    PrimeSeq ps;
    conv(y,D);
    B = sqrt(-y)/(4*atan(1));
    y = sqrt(-y/3);
    for(i=0, x=1; (p = ps.next()) <= M;) {
        if(p==0) Error("too large discriminant");
        k = Kronecker(D,p);
        B /= 1 - k/double(p);
        if(k<0) continue;
        if(p==2) {
            k = trunc_long(D,4);
            if(k==0 || k==12) continue;
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
    K = N + MCC_EXTRA;

    FE.SetDims(M, MCC_EXPMAX);
    for(i=0; i<M; i++) {
        FE[i][0] = F[P0[i]];
        reduce(FE[i][0]);
        for(j=1; j<MCC_EXPMAX; j++)
            mul(FE[i][j], FE[i][j-1], FE[i][0]);
    }
    v.SetLength(N);
    A.SetDims(K,N);
    for(k=0;;) {
        while(k<K) {
            l = k%N;
            f = F[l];
            reduce(f);
            for(i=0; i<M; i++) {
                j = RandomBnd(MCC_EXPMAX);
                f *= FE[i][j];
                A[k][P0[i]] = j+1;
            }
            if(!f.a.SinglePrecision()) continue;
            conv(p, f.a);
            for(i=0; i<N; i++) v[i] = 0;
            for(i=0; i<N && p>=P[i]; i++)
                while(p%P[i]==0) { p/=P[i]; v[i]++; }
            if(p>1) continue;
            for(i=0; i<N; i++) {
                if(v[i]==0) continue;
                if(f.b%(P[i]<<1) <= P[i])
                     A[k][i] -= v[i];
                else A[k][i] += v[i];
            }
            A[k][l]++;
            if(k<N && IsZero(A[k][k])) {
                clear(A[k]);
                continue;
            }
            k++;
        }
        if(IsZero(c)) {
            if(image(a,A) < N) {
                clear(A);
                k=0;
                continue;
            }
            SqrRoot(c,a);
        }
        HNF(A,A,c);
        set(ImQF::h);
        for(i=0; i<N; i++) ImQF::h *= A[i][i];
        conv(x, ImQF::h);
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
