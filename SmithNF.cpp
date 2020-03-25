// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/mat_ZZ.h>
using namespace NTL;

void SmithNF(vec_ZZ& D, mat_ZZ& U, const mat_ZZ& A, const ZZ& m)
// input:
//   A = non-singular square matrix of integers
//   m = det(A)
// output:
//   D = diagonal elements of Smith normal form S of A
//       in decreasing order
//   U = square matrix of integers such that
//       A=VSU, det(U)=1, for some V, det(V)=1
// reference:
//   H. Cohen, "A Course in Computational Algebraic Number Theory"
//     Algorithm 2.4.14
{
    long i,j,k,n(A.NumRows());
    ZZ x,y,u,v,d,t,M(m);
    mat_ZZ B(A);
    D.SetLength(n);
    ident(U,n);
    for(k=n-1; k>0; k--) {
    a:  ;// eliminate column
        for(i=k-1; i>=0; i--) {
            if(IsZero(B[i][k])) continue;
            XGCD(d, x, y, B[i][k], B[k][k]);
            div(u, B[k][k], d);
            div(v, B[i][k], d);
            B[k][k] = d;
            clear(B[i][k]);
            for(j=k-1; j>=0; j--) {
                mul(t, x, B[i][j]);
                MulAddTo(t, y, B[k][j]);
                B[i][j] *= u;
                MulSubFrom(B[i][j], v, B[k][j]);
                B[i][j] %= M;
                rem(B[k][j], t, M);
            }
        }
    b:  ;// eliminate row
        for(j=k-1; j>=0; j--) {
            if(IsZero(B[k][j])) continue;
            XGCD(d, x, y, B[k][j], B[k][k]);
            div(u, B[k][k], d);
            div(v, B[k][j], d);
            B[k][k] = d;
            clear(B[k][j]);
            for(i=k-1; i>=0; i--) {
                mul(t, x, B[i][j]);
                MulAddTo(t, y, B[i][k]);
                B[i][j] *= u;
                MulSubFrom(B[i][j], v, B[i][k]);
                B[i][j] %= M;
                rem(B[i][k], t, M);
            }
            for(i=0; i<n; i++) {
                mul(t, v, U[j][i]);
                MulAddTo(t, u, U[k][i]);
                U[j][i] *= y;
                MulSubFrom(U[j][i], x, U[k][i]);
                U[j][i] %= M;
                rem(U[k][i], t, M);
            }
        }
        if(i>0) goto a;
        for(i=k-1; i>=0; i--) for(j=k-1; j>=0; j--) {
            if(divide(B[i][j], B[k][k])) continue;
            for(j=k-1; j>=0; j--) B[k][j] = B[i][j];
            goto b;
        }
        GCD(D[k], B[k][k], M);
        M /= D[k];
    }
    if(n) GCD(D[0], B[0][0], M);
}
