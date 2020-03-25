#include "ImQF.h"

main() {
    long i,j,k,k1(0),k2(276),n(30);
    ZZ d,d0;
    Vec<Pair<ImQF,ZZ> > G;
    power(d0,10,n);
    negate(d0,d0);
    for(k=k1; k<k2; k++) {
        sub(d,d0,k);
        if(!ImQF::init(d)) continue;
        ClassGroup2(G);
        std::cout << "$-10^{" << n << "}-" << k << "$ & $";
        std::cout << ImQF::ClassNum() << "$ & $";
        for(i=j=0; i<G.length(); i++) {
            if(i && G[i].b == G[i-1].b) j++;
            else {
                if(j) std::cout << "^{" << j+1 << '}';
                if(i) std::cout << "\\cdot ";
                std::cout << G[i].b;
                j=0;
            }
        }
        if(j) std::cout << "^{" << j+1 << '}';
        std::cout << "$\\\\\n";
    }
}
