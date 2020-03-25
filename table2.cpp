#include "ReQF.h"

main() {
    long i,j,k,k1(0),k2(373),n(30);
    double R;
    ZZ d,d0;
    Vec<Pair<ReQF,ZZ> > G;
    power(d0,10,n);
    for(k=k1; k<k2; k++) {
        add(d,d0,k);
        if(!ReQF::init(d)) continue;
        ClassGroup2(G);
        conv(R, ReQF::regulator());
        std::cout << "$10^{" << n << "}+" << k << "$ & ";
        std::cout << ReQF::ClassNum() << " & ";
        std::cout << R << " & $";
        std::cout << ReQF::E << "$ & $";
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
