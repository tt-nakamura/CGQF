#include "ReQF.h"
#include<fstream>

enum { Shanks, Buchmann };

void timer(const char *filename, long C, double d1, double d2, long n) {
    long i;
    double dx(pow(d2/d1, 1./n)), t;
    ZZ d;
    Vec<Pair<ReQF,ZZ> > G;
    std::ofstream ofs(filename);
    for(i=0; i<=n; i++) {
        conv(d, d1*pow(dx,i));
        while(!ReQF::init(d)) d++;
        t = GetTime();
        switch(C) {
            case Shanks:   ClassGroup(G); break;
            case Buchmann: ClassGroup2(G); break;
        }
        t = GetTime() - t;
        std::cout << d << '\t' << t << '\n';
        ofs << d << '\t' << t << '\n';
    }
}

main() {
    timer("fig2a.txt", Shanks,   5e3, 1e15, 50);
    timer("fig2b.txt", Buchmann, 5e3, 1e30, 100);
}
