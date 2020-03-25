#include "ImQF.h"
#include<fstream>

enum { Shanks, McCurley };

void timer(const char *filename, long C, double d1, double d2, long n) {
    long i;
    double dx(pow(d2/d1, 1./n)), t;
    ZZ d;
    Vec<Pair<ImQF,ZZ> > G;
    std::ofstream ofs(filename);
    for(i=0; i<=n; i++) {
        conv(d, d1*pow(dx,i));
        while(!ImQF::init(d)) d--;
        t = GetTime();
        switch(C) {
            case Shanks:   ClassGroup(G); break;
            case McCurley: ClassGroup2(G); break;
        }
        t = GetTime() - t;
        std::cout << d << '\t' << t << '\n';
        ofs << d << '\t' << t << '\n';
    }
}

main() {
    timer("fig1a.txt", Shanks,   -5e3, -5e21, 100);
    timer("fig1b.txt", McCurley, -5e3, -1e30, 100);
}
