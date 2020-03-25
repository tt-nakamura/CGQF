// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/vector.h>
#include<NTL/pair.h>
using namespace NTL;

#define MAP_SIZE_MIN 1
#define MAP_SIZE_MAX 24

template<class T, class S>
struct HashMap {
    long n,l,i,j;
    Vec<Vec<Pair<T,S> > > elem;
    HashMap(long l1=0) : l(l1) {
        if(l < MAP_SIZE_MIN) l = MAP_SIZE_MIN;
        else if(l > MAP_SIZE_MAX) l = MAP_SIZE_MAX;
        n = (1<<l);
        elem.SetLength(n);
    }
    long hash_val(unsigned long v) {
        for(i=0; v;) {
            i ^= v&(n-1);
            v >>= l;
        }
        return i;
    }
    void clear() {
        elem.kill();
        elem.SetLength(n);
    }
    T* lookup(const T& t, S& s) {
        i = hash_val(t.hash_val());
        for(j=0; j < elem[i].length(); j++)
            if(elem[i][j].a == t) break;
        if(j == elem[i].length()) return 0;
        s = elem[i][j].b;
        return &(elem[i][j].a);
    }
    void install(const T& t, const S& s) {
        i = hash_val(t.hash_val());
        j = elem[i].length();
        elem[i].SetLength(j+1);
        elem[i][j].a = t;
        elem[i][j].b = s;
    }
};
