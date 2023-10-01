#include <bits/stdc++.h>
using namespace std;

#define rep(i, a, b) for (int i = a; i < (int)b; i++)
#define mem(a, b) memset(a, b, sizeof(a))
typedef long long ll;

const int maxn = (1 << 17) + 9;
const long double pi = acos(-1.0L);
using C = complex<long double>;
void fft(C a[], int n, int ty) {
    static C w[maxn * 4];
    w[0].real(1);
    for (int i = 0, j = 0; i < n; i++) {
        if (i > j) swap(a[i], a[j]);
        for (int l = n / 2; (j ^= l) < l; l /= 2)
            ;
    }
    for (int i = 1; i < n; i *= 2) {
        C wn(cos(pi / i), ty * sin(pi / i));
        for (int j = (i - 2) / 2; ~j; j--)
            w[j * 2 + 1] = (w[j * 2] = w[j]) * wn;
        for (int j = 0; j < n; j += i * 2)
            for (int k = j; k < j + i; k++) {
                C x = a[k], y = a[k + i] * w[k - j];
                a[k] = x + y, a[k + i] = x - y;
            }
    }
    if (ty != 1)
        for (int i = 0; i < n; i++) a[i] /= n;
}
template <class ty>
void operator*=(vector<ty>& a, vector<ty>& b) {
    static C c[maxn * 4], d[maxn * 4];
    int n = a.size(), m = b.size(), len = 1 << int(ceil(log2(n + m)));
    rep(i, 0, n) c[i] = a[i];
    fill(c + n, c + len, C());
    rep(i, 0, m) d[i] = b[i];
    fill(d + m, d + len, C());
    fft(c, len, 1), fft(d, len, 1);
    rep(i, 0, len) c[i] *= d[i];
    fft(c, len, -1);
    a.resize(len = n + m - 1);
    rep(i, 0, len) a[i] = c[i].real() + 0.5;
}
