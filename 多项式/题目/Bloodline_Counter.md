# Bloodline Counter

出处：2023牛客暑期多校训练营8，B题
https://ac.nowcoder.com/acm/contest/57362/B

题意：给定 $n$ , $k$ ，求 $n$ 个点的编号竞赛图满足最大环长度恰为k的方案数。

令 $G(x)$ 代表x个点的竞赛图的数量所构成的指数生成函数，易得 $[x^n] = \frac{2^{(n*(n-1)/2)}}{n!}$ 。

令 $F(x)$ 代表x个点的单个强连通分量的数量所构成的生成函数，则有 $F(x)+F(x)^2+F(x)^3+... = G(x)$ 。所以 $F(x) = \frac{G(x)}{1+G(x)}$ 。

现在要求最大环长度恰为k的方案数，那么就是长度小于等于k的方案数减去长度小于等于k-1的方案数。

令 $H_{k}(x)$ 代表x个点的的最大环长度小于等于k的数量所构成的生成函数，那么 $H_k(x) = F(x)+F(x)^2+F(x)^3+...$ ，即 $H_k(x) = \frac{F(x)}{1-F(x)}$ ，其中 $F(x)$ 最高项的次数为k。

同理， $H_{k-1}(x) = F(x)+F(x)^2+F(x)^3+...$ ， 即 $H_{k-1}(x) = \frac{F(x)}{1-F(x)}$ ，其中 $F(x)$ 最高项的次数为 $k-1$ 。

代码实现：

```cpp
int main() {
    pre(); // 快速读写
    init(); // 预处理阶乘和阶乘逆元
    ll n, k;
    cin >> n >> k;
    Poly g(k + 1); 
    for (int i = 0; i <= k; i++) //注意i*(i-1)/2可能会爆int，需要开long long
        g[i] = qmi(2, i * (i - 1) / 2, mod) * invfac[i] % mod; //指数生成函数
    Poly f = Inv(g, k + 1); 
    f.resize(n + 1); // 注意求逆前要把f开到deg大小
    ll ans = 0;
    Poly h = Inv(f, n + 1);
    ans += h[n];
    f[k] = 0; // 求k-1的时候要把k的系数置为0
    h = Inv(f, n + 1);
    ans = (ans + mod - h[n]) % mod;
    ans = ans * fac[n] % mod; // 从指数生成函数还原，乘上n!即为答案
    cout << ans << endl;
}
```
