# Exact ground states of a four-site alternating Riesz chain

Research note · 26 November 2026

A one-dimensional periodic cell contains two charges of each species in alternating order. For the interaction $-q_iq_j r^{-p}$ with $p>1$, the four-site problem admits an exact threshold. Let $M$ be the larger charge magnitude divided by the smaller. Equal spacing is the unique ground state when $M<2^{p+1}$. At equality, the two larger charges may have any positive separation, provided each smaller charge bisects its arc between them. Above the threshold the energy is unbounded below.

The starting point is the alternating-chain model of Bétermin, Knüpfer and Nolte. Their Corollary 2.3 gives a sufficient condition for crystallization for power-law interactions across all even periods. This note treats a prescribed period of four sites and gives the minimizers at its sharp boundary. The Gaussian counterexample in the Codex–Jiang manuscript concerns a different potential. [Bétermin, Knüpfer and Nolte (2020)](https://doi.org/10.1007/s10955-020-02603-2) · [Codex–Jiang manuscript](https://github.com/xiujiang4/alternating-particle-chains) · [All-period Riesz calculation](all-period-riesz.md)

## Energy and normalization

Fix the circumference $L>0$. After translation, the two larger charges $+M$ occupy $0,x$ and the smaller charges $-1$ occupy $y,z$, with $0<y<x<z<L$. This convention also covers charges $m,-1$ for $0<m<1$. Negating all charges, dividing them by $m$, and changing the starting site gives $M=1/m$. The energy in the original units equals $m^2$ times the normalized energy.

Use the energy per particle with each pair counted in both directions. Write $X=(0,y,x,z)$, $q=(M,-1,M,-1)$, and

$$
E=-\frac14\sum_{i,j=0}^3
\sum_{\substack{n\in\mathbb Z\\(j,n)\ne(i,0)}}
\frac{q_iq_j}{|X_j+nL-X_i|^p},
\qquad
S(t)=\sum_{n\in\mathbb Z}|t+nL|^{-p},\quad t\notin L\mathbb Z.
$$

Isolating the four self-image sums gives

$$
E=\tfrac12 F(x,y,z)-(M^2+1)L^{-p}\zeta(p),
$$

where $F$ counts each distinct-site pair once and has the following expression.

$$
F(x,y,z)=M\,[S(y)+S(x-y)+S(z-x)+S(L-z)]
       -M^2 S(x)-S(z-y).
$$

Indeed, each self-image sum is $2L^{-p}\zeta(p)$ and $\sum_i q_i^2=2(M^2+1)$. The series for $S$ and its first two derivatives converge locally uniformly away from $L\mathbb Z$. The kernel is even, $L$-periodic, and strictly convex on $(0,L)$.

## The two smaller charges bisect the arcs

Fix $x$. The terms $-M^2 S(x)$ and the self-image contribution are constant. We show that $F$ is strictly convex jointly in $(y,z)$ on $(0,x)\times(x,L)$.

Let $(v,w)\ne(0,0)$ be a displacement of the two smaller charges and put $s=p+2$. The pair between $y$ and $z+nL$ contributes the negative of $p(p+1)(v-w)^2/|z+nL-y|^s$ to the second directional derivative. For each $n$ choose an image $h_n$ of a larger charge by the following rule for $k\in\mathbb Z$.

$$
h_{2k}=x+kL,\qquad h_{2k+1}=(k+1)L.
$$

In each case $h_n$ lies strictly between $y$ and $z+nL$. Thus $a=|y-h_n|$ and $b=|z+nL-h_n|$ are positive and $a+b=|z+nL-y|$. Weighted Cauchy–Schwarz and $s>1$ give

$$
\frac{(v-w)^2}{(a+b)^s}
\ \le\ \frac{a^s+b^s}{(a+b)^s}
          \left(\frac{v^2}{a^s}+\frac{w^2}{b^s}\right)
\ <\ \frac{v^2}{a^s}+\frac{w^2}{b^s}.
$$

The sequence $h_n$ visits every image of the two larger charges exactly once. So does $h_n-nL$, whose even and odd terms are $x-kL$ and $-kL$. On summing, the right side accounts exactly once for all four smaller–larger pair sums in the directional derivative, with their common coefficient $M$ omitted. All sums converge, the inequality is strict, and $M\ge1$. Hence the directional derivative of second order is positive for every nonzero $(v,w)$.

Reflection in the midpoint of the arc from $0$ to $x$ interchanges the larger charges and maps $(y,z)$ to $(x-y,L+x-z)$. These two configurations have the same energy. Their average is

$$
(y_*,z_*)=\left(\frac{x}{2},\frac{L+x}{2}\right).
$$

Strict convexity gives the unique minimizer for this fixed $x$, including when $M=1$.

## The remaining separation

With the smaller charges at the arc midpoints, their separation is $L/2$. Splitting the image sum into even and odd integers gives $S(x/2)+S((L-x)/2)=2^p S(x)$. Consequently

$$
\min_{y,z}F(x,y,z)=M(2^{p+1}-M)S(x)-S(L/2).
$$

The function $S(x)$ has its unique minimum at $L/2$ by strict convexity and reflection symmetry, and it diverges as $x\downarrow0$. The sign of its coefficient now determines the ground states.

| Charge ratio | Ground states at fixed $L$ |
| --- | --- |
| $1\le M<2^{p+1}$ | $x=L/2$ and all four gaps equal $L/4$. |
| $M=2^{p+1}$ | Every $0<x<L$ minimizes with $(y,z)=(x/2,(L+x)/2)$. The gaps are $x/2,x/2,(L-x)/2,(L-x)/2$ in cyclic order. |
| $M>2^{p+1}$ | Sending $x$ to zero along these midpoint configurations gives $E\to-\infty$. |

Ground states are understood up to translation and the choice of starting site. For a normalization check, take $p=2$, $L=4$, and $M=8$. Both the uniform sites $(0,1,2,3)$ and the nonuniform sites $(0,1/2,1,5/2)$ have energy per particle $-17\pi^2/24$. At $M=6$ the four-site chain crystallizes even though sufficiently large cells at the same charge ratio can collapse. [Large-cell construction](all-period-riesz.md#collapse-beyond-the-threshold)

OpenAI Codex derived and drafted this note. The proof is presented for independent checking. Priority and publication status require separate verification.
