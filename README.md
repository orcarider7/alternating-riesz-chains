# The sharp charge-ratio threshold for alternating Riesz chains

Research note · 26 November 2026 · Maintained by Yingbo Sun (orcarider7)

[Literature comparison and priority search](PRIOR-ART.md) · [Citation information](CITATION.cff) · [Companion replication studies](https://github.com/orcarider7/paper-replication-lab)

**The first collapsing even period is exactly 78** for charges $6,-1$ with inverse-square interactions. Up to 76 sites per period the energy attains a finite minimum for every fixed density. Beginning at 78 sites it can decrease without bound. The uniform chain remains a strict local minimum at 78 and first acquires a negative second variation at 96. [Research manuscript](collapse-onset.md) · [PDF](first-collapse.pdf) · [Exact verification script](check_collapse_onset.py)

中文摘要：在一维周期链中，两类粒子的电荷交替为 $6$ 和 $-1$，相互作用能与距离平方成反比，同类相吸、异类相斥。固定周期长度时，每周期最多 76 个粒子（偶数）都有最低能量。从 78 个开始，可以通过压缩一组粒子让能量趋向负无穷。证明用整数运算覆盖了所有较小周期的排法。下载仓库后运行 **python check_collapse_onset.py** 即可复核。

The [four-particle proof](four-sites.md) classifies all ground states at the critical charge ratio. The [inverse-square Hessian and trial configurations](inverse-square.md) can be checked with [check_finite_cells.py](check_finite_cells.py).

For a one-dimensional periodic chain with alternating charges $m,-1$, write the pair interaction as $-q_iq_j|x_i-x_j|^{-p}$, with $m>0$ and $p>1$. Equal spacing uniquely minimizes the energy simultaneously at every density and every even period if and only if

$$
m+m^{-1}\le 2^{p+1}-2.
$$

Outside this range, the energy has infimum $-\infty$ at every density for every sufficiently large even period. The threshold sharpens the sufficient bound in Corollary 2.3 of Bétermin, Knüpfer, and Nolte. The proof adapts the exponential gap estimate in the Codex–Jiang follow-up by retaining the part of the variational problem that depends on the sum of the two charges. [Original paper](https://doi.org/10.1007/s10955-020-02603-2) · [Exponential gap estimate, Theorem 3](https://github.com/xiujiang4/alternating-particle-chains/blob/a48fe4c2623a29dd45edd227410a07baba400250/en/paper.tex)

## Model and theorem

Fix an even integer $N\ge 2$ and a mean gap $\ell>0$. Let $L=N\ell$, and take a strictly ordered chain with $x_{i+N}=x_i+L$ and $x_0=0$. Set $d_i=x_{i+1}-x_i$ and $\sigma_i=(-1)^i$. The charge $q_i$ is $m$ for even $i$ and $-1$ for odd $i$. The energy per particle for a radial kernel $f$ is

$$
E_f(X)=-\frac1N\sum_{i=0}^{N-1}
\sum_{\substack{j\in\mathbb Z\\j\ne i}}
q_iq_j f(|x_i-x_j|).
$$

Each pair contributes in both directions, as in the original paper. For $f(r)=r^{-p}$ with $p>1$ the inner sum is absolutely convergent. Write $E_p$ for this energy and $U_d$ for the equidistant chain of spacing $d$.

Define

$$
C_p(m)=\zeta(p)\big[2m(1-2^{-p})-(1+m^2)2^{-p}\big].
$$

**Theorem.** Every admissible chain satisfies

$$
E_p(X)\ge \frac{C_p(m)}{N}\sum_{i=0}^{N-1}d_i^{-p},
$$

with equality exactly for equidistant chains. Consequently, when $C_p(m)\ge 0$, the unique minimum is $E_p(U_\ell)=C_p(m)\ell^{-p}$. When $C_p(m)<0$, for every sufficiently large even $N$ and every $\ell>0$, the energy is unbounded below at fixed $N$ and $\ell$.

For a uniform chain, the average charge product of two sites $k$ steps apart is $-m$ for odd $k$ and $(1+m^2)/2$ for even $k$. Summing the two directions gives $E_p(U_d)=C_p(m)d^{-p}$. The factorization

$$
C_p(m)=m2^{-p}\zeta(p)\big[2^{p+1}-2-(m+m^{-1})\big]
$$

gives the stated threshold.

## Exponential gap estimate with unequal charges

For $\kappa>0$, abbreviate $E_{e^{-\kappa r}}$ by $E_\kappa$. Define the $L$-periodic kernel

$$
K_\kappa(t)=\sum_{n\in\mathbb Z}e^{-\kappa|t+nL|},
\qquad
Q_\kappa=\sum_{i,j=0}^{N-1}q_iq_jK_\kappa(x_i-x_j).
$$

The periodic kernel counts every image. Subtracting each site's own term yields $E_\kappa(X)=(m^2+1)/2-Q_\kappa/N$. On the circle of length $L$, the derivative jump of $K_\kappa$ at zero gives $(-\partial_t^2+\kappa^2)K_\kappa=2\kappa\delta_0$. The resulting variational identity is

$$
Q_\kappa=\sup_{u\in H^1_{\mathrm{per}}(0,L)}
\left\{
2\sum_i q_i u(x_i)
-\frac1{2\kappa}\int_0^L\big[(u')^2+\kappa^2 u^2\big]\,dt
\right\}.
$$

Indeed, the maximizing field $u_*(t)=\sum_j q_jK_\kappa(t-x_j)$ represents the linear functional $u\mapsto\sum_i q_i u(x_i)$ in the inner product defined by the integral divided by $2\kappa$. Completing the square gives the identity.

Split the integral over the gaps. On gap $i$, let $v=u(x_i)$, $w=u(x_{i+1})$, and $t_i=\tanh(\kappa d_i/2)$. For those endpoints, the minimum of the integral penalty on the gap is

$$
\frac14\left[t_i(v+w)^2+t_i^{-1}(v-w)^2\right].
$$

Center the gap at zero. The minimizing field solves $u''=\kappa^2u$ and is a linear combination of $\cosh(\kappa t)$ and $\sinh(\kappa t)$. Evaluating the integral by $(u u')'=(u')^2+\kappa^2u^2$ gives the displayed formula.

The linear part splits over gaps as $\sum_i(q_i v+q_{i+1}w)$, with $v,w$ the endpoints of gap $i$. Its gap contribution is

$$
q_i v+q_{i+1}w
=\frac{m-1}{2}(v+w)+\frac{(m+1)\sigma_i}{2}(v-w).
$$

Maximizing in $v+w$ and $v-w$ independently therefore gives

$$
Q_\kappa\le \sum_{i=0}^{N-1}\left[
\frac{(m-1)^2}{4t_i}+\frac{(m+1)^2t_i}{4}
\right].
\tag{1}
$$

The unique maximizing endpoint values for gap $i$ obey

$$
v+w=(m-1)/t_i,\qquad
v-w=(m+1)\sigma_i t_i.
\tag{2}
$$

When all gaps equal $d$, these endpoint values agree at shared sites because $\sigma_{i+1}=-\sigma_i$. Joining the minimizing fields on the gaps gives a periodic field attaining (1). Thus the bracket in (1) is the uniform $Q_\kappa/N$ at spacing $d_i$. With the self term accounted for, (1) reads

$$
E_\kappa(X)\ge \frac1N\sum_{i=0}^{N-1}E_\kappa(U_{d_i}).
\tag{3}
$$

For a nonuniform chain there is a pair of unequal neighboring lengths $d,e$ meeting at some site $i$. By (2), the difference between the two individually maximizing values for that site is

$$
\frac{t_d-t_e}{2}\left[(m+1)\sigma_i-\frac{m-1}{t_d t_e}\right],
\qquad t_s=\tanh(\kappa s/2).
$$

For each $\kappa>0$, the unequal lengths give $t_d\ne t_e$. As $\kappa\to\infty$, the bracket tends to $2$ if $\sigma_i=1$ and to $-2m$ if $\sigma_i=-1$. Thus the gap maxima are incompatible for all sufficiently large $\kappa$. Since the global variational maximum is attained, (3) is strict for those $\kappa$.

## From exponentials to Riesz kernels

For $r>0$, the inverse power is a positive mixture of exponentials.

$$
r^{-p}=\frac{1}{\Gamma(p)}\int_0^\infty
e^{-\kappa r}\kappa^{p-1}\,d\kappa.
$$

For every fixed admissible chain and $p>1$, the sum of the absolute values of the integrated pair contributions is finite. Indeed, $\lvert q_iq_j\rvert$ is bounded and the distances to image cells grow linearly with the image index. The same holds for every $U_{d_i}$. Fubini's theorem therefore justifies integrating (3) to obtain the theorem's gap inequality, with strict inequality for any unequal gaps.

When $C_p(m)\ge 0$, convexity of $d\mapsto d^{-p}$, together with $\sum_i d_i=N\ell$, gives $E_p(X)\ge C_p(m)\ell^{-p}$. The strict exponential estimate shows uniqueness even at $C_p(m)=0$.

## Collapse beyond the threshold

Assume $C_p(m)<0$. Place $N$ alternating sites at $x_i=i\delta$, for $0\le i<N$, and repeat the cell at period $L=N\ell$. For $0<\delta<L/[2(N-1)]$, the cell is admissible. The energy contribution from pairs within the same cell is $\delta^{-p}C_{N,p}(m)$, with

$$
C_{N,p}(m)=\sum_{k=1}^{N-1}
\left(1-\frac{k}{N}\right)k^{-p}
\begin{cases}
2m,&k\text{ odd},\\
-(1+m^2),&k\text{ even}.
\end{cases}
$$

The factor $1-k/N$ counts the fraction of sites that have a $k$-step neighbor ahead within the cell. When $k$ is even, the $N-k$ possible starting indices comprise equal numbers of the two species because $N-k$ is even. Since $p>1$, dominated convergence gives $C_{N,p}(m)\to C_p(m)$ through even $N$, so $C_{N,p}(m)<0$ for all sufficiently large even $N$.

For each fixed such $N$ and $\ell>0$, distances to sites in image cell $n\ne 0$ are at least $\lvert n\rvert L/2$. The total energy contribution from those images is bounded uniformly as $\delta\downarrow 0$, since $\sum_{n\ne 0}\lvert n\rvert^{-p}$ converges. The same-cell contribution tends to $-\infty$, proving the collapse statement.

## Convex envelope for positive exponential mixtures

Let $f(r)=\int_{(0,\infty)}e^{-\kappa r}\,d\mu(\kappa)$ for a positive measure $\mu$, with $f(r)<\infty$ for every $r>0$ and $\int_1^\infty f(r)\,dr<\infty$. Set $g(d)=E_f(U_d)$. Then the limit of the energy infimum at fixed mean spacing $\ell$ as the even period size tends to infinity is the lower convex envelope of $g$ at $\ell$, allowing the value $-\infty$.

Monotonicity of $f$ and its integrable tail give $\sum_{k\ge 1} f(kc)<\infty$ for every $c>0$. Thus absolute convergence permits integrating (3) against $\mu$. Every chain has $E_f(X)\ge N^{-1}\sum_i g(d_i)$, which is at least the lower convex envelope at the mean gap.

For the upper bound, fix $a,b>0$ and $0<\lambda<1$ such that $\ell=\lambda a+(1-\lambda)b$. Choose even $n_N$ with $n_N/N\to\lambda$. Set $n_N$ consecutive gaps to $a$ and all other gaps to $b_N=(N\ell-n_Na)/(N-n_N)$, which tends to $b$ and is positive for sufficiently large $N$. The mean gap is exactly $\ell$. For any fixed $k$, the interaction between sites $k$ gaps apart has at most $2(k-1)$ starting sites whose path crosses a gap-block boundary. The other starting sites give the uniform contributions at $a$ and $b_N$, with species counts balanced up to one site per block. For all sufficiently large $N$, every gap is bounded below by a common $c>0$. The absolute mean contribution of shell $k$, counting both directions, is therefore bounded by $2\max(1,m)^2 f(kc)$. Dominated convergence, and continuity of $f$ on $(0,\infty)$ from its integral representation, give the limiting energy $\lambda g(a)+(1-\lambda)g(b)$. Taking the infimum of all such chord values, including $g(\ell)$, gives the lower convex envelope. It suffices to use two points because minimizing a finite weighted average with prescribed mean takes place at a vertex of the weight polytope, where at most two weights are positive.

## Example and relation to prior work

For $p=2$, the threshold is $m+m^{-1}\le 6$, or $m^2-6m+1\le 0$. Equal spacing is the unique optimizer throughout

$$
3-2\sqrt 2\le m\le 3+2\sqrt 2.
$$

The upper endpoint is about $5.8284$. Corollary 2.3 of the original paper covers $m+m^{-1}<2^{p+1}/\zeta(p)$ for $p>p_1\approx 1.465$, which at $p=2$ has an upper endpoint of about $4.6483$. Thus $m=5$ is a concrete case covered by the sharpened bound. At $m=6$ the energy is unbounded below for all sufficiently large even cells. [Bétermin, Knüpfer, and Nolte, Corollary 2.3 and its proof](https://doi.org/10.1007/s10955-020-02603-2)

The exponential variational method reduces to the Codex–Jiang draft's neutral bound when $m=1$. The sum contribution in (1) vanishes in that case. The draft also relates its neutral result to earlier sharp-interface estimates. [Codex–Jiang draft, exponential kernels and relation to sharp-interface results](https://github.com/xiujiang4/alternating-particle-chains/blob/a48fe4c2623a29dd45edd227410a07baba400250/en/paper.tex)

### How large must a uniformly compressed block be?

At $p=2$ and $m=6$, the open-block coefficient $C_N=C_{N,2}(6)$ is negative exactly for even $N\ge162$. Thus, in this family, 162 is the first even cell size where reducing the internal gap to zero drives the energy to $-\infty$. Other choices of internal gaps are a separate finite-cell optimization problem.

To prove the cutoff, let $w_k=12$ for odd $k$ and $w_k=-37$ for even $k$. Set $B_N=\sum_{k<N,\,k\ge1} w_k/k$. Subtracting the two finite sums defining $C_{N+2}$ and $C_N$ gives, for even $N$,

$$
(N+2)(C_{N+2}-C_N)=2B_N/N-74/N^2+12/(N+1)^2.
$$

Here $B_4=-5/2$ and $B_{N+2}-B_N=-37/N+12/(N+1)<0$. The right-hand side is negative for $N\ge4$, and direct substitution gives the same sign at $N=2$. Exact rational summation in the accompanying script gives $0.002347<C_{160}<0.002349$ and $-0.001801<C_{162}<-0.001798$. Monotonicity proves the cutoff. At any fixed cell size the contribution from other image cells stays bounded during this compression, as shown above.

The accompanying [verification script](verify.py) evaluates numerical examples alongside the analytic proof.

OpenAI Codex derived and drafted these notes. Astra found the original 78-site rational witness, and Sol independently evaluated the signed pair energies as exact fractions. Terra audited the proof and the sign of the integer ceiling bound. Gemini reviewed the scaling argument conditional on the certificate. For the exact finite-period cutoff, Astra and Sol separately checked the 75-site integer certificate, Terra reviewed the reduction, and Gemini reviewed the analytic argument conditional on the integer checks. All verification scripts use the Python standard library. The optional flow generator uses SciPy. Formal peer review and literature priority remain outstanding.
