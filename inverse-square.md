# Finite-period stability and collapse for an alternating inverse-square chain

Consider ordered positions on a line with period length $L=Nh$, where $N$ is even and $h>0$. Write $x_{j+N}=x_j+L$, $q_{2j}=m$, and $q_{2j+1}=-1$, with $m>0$. The energy per particle is

$$
E(X)=-\frac1N\sum_{j=0}^{N-1}\sum_{k\in\mathbb Z\setminus\{0\}}
\frac{q_jq_{j+k}}{|x_{j+k}-x_j|^2}.
$$

This counts pairs in both directions and includes every image cell. The interaction attracts like charges and repels unlike charges. The uniform configuration is $x_j=jh$. We keep $N$ and $L$ fixed when varying the positions. The model is the $p=2$ case of the non-neutral inverse-power problem of Bétermin, Knüpfer and Nolte. [Original paper, Corollary 2.3](https://doi.org/10.1007/s10955-020-02603-2)

At $m=6$, period 96 is the first even period at which the uniform chain has a direction of negative second variation. A separate rational configuration already forces the energy to be unbounded below at period 78. At that same period the uniform chain is a strict local minimum modulo translation.

## The uniform-chain Hessian

**Theorem 1.** For even $N\ge4$ set $u_N=(2/N)(1-2/N)$. The Hessian of $E$ at the uniform chain is positive semidefinite if and only if

$$
m+m^{-1}\le H_N,\qquad H_N=6+8u_N+2u_N^2.
$$

Under strict inequality the uniform chain is a strict local minimizer modulo translation. If $m+m^{-1}>H_N$, arbitrarily small admissible displacements decrease its energy. At equality the Hessian has an additional zero direction and the quadratic calculation alone does not decide local minimality. At $N=2$ the relative displacement always has positive second variation.

**Proof.** Opposite neighbors at distances $kh$ and $-kh$ cancel each other's forces at each site, so the uniform chain is stationary. For a real $N$-periodic displacement $v$, its second variation is

$$
-\frac{6}{Nh^4}\sum_{j=0}^{N-1}\sum_{k\ne0}
\frac{q_jq_{j+k}}{k^4}(v_{j+k}-v_j)^2 .
$$

Termwise differentiation near the uniform chain is justified by uniform convergence of the derivative sums. The coefficients repeat every two sites. A complex Fourier block therefore has $v_j=a\exp(i\theta j)$ at even sites and $v_j=b\exp(i\theta j)$ at odd sites, with $\theta=2\pi l/N$. Adding $\pi$ to $\theta$ changes the sign convention for $b$ and leaves the block unchanged. Conjugation reverses the sign of $\theta$. It suffices to take $t=\theta/\pi$ in $[0,1/2]$, with $t=2l/N$ for $0\le l\le\lfloor N/4\rfloor$. Complexifying the real quadratic form preserves its sign.

Let $O=\sum_{k>0,\ k\ {\rm odd}}k^{-4}=\pi^4/96$, and put

$$
e=\frac1O\sum_{k>0,\ k\ {\rm even}}\frac{1-\cos k\theta}{k^4},
\qquad
d=\frac1O\sum_{k>0,\ k\ {\rm odd}}\frac{1-\cos k\theta}{k^4}.
$$

For even $k$, the terms contribute separately to $|a|^2$ and $|b|^2$ with charge factors $-m^2$ and $-1$. For odd $k$, summing the two orientations gives the real cross term proportional to $-2m(1-d)\operatorname{Re}(\bar a b)$. After removal of a common positive factor, the matrix is

$$
M(t)=
\begin{pmatrix}
m-m^2 e & -m(1-d)\\
-m(1-d) & m-e
\end{pmatrix}.
$$

The following cosine sum evaluates the entries for $0\le z\le2\pi$:

$$
F(z)=\sum_{k=1}^{\infty}\frac{1-\cos(kz)}{k^4}
=\frac{\pi^2z^2}{12}-\frac{\pi z^3}{12}+\frac{z^4}{48}.
$$

Indeed, integrating the cosine series of $(z-\pi)^2$ twice gives this formula, using $F(0)=F'(0)=0$ and $F''(z)=\pi^2/6-\pi z/2+z^2/4$. In particular $O=F(\pi)/2$. The even numerator is $F(2\theta)/16$ and the odd numerator is $F(\theta)-F(2\theta)/16$. With $u=t(1-t)$ this gives

$$
e=2u^2,\qquad d=2t^2(3-2t),\qquad d(2-d)=4u^2(3+4u).
$$

Consequently

$$
\det M(t)/m^2
=d(2-d)+e^2-(m+m^{-1})e
=e\big[6+8u+2u^2-(m+m^{-1})\big].
$$

One diagonal is strictly positive throughout the interval because $0\le e\le1/8$. Use $m-e$ when $m\ge1$ and $m-m^2e$ when $m\le1$. For $t>0$, positive semidefiniteness is therefore equivalent to nonnegativity of the bracket. Both $u$ and $6+8u+2u^2$ increase on $0\le t\le1/2$, so the smallest nonzero frequency $t=2/N$ decides the sign. At $t=0$ the matrix has eigenvalues $0$ and $2m$; the kernel is the common translation. This is the only block when $N=2$. Positive definiteness on the translation-free finite-dimensional displacement space proves strict local minimality. A negative direction gives a decrease by Taylor expansion, and a sufficiently small perturbation preserves ordering. $\square$

For $m=6$, $m+m^{-1}=37/6$, and exact arithmetic gives

$$
H_{94}-37/6=23183/29278086>0,\qquad
H_{96}-37/6=-7007/2654208<0 .
$$

Since $H_N$ decreases with even $N\ge4$, these comparisons prove the claimed first unstable period. For a particular decreasing direction at $N=96$, use $t=2/96$, $a=1$, $b=m(1-d)/(m-e)$ and the real part of the Fourier vector above. The quadratic form on $(1,b)$ is $\det M/(m-e)<0$.

## A finite cluster with an exact collapse coefficient

**Theorem 2.** Let $N\ge6$ be even, $0<\alpha<1$, $A=(m-1)/2$, and $B=(m+1)/2$. Define the ordered real positions

$$
y_j=-\cot\theta_j,\qquad \theta_j=\frac{\pi(j+\alpha)}{N},\qquad 0\le j<N .
$$

Their finite, unordered-pair energy is

$$
S_N=-\sum_{0\le i<j<N}\frac{q_iq_j}{(y_i-y_j)^2}
=-\frac N{16}\left[A^2(N-1)^2-B^2(N^2/2-1)\right].
$$

In particular $A^2(N-1)^2>B^2(N^2/2-1)$ suffices for $\inf E=-\infty$ at fixed $N$ and any $h>0$.

**Proof.** The cotangent difference formula gives

$$
(y_i-y_j)^{-2}=\sin^2\theta_i\sin^2\theta_j\,\csc^2(\pi(i-j)/N).
$$

Let $K_0=0$ and $K_k=\csc^2(\pi k/N)$ for $1\le k<N$. The real Fourier eigenvalues of this circulant kernel are

$$
D_l=\sum_{k=1}^{N-1}K_k\cos(2\pi lk/N)
=\frac{N^2-1}{3}-2l(N-l),\qquad 0\le l<N.
$$

For a direct verification, the second difference in $l$ of the cosine sum is $-4\sum_{k=1}^{N-1}\cos(2\pi lk/N)=4$ for $1\le l<N$, with periodic indices. This fixes the quadratic up to the additive constant, using periodicity. Its sum over all $l$ is zero because $K_0=0$, giving $(N^2-1)/3$.

Write $w_j=q_j\sin^2\theta_j$ and $\widehat w_l=N^{-1}\sum_j w_j\exp(-2\pi ilj/N)$. Since $q_j=A+B(-1)^j$ and $\sin^2\theta_j=(1-\cos(2\theta_j))/2$, the Fourier magnitudes at modes $0,\pm1$ are respectively $|A|/2, |A|/4, |A|/4$. Those at $N/2,N/2\pm1$ are $|B|/2, |B|/4, |B|/4$. These six modes are distinct for $N\ge6$ and the magnitudes are independent of $\alpha$. Orthogonality now gives

$$
S_N=-\frac N2\sum_l D_l|\widehat w_l|^2
=-\frac N{16}\left[A^2(2D_0+D_1)+B^2(2D_{N/2}+D_{N/2-1})\right].
$$

Substituting the eigenvalues proves the formula. For collapse, place $x_j=\delta(y_j-y_0)$ in a cell of length $L$ and take $\delta>0$ small enough that its span is less than $L/2$. The same-cell portion of $E$ equals $(2S_N/N)\delta^{-2}$. Pairs separated by $s\ne0$ cells have distance at least $|s|L/2$, hence their total contribution is bounded uniformly as $\delta\downarrow0$. A negative $S_N$ therefore gives the conclusion. $\square$

At $m=6$ the coefficient simplifies to $S_N=-N(N^2-100N+148)/128$. In the stated family its first negative value over even $N\ge6$ occurs at $100$, with $S_{98}=147/4$ and $S_{100}=-925/8$. Thus the doubly counted per-particle energy at period 100 has leading term $-(37/16)\delta^{-2}$. No conclusion about boundedness for all possible configurations at smaller periods follows from the sign of this trial coefficient.

## A rational certificate for a 78-site cell

Let $z_0<\cdots<z_{77}$ be the 78 integers in the [verification script](check_finite_cells.py), put $y_i=z_i/10^6$ and let $q_i=6,-1,6,-1,\ldots$. The unordered open-block energy is $S=\sum_{i<j}-q_iq_j/(y_j-y_i)^2$. Summing an integer upper bound on each signed term gives

$$
10^6 S \le -\sum_{i<j}\left\lfloor\frac{10^{18}q_iq_j}{(z_j-z_i)^2}\right\rfloor=-6\,134\,970.
$$

A numerical search suggested the coordinates. Once the integers are fixed, checking their order and summing 3,003 integer quotients suffices. Compress these sites into an interval of length less than $L/4$ by multiplying their translated positions by $\delta>0$. Their contribution is $(2S/78)\delta^{-2}$. The image sums stay bounded. For any greater even period $N$, place the additional alternating sites at distinct fixed positions in $(L/2,3L/4)$. All new distances to the compressed cluster and its images stay bounded away from zero, and the inverse-square image tails converge. The leading coefficient is then $2S/N<0$. This proves collapse for every even $N\ge78$.

At 78 the positive Hessian margin is exactly $159397/4626882$. The local minimum therefore coexists with a globally unbounded energy. The [exact cutoff proof](collapse-onset.md) establishes boundedness for every smaller even period, using a certificate valid for arbitrary ordered coordinates.

## Scope and reproducibility

The local stability test, cotangent trial and integer certificate can be checked independently of the all-exponent argument in the [all-period appendix](all-period-riesz.md). Codex agents derived and drafted the formulas. An Astra review produced the 78-site witness, which a separate Sol review checked by exact rational summation. Literature priority and formal peer review are outstanding.
