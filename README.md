# When do alternating Riesz chains collapse?

<a name="the-sharp-charge-ratio-threshold-for-alternating-riesz-chains"></a>

Research manuscript · Yingbo Sun (orcarider7), maintainer

**With alternating charges $+6,-1$ and inverse-square interactions, 78 is the smallest even number of sites per period that can collapse.** At any fixed positive cell length, the energy attains a finite minimum through 76 sites. At 78 sites a group of particles can move arbitrarily close together and drive the energy to $-\infty$. At that same cell size, the equally spaced configuration is still a strict local minimum up to translation.

The cutoff of 78 is one instance of a general rule. Increasing the charge ratio makes collapse possible in smaller cells. In the inverse-square model every even period from four onward is the first to collapse for some interval of charge ratios. The finite-period critical ratios decrease strictly from $8$ to $3+2\sqrt2$. Near charge ratio six, the first collapsing period decreases from 80 at $m=5.995$ to 78 at $m=6$ and 76 at $m=6.002$.

[Reproduce the certificates](#reproducing-the-results) · [Literature and priority](PRIOR-ART.md) · [Citation information](CITATION.cff)

中文摘要：在固定周期长度和反平方相互作用下，电荷交替为 $+6,-1$ 时，78 是能量可降到负无穷的最小偶数周期。每周期最多 76 个粒子时存在最低能量，78 个粒子时，扣除整体平移后，等距构型仍是严格局部极小值。只改变电荷大小之比，就可以让从 4 开始的任意偶数成为第一个可塌缩的周期。运行 `python check_collapse_onset.py` 和 `python check_charge_windows.py` 可以复核相关证书。

## Model and theorem

A cell contains an even number $N\ge2$ of particles in strict order, with $x_{i+N}=x_i+L$ for cell length $L>0$. The charges alternate between $q_{2j}=m$ and $q_{2j+1}=-1$. The pair energy is $-q_iq_j/|x_i-x_j|^p$. Like charges attract and opposite charges repel in this convention. For $p>1$, the energy per particle is the absolutely convergent sum

$$
E_p(X)=-\frac1N\sum_{i=0}^{N-1}
\sum_{\substack{j\in\mathbb Z\\j\ne i}}
\frac{q_iq_j}{|x_i-x_j|^p}.
$$

It counts each pair in both directions and includes all image cells. The mean gap is $\ell=L/N$. With $N,L$ fixed, collapse means $\inf_X E_p(X)=-\infty$ over the admissible ordered configurations. The transformation $x_i\mapsto a x_i$ for $a>0$ scales the energy by $a^{-p}$, so collapse at a given period is independent of the positive cell length. We take $m\ge1$ for the finite-period classification. Exchanging species, reversing the charge signs and rescaling the magnitudes maps a ratio $0<m<1$ to $1/m$ and changes the energy by a positive factor. [Original periodic model](https://doi.org/10.1007/s10955-020-02603-2) · [Four-site normalization](four-sites.md#energy-and-normalization)

### A finite block decides collapse

For $p=2$, write $S_n(x;m)$ for the open energy of $n$ ordered sites, with pairs counted once. In a heavy-ended block the first and last sites both have charge $m$. A periodic cell of even size $N$ collapses exactly when some heavy-ended open block of $N-1$ sites has negative $S_{N-1}$.

To see the lower bound, cut the periodic cell at its largest gap, which is at least $L/N$. The terms involving other image cells are then uniformly bounded. Orient the remaining cell from a heavy site to a light site. From the terminal light, each earlier light can be paired with a closer heavy, so removing the terminal light lowers the open energy. Nonnegative energy for every heavy-ended block of $N-1$ sites therefore bounds the periodic energy below. For the converse, compress a negative open block by a factor $\delta\downarrow0$ and keep the terminal light at a fixed positive distance. The internal block contributes $2S_{N-1}/(N\delta^2)$ to $E_2$ and all other terms stay bounded. [Detailed reduction](collapse-onset.md#1-reduction-to-an-open-chain) · [Variable-charge version](collapse-ladder.md#periodic-consequence)

## The finite-period thresholds

For $k\ge1$ let $\mu_k$ be the infimum of charge ratios $m\ge1$ at which some heavy-ended open block of $2k+1$ sites has negative energy, and set $\mu_0=+\infty$. Then

$$
8=\mu_1>\mu_2>\cdots\downarrow 3+2\sqrt2.
$$

Every such open block has positive energy for $1\le m<\mu_k$. At $m=\mu_k$, all such blocks have nonnegative energy and a finite block attains zero. Above $\mu_k$, a negative block exists. Together with the reduction, this gives the classification

$$
\text{First collapsing period }N=2k+2
\quad\Longleftrightarrow\quad
\mu_k<m\le\mu_{k-1}\qquad(k\ge1).
$$

Period two is always bounded below. At the exact threshold $\mu_k$, period $2k+2$ is still bounded below. Attainment of the periodic infimum at a critical ratio requires a separate argument. For $1\le m\le3+2\sqrt2$, no finite period collapses. [Threshold theorem and proof](collapse-ladder.md)

For a fixed block shape, the open energy is the quadratic $mC-m^2A-B$, where $C,A,B$ are the inverse-square distance sums for mixed, heavy-heavy and light-light pairs respectively. Dividing each same-charge interval into blocks of odd index length gives a strict finite-block lower bound at $m=3+2\sqrt2$. At each critical ratio, a finite zero-energy block can be extended by a distant light and heavy whose net energy increment is negative. This lowers the next threshold strictly. To obtain a finite zero at that next threshold, normalize the smallest gap of a minimizing sequence to one. Any gap that tends to infinity would split the block into shorter clusters with a strictly positive total energy, contradicting the zero limit. Finally, uniformly spaced open blocks have negative limiting energy per site whenever $m>3+2\sqrt2$, giving the limit of the decreasing sequence. [Dyadic bound](collapse-ladder.md#charge-dependence-and-the-limiting-threshold) · [Strict decrease and attainment of open-block zeros](collapse-ladder.md#strict-decrease-and-attainment)

### Four and six sites

- **Four sites.** For three open sites with successive gaps $a,b$, the energy is $m(a^{-2}+b^{-2})-m^2(a+b)^{-2}$. Convexity gives $(a+b)^2(a^{-2}+b^{-2})\ge8$, with equality at $a=b$. Hence $\mu_1=8$. More generally, the exact four-site threshold is $2^{p+1}$ for $p>1$. Below it the equally spaced chain is the unique ground state up to translation. At the threshold the two heavy sites may have any positive separation, and each light site bisects its arc between them. [Four-site ground states](four-sites.md)
- **Six sites.** The critical ratio $\mu_2=16(1+3u)/(9(1-u)^3)$ has certified enclosure $7.44624222875<\mu_2<7.44624223510$. Its parameter $u$ is the unique root in $(1/4,13/50)$ of

  $$
  256u^2(3+u)(1+3u)-128(1+3u)(1-u)^3+9(1-u)^6=0.
  $$

  The proof uses a stationary symmetric five-site block at the threshold and a convexity certificate valid for every five-site shape. Below $\mu_2$ the periodic minimum is attained. At $\mu_2$ the periodic infimum is finite. [Algebraic threshold and global certificate](six-site-threshold.md)

## Why the cutoff is exactly 78 at charge ratio six

Two exact computations settle the cases on either side of the cutoff.

1. **Through period 76.** An integer flow proves that the absolute same-charge contribution is at most $\rho$ times the mixed-charge contribution for every heavy-ended 75-site block, with $\rho=228876486320393532885019483/228881835937500000000000000<1$. For increasing reference coordinates $B$ and arbitrary ordered coordinates $x$, convexity of the cube bounds $(B_v-B_s)^3/(x_v-x_s)^2$ by the sum of the corresponding expressions along a path from $s$ to $v$. The flow routes same-charge terms along mixed edges. Its exact checked capacities give the strict margin. Padding with sites that recede to infinity extends the bound to shorter blocks. The margin also controls all interior adjacent gaps. The terminal-light contribution controls its adjacent gap, so the energy tends to $+\infty$ whenever any cyclic gap tends to zero. This gives attainment by compactness of the closed gap simplex. [Universal lower certificate and attainment](collapse-onset.md#2-a-certificate-valid-for-every-choice-of-the-75-positions) · [Integer flows](flow-75.json)
2. **From period 78 onward.** A 77-site block has evenly indexed sites at $g((j-38)/38)$, where $g(t)=2t+t^5$, and odd sites at the midpoints of their neighbors. Signed integer arithmetic gives $10^6 S_{77}\le -36{,}340{,}297$. Compressing this block drives the energy to $-\infty$ in every even period $N\ge78$, with the remaining sites held at distinct positive distances. [Rational upper witness](collapse-onset.md#3-a-negative-77-site-block)

The same methods give strict enclosures of width $10^{-7}$ for the four critical ratios for periods 74, 76, 78 and 80. A flow checks positivity for every open shape at the lower endpoint of each enclosure. A rational shape has negative energy at the upper endpoint. These enclosures and the strict ordering give the following certified ranges, including the endpoints shown.

- From $5.9921304$ to $5.9961805$, the first collapsing period is 80.
- From $5.9961806$ to $6.0004364$, the first collapsing period is 78.
- From $6.0004365$ to $6.0049142$, the first collapsing period is 76.

The transition between each consecutive range lies strictly within the intervening gap. [Exact endpoint data and verification](charge-windows/README.md)

## When equal spacing loses local stability

At $p=2$, for even $N\ge4$ set $u_N=(2/N)(1-2/N)$. The Hessian at equal spacing is positive semidefinite exactly when

$$
m+m^{-1}\le H_N,\qquad H_N=6+8u_N+2u_N^2.
$$

Strict inequality gives a strict local minimum modulo translation. Above the bound there is a negative second variation. At equality an additional zero direction leaves the local question open at quadratic order. At $m=6$, the inequality is strict at 78 and the first negative second variation occurs at 96. Thus a chain at period 78 can resist every sufficiently small displacement of the equally spaced state even though a compressed block gives arbitrarily low energy at the same period and density. [Fourier derivation and exact checks](inverse-square.md#the-uniform-chain-hessian)

## Collapse beyond the threshold

For general $p>1$, equal spacing uniquely minimizes the energy up to translation at every positive density and every even period exactly when

$$
m+m^{-1}\le 2^{p+1}-2.
$$

If $d_i=x_{i+1}-x_i$, the exponential-kernel argument gives

$$
E_p(X)\ge\frac{C_p(m)}{N}\sum_{i=0}^{N-1}d_i^{-p},\qquad
C_p(m)=m2^{-p}\zeta(p)\big[2^{p+1}-2-(m+m^{-1})\big].
$$

The inequality is strict for unequal gaps, even when the coefficient is zero. A uniformly spaced chain with gap $d$ has energy $C_p(m)d^{-p}$. When $C_p(m)\ge0$, convexity in the gaps gives the minimum value $C_p(m)\ell^{-p}$. When $C_p(m)<0$, sufficiently long uniformly spaced open blocks have negative energy. Compressing such a block proves collapse for all sufficiently large even periods. At $p=2$ this gives the limiting charge ratio $3+2\sqrt2$ for $m\ge1$. At $m=6$, the first even period that collapses by uniform compression of the entire open cell is 162. The 77-site certificate illustrates how varying the gaps can substantially lower the first collapsing period. <a name="how-large-must-a-uniformly-compressed-block-be"></a> [All-exponent proof](all-period-riesz.md) · [Uniform compression calculation](all-period-riesz.md#how-large-must-a-uniformly-compressed-block-be)

<a name="exponential-gap-estimate-with-unequal-charges"></a>
<a name="from-exponentials-to-riesz-kernels"></a>
<a name="convex-envelope-for-positive-exponential-mixtures"></a>

For positive mixtures of exponential kernels with integrable tails, the proof also identifies the limit of the energy infima as the even period grows. It equals the lower convex envelope of the corresponding uniform-chain energy evaluated at the mean gap, allowing the value $-\infty$. [Exponential gap estimate](all-period-riesz.md#exponential-gap-estimate-with-unequal-charges) · [Passage to Riesz kernels](all-period-riesz.md#from-exponentials-to-riesz-kernels) · [Convex-envelope proof](all-period-riesz.md#convex-envelope-for-positive-exponential-mixtures)

<a name="example-and-relation-to-prior-work"></a>

Bétermin, Knüpfer and Nolte's Corollary 2.3 proves crystallization under the sufficient condition $m+m^{-1}<2^{p+1}/\zeta(p)$ for $p>p_1\approx1.465$. At $p=2$, the sharp all-period bound here includes $m=5$, which lies outside that sufficient range. The proof extends the neutral exponential gap estimate in the Codex–Jiang draft by retaining the term proportional to $(m-1)^2$. The [priority search](PRIOR-ART.md) records the dated literature and upstream comparison. It found no equivalent finite-period classification in the sources checked. Literature priority remains open. [Corollary 2.3 and its proof](https://arxiv.org/abs/1804.05743) · [Pinned Codex–Jiang draft, Theorem 3](https://github.com/xiujiang4/alternating-particle-chains/blob/a48fe4c2623a29dd45edd227410a07baba400250/en/paper.tex)

## Reproducing the results

Run these commands from the repository root with Python 3. The verification scripts use the standard library.

```sh
python check_collapse_onset.py
python check_charge_windows.py
python check_six_site.py
python check_finite_cells.py
python verify.py
```

- `check_collapse_onset.py` checks 4,078 integer flow records and the negative 77-site witness for the cutoff at 78. [Proof and verifier](collapse-onset.md#verification-and-credit)
- `check_charge_windows.py` checks both rational endpoints of all four charge enclosures using the compressed JSON data. [Endpoint data and method](charge-windows/README.md#reproducing-the-bounds)
- `check_six_site.py` verifies the rational enclosure for $\mu_2$ from the polynomial. The derivation and the proof of the universal lower bound are given separately. [Six-site proof](six-site-threshold.md)
- `check_finite_cells.py` checks the four-site identity, the Hessian thresholds, and the cotangent and rational collapse witnesses. `verify.py` evaluates the gap bound on examples and checks the exact cutoff for uniform compression. [Hessian and trial configurations](inverse-square.md) · [All-period calculations](all-period-riesz.md#how-large-must-a-uniformly-compressed-block-be)

The optional flow generator uses SciPy. Run `python generate_flow.py alternative.json` for the fixed-ratio reference, or `python generate_flow.py alternative.json --reference charge-windows/open-75.json.gz` for one of the variable-ratio references. The stored certificates can be verified independently of how the reference positions and flows were generated. [Flow generator](generate_flow.py)

## Credit and version history

Yingbo Sun (orcarider7) maintains this research. OpenAI Codex contributed to the proofs, code and writing. Astra and Sol independently recomputed the integer checks for the fixed-ratio cutoff. Astra developed the strict-decrease argument and Sol developed the algebraic six-site certificate. Each checked the other's argument and the charge enclosures. Terra reviewed the analytic reduction and the finite-threshold reasoning. Gemini, accessed in Auto mode, reviewed the analytic arguments conditional on the exact computations. The linked notes record the scope of these reviews and the changes they prompted. Formal peer review remains outstanding. [Fixed-ratio review record](collapse-onset.md#verification-and-credit) · [Variable-ratio review record](collapse-ladder.md#attribution-and-checks) · [Citation metadata](CITATION.cff)

The earlier [fixed-ratio cutoff PDF](first-collapse.pdf) and [finite-threshold account PDF](finite-collapse-thresholds.pdf) record the original stages of the work. The [priority record](PRIOR-ART.md#search-record-and-limits) links to the dated releases and pinned source versions. [Companion replication studies](https://github.com/orcarider7/paper-replication-lab)
