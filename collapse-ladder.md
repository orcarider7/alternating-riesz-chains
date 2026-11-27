# A strictly decreasing sequence of finite collapse thresholds

Research note · Yingbo Sun (orcarider7), repository maintainer · 26 November 2026. [PDF summary](finite-collapse-thresholds.pdf) · [Six-site formula](six-site-threshold.md) · [Charge intervals near six](charge-windows/README.md)

Fix heavy charge $m\ge1$, light charge $-1$, and interaction $-q_iq_j/r^2$. For $k\ge1$ write $S_{2k+1}(x;m)$ for the once-counted open energy of $2k+1$ strictly ordered sites beginning and ending with a heavy. Define
$$
\mu_k=\inf\{m\ge1:\ S_{2k+1}(x;m)<0\ \text{for some }x\},\qquad \mu_0=+\infty .
$$

**Theorem.** The thresholds satisfy $8=\mu_1>\mu_2>\cdots\downarrow 3+2\sqrt2$. For each $k$, every finite configuration has positive energy for $1\le m<\mu_k$. At $m=\mu_k$ every configuration has nonnegative energy and at least one finite configuration has zero energy. For every $m>\mu_k$ a finite configuration has negative energy. Consequently the first collapsing even period is $N=2k+2$ exactly when $\mu_k<m\le\mu_{k-1}$. Here collapse means that the energy at fixed positive cell length has infimum $-\infty$.

### Charge dependence and the limiting threshold

In particular, each even period $N\ge4$ is first to collapse for a nonempty range of charges. For $1\le m\le3+2\sqrt2$ no finite period collapses. The proof below gives existence of all the finite thresholds. The [period-six calculation](six-site-threshold.md) identifies one of them algebraically.

For a fixed shape put $S=mC-m^2A-B$, where $C,A,B$ sum inverse-square distances over mixed, heavy-heavy and light-light pairs respectively. A dyadic subdivision gives a convenient form of the all-period bound for open chains. Write each even index distance as $2^a o$, $o$ odd, $a\ge1$. Split the interval into $2^a$ consecutive pieces of index length $o$. Jensen's inequality bounds its inverse-square distance by $2^{-3a}$ times the sum of inverse-square piece lengths. At a fixed level $a$, a given piece has at most $2^{a-1}$ starting sites of each parity for the original interval, since the possible offsets are consecutive multiples of the odd integer $o$. For $J=\lfloor\log_2(2k)\rfloor$ this gives
$$
m^2 A+B\le(m^2+1)\sum_{a=1}^J 2^{-2a-1} C
=\frac{m^2+1}{6}(1-4^{-J})C .
$$
Write $m_*=3+2\sqrt2$, so $m_*+m_*^{-1}=6$. The strict margin in the finite bound holds on an open neighborhood of $m_*$ and throughout $[1,m_*]$. Thus $\mu_k>m_*$ for every finite $k$.

In particular $S(1)>0$. Since $A>0$, each fixed shape has a unique root above 1, with positive energy from 1 up to that root and negative energy thereafter. The threshold is the infimum of these roots. This proves all the sign assertions except existence of a zero-energy configuration at the threshold. It also shows every threshold is finite.

For equally spaced sites at unit spacing, divide the open energy by $2k+1$ and let $k\to\infty$. For each fixed odd index distance the normalized mixed-pair count tends to 1, and for each even distance each same-parity count tends to $1/2$. The summands are dominated by $m^2/d^2$ at index distance $d$, so the limiting energy per site is
$$
m\sum_{d\ \mathrm{odd}}d^{-2}-\tfrac12(m^2+1)\sum_{d\ \mathrm{even}}d^{-2}
=\frac{\pi^2}{48}(6m-m^2-1),
$$
where both sums run over positive integers. For every $m>m_*$ the limit is negative. Hence $\limsup_k\mu_k\le m_*$, establishing the value of the limit.

### Strict decrease and attainment

For three sites separated by gaps $a,b$, convexity gives $a^{-2}+b^{-2}\ge 8/(a+b)^2$, with equality precisely at $a=b$. Thus the threshold is 8 and equal gaps attain zero energy there.

Assume thresholds through $\mu_k$ are strictly decreasing and a finite shape attains zero at $m=\mu_k$. Translate it so its first site is at zero. Its total charge is $Q=(k+1)m-k$. Append a light at $R$ and a heavy at $cR$ with fixed $c>1$. Expanding each of the finitely many new pair terms gives the energy increment
$$
R^{-2}\left[Q-mQ/c^2+m/(c-1)^2\right]+O(R^{-3}).
$$
Taking $c=Q^{1/3}/(Q^{1/3}-1)$ makes the bracket $Q-m(Q^{1/3}-1)^3$. It is negative if $Q>(1-m^{-1/3})^{-3}$. For $k=1$ we have $m=8$ and $Q=15>8$. For $k\ge2$, $m\ge m_*>27/5$, so $m^{1/3}>7/4$ (since $27/5>343/64$). Therefore $(1-m^{-1/3})^{-3}<(7/3)^3<14$, whereas $Q\ge3m-2>14$. For sufficiently large finite $R$ the extended shape has negative energy at $\mu_k$. By continuity it already has negative energy at a slightly smaller charge. Hence $\mu_{k+1}<\mu_k$.

To prove attainment at $m=\mu_{k+1}$, take a sequence of $(2k+3)$-site shapes with respective roots tending to $m$. Scale each shape to have smallest gap 1, preserving its zero at its own root. Pass to a subsequence along which the index of a smallest gap is fixed and every gap converges in $[1,+\infty]$. The energies evaluated at $m$ tend to zero: all inverse-square pair terms are at most 1, so the charge change from the corresponding root changes the energy by a quantity tending uniformly to zero.

If any gap tends to infinity, the limiting energy is the sum of open energies of the resulting proper contiguous clusters, since all cross-cluster terms vanish. Each singleton contributes zero, and at least one cluster contains a gap of length 1. Each nonsingleton proper cluster has strictly positive energy, as follows. A light endpoint makes a strictly positive contribution to interactions with all remaining sites: moving inward, pair each light with the preceding, closer heavy, with a possible unpaired heavy contributing positively. The inequality uses $m\ge1$ and strictly ordered sites. Removing light endpoints thus costs strictly positive energy. The remaining heavy-ended core is a singleton of energy zero or has length $2j+1$ with $1\le j\le k$. Such a core has positive energy since $m<\mu_k\le\mu_j$. The same argument covers clusters that already have two heavy endpoints. The positive limiting sum contradicts convergence to zero. All gaps therefore have finite limits, giving a strictly ordered zero-energy shape at $m$. This completes the induction.

### Periodic consequence

For $N=2k+2$, the [largest-gap reduction](collapse-onset.md) applies for every $m\ge1$. Cut across a gap at least $L/N$ and orient the cell from a heavy to a light. Image terms outside the open cell are uniformly bounded, and the terminal light contributes positively by the pairing argument above. Thus nonnegativity of the heavy-ended $(N-1)$-site energy bounds the periodic energy below. A negative heavy-ended block can instead be compressed with the remaining light at fixed positive separation to force collapse. Period 2 is bounded below directly by the same cut. This proves the claimed first-period classification, including the upper endpoint of each charge interval. At $m=\mu_k$ the period $2k+2$ has finite infimum. Boundedness here does not settle whether that periodic infimum is attained.

### Attribution and checks

The periodic model comes from Bétermin, Knüpfer and Nolte. The finite-size questions here extend the project's earlier exact onset at charge ratio six. Codex generated the arguments, code and draft. Astra developed the strict-decrease proof, and Sol developed the algebraic six-site certificate. They independently checked each other's arguments and the four rational enclosures. Terra and Gemini reviewed the analytic reasoning. Gemini was accessed in Auto mode and did not recompute the data files. A proposed monotonicity shortcut failed on a five-site configuration, so the proof uses the two roots of the fixed-shape energy quadratic. Formal peer review and independent confirmation of literature priority remain open. [Original article](https://doi.org/10.1007/s10955-020-02603-2) · [Literature search](PRIOR-ART.md)
