# The first collapsing period of an alternating inverse-square chain

Research note, 26 November 2026. Maintained by Yingbo Sun (orcarider7).

**The first collapsing even period is 78.** Consider a periodic one-dimensional chain with alternating charges $6,-1$, pair energy $-q_iq_j/r^2$, and a fixed positive cell length. For each even period $N\le76$ the energy attains a finite minimum. For every even $N\ge78$ the infimum is $-\infty$. At $N=78$ equal spacing is still a strict local minimum up to translation. The first even period where equal spacing has a negative second variation is 96. The last two statements follow from the [uniform-chain Hessian](inverse-square.md).

[Download the research note](first-collapse.pdf) · [Integer certificate](flow-75.json) · [Python verifier](check_collapse_onset.py) · [Literature comparison](PRIOR-ART.md)

## 1. Reduction to an open chain

Fix $L>0$, $x_0<\cdots<x_{N-1}<x_0+L$, $x_{i+N}=x_i+L$, and $q_{i+N}=q_i$. The periodic energy per site counts both directions,

$$
E_N=-\frac1N\sum_{i=0}^{N-1}\sum_{j\in\mathbb Z,\ j\ne i}\frac{q_iq_j}{(x_i-x_j)^2}.
$$

Write $S_n=-\sum_{0\le i<j<n}q_iq_j/(x_j-x_i)^2$ for the open energy counted once. Call a site of charge $6$ heavy and one of charge $-1$ light. Cut the periodic chain across its largest gap $g\ge L/N$. The two ends have opposite charges. Read the remaining sites in the direction that starts at a heavy site. Then $E_N=2S_N/N+R$, where $|R|$ is bounded by a constant depending only on $N,L$. The nearest image cells are at distance at least $g$ and the cells $k$ periods away are at distance at least $(|k|-1)L+g$. The square-inverse series converges. This uniform bound continues to hold along any sequence of ordered configurations at fixed $N,L$.

The final site is light. Pair each preceding light with the heavy immediately **to its right**. Seen from the final light, that heavy is closer, so its repulsion exceeds the paired attraction. The initial heavy supplies an additional positive term. In particular $S_N\ge S_{N-1}$. For $N\ge4$ the last such pair contributes at least $5/(x_{N-1}-x_{N-2})^2$, and for $N=2$ the terminal contribution is $6/(x_1-x_0)^2$.

If every heavy-ended open block of length $N-1$ has nonnegative $S_{N-1}$, the periodic energy has a lower bound. Conversely, compress any such block with $S_{N-1}<0$ by a factor $\delta\downarrow0$, leaving a light at a fixed positive distance from the block and its periodic images. Its internal contribution to $E_N$ is $2S_{N-1}/(N\delta^2)$ and all remaining terms stay bounded. Thus these alternatives characterize whether the periodic energy is unbounded below.

## 2. A certificate valid for every choice of the 75 positions

Let $x_0<\cdots<x_{74}$ be arbitrary real positions, beginning and ending with heavy sites. The integers $B_0<\cdots<B_{74}$ in the certificate are auxiliary reference coordinates. For a source $s$ and nodes $s\le u<v$, set

$$
P(u,v)=\frac{(B_v-B_u)^3}{(x_v-x_u)^2},\qquad P(s,s)=0.
$$

Then $P(s,v)-P(s,u)\le P(u,v)$. For $s<u$, apply convexity of the cube on positive numbers to the two slopes $\Delta B/\Delta x$, weighted by their respective $\Delta x$. Their weighted mean is the slope on the entire interval. Multiplying by $x_v-x_s$ gives the inequality. For $s=u$ it is an equality.

The certificate records nonnegative integers $f_{suv}$ on directed edges $s\le u<v$ with $v-u$ odd. All omitted records are zero. Put $K=10^{13}(10^6)^3$. Let $b_{sv}$ be the flow into $v$ minus the flow out of $v$ for source $s$, and $F_{uv}=\sum_s f_{suv}$. The verifier checks the following exact inequalities.

* When $v>s$ has the same parity as $s$, $b_{sv}(B_v-B_s)^3\ge K q_sq_v$. For opposite parity it checks $b_{sv}\ge0$.
* On every mixed edge, $F_{uv}(B_v-B_u)^3\le6K\rho$, where the largest capacity fraction is

$$
\rho=\frac{228876486320393532885019483}{228881835937500000000000000}<1.
$$

To see why this proves a statement about the unknown positions, multiply $P(s,v)-P(s,u)\le P(u,v)$ by $f_{suv}$ and sum. The left side collects by node to $\sum_{s,v>s} b_{sv}P(s,v)$, since $P(s,s)=0$. The right side is $\sum_{u<v}F_{uv}P(u,v)$. Applying the two checks and cancelling $K$ gives

$$
\sum_{\substack{s<v\\v-s\ {\rm even}}}\frac{q_sq_v}{(x_v-x_s)^2}
\ \le\ \rho\sum_{\substack{u<v\\v-u\ {\rm odd}}}\frac{6}{(x_v-x_u)^2}.
$$

Thus $S_{75}$ is at least $1-\rho$ times the positive mixed-pair sum. The check covers 1,369 same-charge pairs, 1,406 mixed edges and 4,078 nonzero flow records. The reference was found by approximately minimizing open energy at fixed width, followed by a linear program for the flows and integer rounding and repair. The inequality above follows entirely from the recorded integers.

Every shorter odd heavy-ended block obeys the same bound. Append alternating sites at $x_{\rm last}+R,x_{\rm last}+2R,\ldots$ until there are 75. All pair terms involving appended sites tend to zero as $R\to\infty$. Taking the limit preserves the bound on the original sites. Section 1 now covers every even $N\le76$. Moreover, as any cyclic gap tends to zero, the positive factor $1-\rho$ controls every interior adjacent gap in the shortened block. The terminal light estimate controls the last gap, and the cut gap is at least $L/N$. Thus $E_N\to+\infty$ at the boundary of the fixed-length gap simplex. Continuity inside and compactness of its closure give attainment of the minimum. This argument leaves the positions of the minimizers undetermined.

## 3. A negative 77-site block

Let $t_j=(j-38)/38$ for $j=0,\ldots,76$ and $g(t)=2t+t^5$. Set $y_j=g(t_j)$ at even $j$, and place each odd site at the midpoint of its two even neighbors. The points are strictly ordered because $g$ is strictly increasing. All $z_j=D y_j$ are integers for $D=2\cdot38^5$. For these 77 sites the verifier checks, using signed integer floor division,

$$
10^6 S_{77}\le-\sum_{i<j}\left\lfloor\frac{10^6D^2 q_iq_j}{(z_j-z_i)^2}\right\rfloor=-36\,340\,297<0.
$$

Compress this block near the origin of a length-$L$ cell and place the remaining $N-77$ alternating sites at distinct fixed positions in $(L/2,3L/4)$, starting with a light. For fixed even $N\ge78$, the distances to all spectators and neighboring images stay bounded away from zero, and the distant image tail converges. The leading energy is $2S_{77}/(N\delta^2)$, proving collapse. Together with Section 2 this establishes the cutoff.

## Verification and credit

Run **python check_collapse_onset.py** next to **flow-75.json**. The checker uses only Python's standard library. To construct another certificate from the same reference with SciPy installed, run **python generate_flow.py alternative.json**. The generator checks every inequality exactly before writing its output. The proof permits different valid flows.

Bétermin, Knüpfer and Nolte developed the periodic two-species model and proved a sufficient crystallization condition. Their corollary and the finite-period result here concern different parameter ranges. This paper studies the case $p=2,m=6$. Generic convexity and network flow arguments have a long history. The dated search found no source giving this exact threshold, and literature priority remains open. [Original article, especially Corollary 2.3](https://doi.org/10.1007/s10955-020-02603-2) · [Search record](PRIOR-ART.md)

Yingbo Sun (orcarider7) maintains the project. AI assistance generated the proof, code and draft. Astra and Sol independently recomputed the integer checks. Terra reviewed the analytic reduction. Gemini, accessed in Auto mode, reviewed the argument conditional on the integer checks and accepted corrections to its initial objections about endpoint pairing and padding. It did not recompute the certificate. This version is a public research manuscript and has not undergone formal journal peer review.
