# The collapse threshold for period six

Consider a periodic chain with six sites, alternating charges $m,-1$ with $m\ge 1$, pair energy $-q_iq_j/r^2$, and fixed positive cell length. There is a critical ratio $m_*$ with

$$
7.44624222875 < m_* < 7.44624223510.
$$

The periodic energy attains a finite minimum for $1\le m<m_*$. It remains bounded below at $m=m_*$, and its infimum is $-\infty$ for $m>m_*$. The exact ratio is

$$
m_* = \frac{16(1+3u_*)}{9(1-u_*)^3},
$$

where $u_*$ is the unique root in $(1/4,13/50)$ of

$$
256u^2(3+u)(1+3u)-128(1+3u)(1-u)^3+9(1-u)^6=0.
$$

Run `python check_six_site.py` for a rational verification of the enclosure.

## Reduction to five sites

For a heavy-ended open block of five sites, let $R$ be the sum of the inverse squares of the mixed-pair distances, $H$ the corresponding heavy-pair sum, and $L$ the light-pair sum. The open energy counted once is $S_5=mR-m^2H-L$. Cutting the periodic chain across its largest gap leaves an open chain of six sites, with one heavy and one light end. The remaining periodic terms are uniformly bounded at fixed cell length. This is the [same cut used for the ratio-six result](collapse-onset.md#1-reduction-to-an-open-chain).

Keeping the heavy end and omitting the terminal light leaves the block of five. Seen from the terminal light, each preceding light has a heavy site immediately to its right that is closer. Since $m\ge 1$, each such pair has nonnegative interaction with the terminal light. The first heavy contributes an additional positive term. Thus the open energy of all six sites is at least $S_5$. Conversely, a block with $S_5<0$ can be compressed inside the period, keeping the remaining light at a fixed positive distance from the block and its periodic images. The internal energy then tends to $-\infty$.

## A stationary block of zero energy

Take reference positions $(0,1-t,1,1+t,2)$ for $0<t<1$. On this symmetric family, put $u=t^2$ and write its energy divided by $m$ as

$$
F(m,u)=\frac{4(1+u)}{(1-u)^2}+\frac{2-1/(4m)}{u}-\frac{9m}{4}.
$$

For $m\ge6$, this function is strictly convex in $u$, diverges at both endpoints, and has negative derivative at $u=1/4$ because $2-1/(4m)>52/27$. Its unique minimizing $u$ therefore exceeds $1/4$. The stationarity equation below expresses $u$ as the inverse of a continuous, strictly increasing function evaluated at $2-1/(4m)$. Thus the minimum depends continuously on $m$. It decreases strictly for $m\ge 6$ since $\partial_m F=1/(4m^2u)-9/4<0$ for $u\ge1/4$. At $m=6$ the minimum is positive. Indeed, the first two terms in the equivalent expression $2/(1-t)^2+(2-1/(4m))/t^2+2/(1+t)^2-9m/4$ sum to at least $14$ by convexity of the inverse square. At $m=8$, evaluating at $u=1/4$ gives a negative value. Hence there is a unique stationary zero with $6<m_*<8$ and $u_*>1/4$.

The stationarity condition is $2-1/(4m)=4u^2(3+u)/(1-u)^3$. Substituting it into $F=0$ gives the formula for $m_*$, and eliminating $m$ gives the displayed polynomial. For $u>1/4$, the formula for $m$ is increasing and gives $m>6$. Thus every root in $(1/4,1)$ corresponds to the unique stationary zero just established. The verifier checks opposite polynomial signs at $1/4$ and $13/50$.

## A global certificate at the critical ratio

Set $t=\sqrt{u_*}$, $m=m_*$ and let $B=(0,1-t,1,1+t,2)$ denote the reference positions. For each same-charge pair $i<j$, assign a demand $q_iq_j/(B_j-B_i)^3$ and route it from $i$ to $j$ along strictly increasing paths of mixed edges. The capacity of a mixed edge $(a,b)$ is $m/(B_b-B_a)^3$. Use the following paths.

- The pairs $(0,2)$, $(2,4)$ and $(1,3)$ each use the path through their single intervening site.
- The pair $(0,4)$ sends $m/(1+t)^3$ along each of the paths $(0,3,4)$ and $(0,1,4)$. It sends the remaining $m^2/8-2m/(1+t)^3$ along $(0,1,2,3,4)$. This is positive since $t>1/2$ and $m>6>128/27$.

Both mixed edges of index span three are now saturated. To check the adjacent edges, write the energy in the four independent gaps. At the reference, reflection makes the first and fourth gap derivatives equal and the second and third equal. Stationarity in $t$ makes these two derivative values equal because the symmetric variation has gap vector $(-1,1,1,-1)$. Euler's identity for a homogeneous function of degree $-2$ gives that all four derivatives vanish, since the energy is zero and the gaps sum to two.

For each gap, its energy derivative is $-2$ times the sum of mixed capacities crossing the gap minus the sum of same-charge demands crossing it, with the difference taken before multiplying by $-2$. The flow load crossing a gap equals that sum of demands because each routed path crosses precisely the gaps between its endpoints. All crossing mixed edges other than the adjacent edge are already saturated, so the vanishing derivative forces the adjacent edge to be saturated as well.

For any actual positions $x_0<\cdots<x_4$ and any route from $i$ to $j$, convexity of the cube gives

$$
\frac{(B_j-B_i)^3}{(x_j-x_i)^2}
\le
\sum_{(a,b)\text{ on route}}\frac{(B_b-B_a)^3}{(x_b-x_a)^2}.
$$

The convexity weights are the actual edge lengths divided by $x_j-x_i$. Multiply by each routed amount and sum. Every same-charge demand is served and the mixed loads equal the capacities, so $m^2H+L\le mR$ for every five-site block. Equality holds at the reference. Moreover, the path $(0,1,2,3,4)$ carries strictly positive flow. Equality in its strictly convex cube inequality forces all four ratios of reference gap to actual gap to agree. Thus every zero-energy open five-site block at the critical ratio is a translate and positive rescaling of the displayed configuration.

## Below and above the threshold

At $m=1$, each same-charge pair $(i,j)$ can be compared to the two shorter mixed edges $(i,j-1)$ and $(i+1,j)$. Both maps are injective, so averaging the corresponding bounds spends at most one unit of each mixed edge. For the span-two same-charge pairs, convexity gives the sharper bound $(a+b)^{-2}\le(a^{-2}+b^{-2})/8$ in terms of their two adjacent gaps $a,b$. Each adjacent edge belongs to at most two such pairs, and the longer same-charge pair uses only the span-three mixed edges. Consequently $S_5$ at $m=1$ is at least $3/4$ times the sum of the inverse squares of the four gaps. For fixed positions, $S_5/m=R-mH-L/m$ is concave in $m>0$. Interpolation between $1$ and $m_*$ gives a positive multiple of the same gap bound for each $1\le m<m_*$. The terminal-light comparison also controls its final adjacent gap since its nearest heavy contributes, after pairing with the preceding light, at least $(m-1)$ times the inverse square of that gap. At $m=1$, writing the final gap as $b$ and the previous gap as $a$, this pair contributes $b^{-2}-(a+b)^{-2}$. Together with the five-site control of $a^{-2}$ this diverges if $b\to0$. Hence the periodic energy diverges to $+\infty$ at the boundary of the fixed-length gap simplex for $1\le m<m_*$, which gives attainment by continuity and compactness.

For $m>m_*$ use the same reference block. Its $F(m,u_*)$ is strictly decreasing for $m\ge m_*$ and is zero at $m_*$. Thus the block has negative energy and yields collapse by compression.
