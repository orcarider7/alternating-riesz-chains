# Which period collapses first near charge ratio six?

In the inverse-square chain with alternating charges $m,-1$ and $m\ge1$, fix the cell length and allow the sites to move while preserving their order. The first collapsing period is the smallest even cell population for which the energy has infimum $-\infty$. The following ranges include both printed endpoints.

| Charge ratio $m$ | First collapsing period | Example |
| --- | ---: | ---: |
| $5.9921304\le m\le5.9961805$ | 80 | 5.995 |
| $5.9961806\le m\le6.0004364$ | 78 | 6 |
| $6.0004365\le m\le6.0049142$ | 76 | 6.002 |

The gaps between rows are each $10^{-7}$ wide. The [strict threshold theorem](../collapse-ladder.md) identifies one transition in each gap. Its location is enclosed to the precision below.

| Period $N$ | Open sites $n=N-1$ | Lower charge | Upper charge | Integer flows | Upper bound on $10^{15}S_n$ at the upper charge |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 74 | 73 | 6.0049142 | 6.0049143 | 4,577 | −4,620,644,506 |
| 76 | 75 | 6.0004364 | 6.0004365 | 4,920 | −5,216,550,989 |
| 78 | 77 | 5.9961805 | 5.9961806 | 5,071 | −10,000,453,899 |
| 80 | 79 | 5.9921303 | 5.9921304 | 5,446 | −5,592,164,452 |

## What the two endpoints prove

The transition for period $N$ is $\mu_{(N-2)/2}$ in the threshold theorem. At the lower charge, we use a certificate valid for **every ordered open configuration** of $n$ sites. In the [flow proof](../collapse-onset.md#2-a-certificate-valid-for-every-choice-of-the-75-positions), replace 6 by $m$ on mixed edges and 36 by $m^2$ at heavy-heavy sinks. At rational $m=p/d$ the [verifier](../check_collapse_onset.py) multiplies through by the denominators and computes the largest capacity fraction $\rho$ exactly. In each case $\rho<1$, so every strictly ordered configuration has positive open energy. Since a finite open configuration has energy zero at its critical charge, the lower endpoint is strictly below the threshold.

At the upper charge, take the reference integers $B_j$ from the same file as actual positions $x_j=B_j/U$, with $U=10^{10}$. Put $q'_j=p$ for heavy sites and $q'_j=-d$ for light sites. The table bounds the once-counted open energy by signed integer floor division:

$$
10^{15} S_n(x;m)
\le -\sum_{i<j}\left\lfloor\frac{10^{15} U^2 q'_iq'_j}{d^2(B_j-B_i)^2}\right\rfloor<0.
$$

A negative value at that charge persists at some slightly smaller charge by continuity. Thus the threshold is strictly inside the printed enclosure. When one period has already crossed its threshold and the preceding even period is still below its own, the strict order of all thresholds identifies the first collapsing period. This gives the closed ranges in the first table.

## Reproducing the bounds

From the repository root run **python check_charge_windows.py**. Python's standard library suffices. The [four compressed JSON files](./) each contain integer coordinates, scales, positive integer flows and the two exact rational charges. The program verifies the reference size and order, every required sink balance, every mixed-edge capacity, and each signed energy bound. It checks the endpoints against the stated rational numbers.

To find candidates, we fixed the endpoints at $0,n-1$, imposed reflection symmetry for the search, and solved the interior force equations and $S_n=0$ jointly in the remaining coordinates and charge. Coordinates were rounded to the reference grid. Reflection symmetry is used only to find witnesses. The lower certificates apply to asymmetric positions as well. An optional sparse linear program generates alternative flows from a stored reference with SciPy. For example, run **python generate_flow.py alternate.json --reference charge-windows/open-75.json.gz**. The choice of reference coordinates or optimizer need not be trusted to check either side of the enclosure.
