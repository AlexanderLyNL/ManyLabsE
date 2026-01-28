# Two-Sided Fisher's Exact Test

- By Order: Study questions sequences (we dont care)
- By Site: Different institutes ('uva' or 'uiuc')
- Global: Aggregate all data together

log odds ratio and d
$$
d = \frac{\text{odds ratio}}{\pi / \sqrt{3}} \\
\delta = \log(d \cdot \pi / \sqrt{3})
$$


## [16. Framing decisions (Tversky & Kahneman, 1981) Problem 10](../OSFdata/Framing%20(Tversky%20&%20Kahneman,%201981)/)

**One-sided, Greater**

**$H_0$:** Price of the item does **not** affect travel choice.

**$H_1$:** People travel **more** for the cheaper item.

|           | $30 Item (Vase) | $250 Item (Hanging) | Total |
|-----------|-----------------|---------------------|-------|
| **Yes**   | ya              | yb                  | y1    |
| **No**    | na - ya         | nb - yb             |       |
| **Total** | na              | nb                  | n     |

**Original Study**: Z = 5.14, p = 7.4e-7, OR = 4.96, 95% CI [2.55, 9.90]

**Effect size**: $\delta$ = 1.601406

---

## [5. Affect and Risk (Rottenstreich & Hsee, 2001) Study 1](../OSFdata/Affect%20&%20Risk%20(Rottenstreich%20&%20Hsee,%202001))

**One-sided, Greater**

**$H_0$:** Choice is **independent** of probability level.

**$H_1$:** "Kiss" is preferred **more** at 1% than at 100%.

|           | 1% Chance | 100% Certain | Total |
|-----------|-----------|--------------|-------|
| **Kiss**  | ya        | yb           | y1    |
| **$50**   | na - ya   | nb - yb      |       |
| **Total** | na        | nb           | n     |

**Original Study**: χ2(1, N=40) = 4.91, p = 0.0267, Kramer φ = 0.35, d = 0.74, 95% CI [<0.001, 1.74].

**Effect size**: $\delta$ = 0.2943186

---

## [11. Trolley Dilemma 1 (Hauser et. al. 2007) Scenarios 1+2](../OSFdata/Trolley%20Dilemma%201%20(Hauser%20et%20al.,%202007))

**Two-sided**
the sites (different studies) are column `Source.Secondary`

**$H_0$:** Acting is the **same** for Denise vs Frank

**$H_1$:** Acting is **favored** for Denise over Frank.

|           | Denise  | Frank   | Total |
|-----------|---------|---------|-------|
| **Yes**   | ya      | yb      | y1    |
| **No**    | na - ya | nb - yb |       |
| **Total** | na      | nb      | n     |

**Original Study**: χ2(1, N = 2646) = 1615.96, p < 0.001, w = 0.78, d = 2.50, 95% CI [2.22, 2.86]

**Effect size**: $\delta$ = 1.511714

---

## [17. Trolley Dilemma 2 (Hauser et. al. 2007) Scenarios 3+4](../OSFdata/Trolley%20Dilemma%202%20(Hauser%20et%20al.,%202007))

**One-sided, Greater**

**$H_0$:** Acting is the **same** for Oscar vs Ned.

**$H_1$:** Acting is **favored** for Oscar over Ned.

|           | Oscar   | Ned     | Total |
|-----------|---------|---------|-------|
| **Yes**   | ya      | yb      | y1    |
| **No**    | na - ya | nb - yb |       |
| **Total** | na      | nb      | n     |

**Original study**: χ2(1, N = 2612) = 72.35, p < 0.001, w = 0.17, d = 0.34, 95% CI [0.26, 0.42].

**Effect size**: -0.4833859

![Denise](Hauser_Denise.png)

![Frank](Hauser_Frank.png)

![Ned](Hauser_Ned.png)

![Oscar](Hauser_Oscar.png)
