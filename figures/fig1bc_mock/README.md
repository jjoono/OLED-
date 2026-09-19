# Fig. 1(b) and 1(c) mock-up

Both panels are eq. (2) and take no input but the single-pass escape probability p and the
round-trip loss A′. `fig1bc.py` draws them; there is no simulation behind it.

(b) bars are the fraction escaping at each pass, p[(1−p)(1−A′)]ⁿ, on the left axis; lines are
the running total on the right axis, and n = 0 is the first pass, before any round trip. Each
line saturates at η_ext, and those three values are the three marked points in (c), so the two
panels share the same three cases and the same p.

(c) η_ext against R_LED = 1 − A′ at p = 0.4, with p = 0.25 and 0.6 in grey to show that p only
shifts the family. Two points are named, the Al electrode at A′ = 0.25 and a lossless reflector
at 0.02; A′ = 0.1 is drawn unnamed because no device in this work sits there.

| A′ | R_LED | η_ext |
|---|---|---|
| 0.02 | 0.98 | 0.97 |
| 0.10 | 0.90 | 0.87 |
| 0.25 | 0.75 | 0.73 |

No material names appear in (b): Al and Ag are introduced in Fig. 2, and at this point the
comparison is meant to stay qualitative. p = 0.4 is the escape probability of ordinary glass,
roughly 1/n_sub², as the text states.
