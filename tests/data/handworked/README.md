# Holes across bends and box corners

Sheet metal parts built in OCP by `tests/handworked_parts.py`, each with its flat pattern
worked out by hand in `<name>.json`: thickness, K (0.44), the bends with their allowance
θ(R + K·t), the flat area, size and number of inner loops. `tests/test_handworked.py`
checks the generated parts at K 0.33, 0.44 and 0.5, and these files at K 0.44.

| File | Part | Flat by hand |
|---|---|---|
| `hb_window` | L, t 2, R 3: a 20 wide window across the bend, 8 into one flange, 10 into the other | hole 20 × (8 + BA + 10) |
| `hb_notch` | the same L, a 12 deep notch cut in from the edge across the bend | outline less 12 × (6 + BA + 6) |
| `hb_slot` | L, t 2, R 2: an 8 wide slot across the bend, round ends 6 into each flange | 8 × (12 + BA) + π·4² |
| `hb_radial` | L, t 3, R 8: a Ø10 hole drilled radially through the middle of the bend | 4·rn·r_s·(E(m) − (1 − m)·K(m)), m = (a/r_s)²; depends on the skin developed (both given) |
| `hb_refold_mid` | L, t 2, R 10: a Ø8 hole cut in the flat at mid-bend, folded at K 0.44 (radial wall) | the Ø8 circle at K 0.44; π·4²·(R + K·t)/(R + 0.44·t) at another K |
| `hb_refold_tangent` | L, t 2, R 6: a Ø6 hole cut in the flat on the bend's tangent line, folded at K 0.44 | the Ø6 circle at K 0.44: half arc, half spline |
| `box_open` | box 80 × 60 × 30, t 2, R 2: short walls stop 2 short of the long walls' bends | base + 4 bends + 4 walls |
| `box_relief` | `box_open` with 5 × 5 corner reliefs in the base | less 4 × 5 × (2 + 5) |
| `box_closed` | `box_relief` with tabs closing the corners to 1 mm | plus the 8 tabs |
| `box_thin` | box 120 × 70 × 25, t 1.5, R 1, gap 1, reliefs 3 | as `box_relief` |
| `box_deep_closed` | box 60 × 40 × 40, t 1, R 1.5, closed corners 0.5 apart, reliefs 3 | as `box_closed` |
| `box_window` | `box_relief` with a 20 wide window across a long bend, 10 into the wall, 15 into the base | less 20 × (10 + BA + 15) |

A hole cut in the flat and then bent has a radial wall: a ruled B-spline surface between
the hole's images on the two skins. Flatmat reports it as `SURFACE_CONVERTED` (its
straight edges are lines) and unfolds it to the circle within 1e-7 mm.

Regenerate with `python tests/handworked_parts.py tests/data/handworked`.
