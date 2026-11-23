# Holes across bends, box corners, hems, joggles and countersinks

Sheet metal parts built in OCP by `tests/handworked_parts.py`, each with its flat pattern
worked out by hand in `<name>.json`: thickness, K (0.44), the bends with their allowance
θ(R + K·t), the flat area, size and number of inner loops. `tests/test_handworked.py`
checks the generated parts at K 0.33, 0.44 and 0.5, and these files at K 0.44.
A countersink or counterbore is checked for its kind, side, diameters, angle and depth;
the flat cuts its through hole, so flat area × t is not the volume there.

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
| `hem_open_U` | U channel, t 1.5, an open hem (R 1, 2 mm gap) folded into the channel | width × (flats + 2 BA + BA_hem) |
| `hem_closed_both` | plate, t 1, closed hems (R 0.1, 0.2 mm gap) at both ends on opposite sides | width × (flats + 2 BA_hem) |
| `hem_rolled_270` | a rolled edge: one 270° bend, R 3, a 2 mm leg | width × (40 + 2 + BA_270) |
| `hem_hole_through_both` | L with a hem (R 0.5); a Ø5 hole drilled through the hem leg and the flange behind it | two Ø5 holes BA_hem + 10 apart |
| `joggle_offset` | 45° joggle, 5 mm between its bends, a Ø6 hole either side | holes 20 + BA + 5 + BA + 20 apart along, 20 across |
| `joggle_z_tight` | Z joggle: two 90° bends, R 1, 2 mm apart, t 1.5 | width × (flats + 2 BA) |
| `joggle_on_flange` | L whose flange is joggled: two tangent 20° bends, R 4, no flat between | width × (flats + 3 BA) |
| `csk_U` | U channel, t 3: 90° countersinks on the base and both flanges, both sides | through holes cut; centres at s = flats + BA + f·flat |
| `csk_angles` | L, t 3: 82°, 100° and 120° countersinks | as `csk_U`; each angle reported |
| `cbore_Z` | Z, t 5: counterbores Ø6.6 / Ø11, 3 deep, on the outer flats on opposite sides; a countersink in the web | as `csk_U`; depth reported |
| `panel_hem_joggle_csk` | t 1.5: a hem at one end, a joggle at the other, two countersinks between | as `csk_U` |

The JSON of the parts with holes gives each hole's centre as (s along the flat, z across):
the tests compare their distances to each other, which do not depend on where the flat is
placed. `features` lists the countersinks and counterbores, `side` +1 for one skin and −1
for the other (which one flatmat calls `top` is its reference skin's choice).

A hole cut in the flat and then bent has a radial wall: a ruled B-spline surface between
the hole's images on the two skins. Flatmat reports it as `SURFACE_CONVERTED` (its
straight edges are lines) and unfolds it to the circle within 1e-7 mm.

Regenerate with `python tests/handworked_parts.py tests/data/handworked`.
