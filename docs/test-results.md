# Test results

Written by `pytest tests --results-page docs/test-results.md` (in CI on every push to `main`).
Do not edit by hand.

- **Run:** 2026-11-23 16:12 UTC, commit `421234cb92`
- **Environment:** Python 3.13.18, OCP 8.0.1.1, Linux x86_64
- **Tests:** 723 passed, 0 failed, 2 skipped

## Test files

| File | Passed | Failed | Skipped |
|---|---:|---:|---:|
| `test_bendmodel.py` | 14 | 0 | 0 |
| `test_canonical.py` | 43 | 0 | 0 |
| `test_cli.py` | 3 | 0 | 0 |
| `test_corners.py` | 6 | 0 | 0 |
| `test_cut_forms.py` | 18 | 0 | 0 |
| `test_dfm.py` | 31 | 0 | 0 |
| `test_engrave.py` | 31 | 0 | 0 |
| `test_export.py` | 7 | 0 | 0 |
| `test_external.py` | 24 | 0 | 2 |
| `test_features.py` | 5 | 0 | 0 |
| `test_fold.py` | 20 | 0 | 0 |
| `test_formed.py` | 9 | 0 | 0 |
| `test_generated.py` | 38 | 0 | 0 |
| `test_geom2d.py` | 3 | 0 | 0 |
| `test_handworked.py` | 163 | 0 | 0 |
| `test_io_xcaf.py` | 24 | 0 | 0 |
| `test_joggle.py` | 27 | 0 | 0 |
| `test_kinds.py` | 2 | 0 | 0 |
| `test_oracle.py` | 63 | 0 | 0 |
| `test_pmi_map.py` | 90 | 0 | 0 |
| `test_pmi_read.py` | 51 | 0 | 0 |
| `test_property.py` | 1 | 0 | 0 |
| `test_refusals.py` | 5 | 0 | 0 |
| `test_sequence.py` | 8 | 0 | 0 |
| `test_sharp.py` | 19 | 0 | 0 |
| `test_spline.py` | 8 | 0 | 0 |
| `test_tooling.py` | 10 | 0 | 0 |

## SolidCraft reference parts

STEP files and hand-worked flats in [`tests/data/solidcraft_oracle`](../tests/data/solidcraft_oracle/README.md), checked by `tests/test_oracle.py`. Each part unfolded at its K: flat area and size against the hand-worked flat (tolerance 1e-09 relative, 1e-06 mm). *Tests*: the session's tests that name the part, passed / run.

| | Part | What it is | t | Bends | K | Flat area: hand-worked | flatmat | Rel. error | Flat size (mm) | Size error (mm) | Holes | Tests |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| ✅ | `sc67_base_flange` | plate 100 × 60 × 2.5, no bend | 2.5 | 0 | 0.44 | 6000.0000 | 6000.0000 | 0.0e+00 | 60.000 × 100.000 | 0.0e+00 | 0 | 7 / 7 |
| ✅ | `sc68_edge_flange` | one 90° flange, R 2.5, t 2.5 | 2.5 | 1 | 0.44 | 7565.4867 | 7565.4867 | 4.9e-15 | 75.655 × 100.000 | 3.7e-13 | 0 | 7 / 7 |
| ✅ | `sc70_flange_45` | one 45° flange | 2.5 | 1 | 0.44 | 7368.5298 | 7368.5298 | 1.4e-14 | 73.685 × 100.000 | 1.1e-12 | 0 | 7 / 7 |
| ✅ | `sc71_custom_rule` | U channel, t 1.5, R 3, K 0.33 | 1.5 | 2 | 0.33 | 6638.3893 | 6638.3893 | 8.6e-15 | 80.000 × 82.980 | 7.1e-13 | 0 | 7 / 7 |
| ✅ | `sc72_contour_flange` | Z profile, two bends turning opposite ways | 2.5 | 2 | 0.44 | 4815.4867 | 4815.4867 | 7.7e-15 | 50.000 × 96.310 | 7.4e-13 | 0 | 7 / 7 |
| ✅ | `sc73_hem` | flange with a flat hem: 180°, inner radius 0.005, leg 0.01 mm from the flange | 2.5 | 2 | 0.44 | 8462.1327 | 8462.1327 | 7.1e-15 | 84.621 × 100.000 | 6.0e-13 | 0 | 7 / 7 |
| ✅ | `sc74_unfold_refold` | U channel with a 10 × 22 slot cut across a bend in the flat | 2.5 | 2 | 0.44 | 6284.7787 | 6284.7787 | 9.3e-15 | 80.000 × 81.310 | 7.4e-13 | 1 | 7 / 7 |
| ✅ | `sc75_bracket_holes` | two flanges meeting at a mitred corner, two Ø6 holes | 2.5 | 2 | 0.44 | 5441.1148 | 5441.1148 | 5.2e-15 | 65.655 × 90.655 | 3.7e-13 | 2 | 7 / 7 |
| ✅ | `sc76_thickness_variant` | one flange, t 3, R 3 | 3 | 1 | 0.44 | 7478.5840 | 7478.5840 | 5.8e-15 | 74.786 × 100.000 | 4.4e-13 | 0 | 7 / 7 |

## Hand-worked parts: holes across bends, boxes, hems, joggles, countersinks

STEP files and hand-worked flats in [`tests/data/handworked`](../tests/data/handworked/README.md), checked by `tests/test_handworked.py`. Each part unfolded at its K: flat area and size against the hand-worked flat (tolerance 1e-09 relative, 1e-06 mm). *Tests*: the session's tests that name the part, passed / run.

| | Part | What it is | t | Bends | K | Flat area: hand-worked | flatmat | Rel. error | Flat size (mm) | Size error (mm) | Holes | Tests |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| ✅ | `box_closed` | `box_relief` with tabs closing the corners to 1 mm | 2 | 4 | 0.44 | 14890.4990 | 14890.4990 | 2.9e-15 | 129.048 × 149.048 | 6.3e-13 | 0 | 6 / 6 |
| ✅ | `box_deep_closed` | box 60 × 40 × 40, t 1, R 1.5, closed corners 0.5 apart, reliefs 3 | 1 | 4 | 0.44 | 11497.1849 | 11497.1849 | 1.7e-15 | 126.095 × 146.095 | 4.0e-13 | 0 | 6 / 6 |
| ✅ | `box_open` | box 80 × 60 × 30, t 2, R 2: short walls stop 2 short of the long walls' bends | 2 | 4 | 0.44 | 14190.4990 | 14190.4990 | 4.7e-15 | 129.048 × 149.048 | 6.0e-13 | 0 | 6 / 6 |
| ✅ | `box_relief` | `box_open` with 5 × 5 corner reliefs in the base | 2 | 4 | 0.44 | 14050.4990 | 14050.4990 | 4.1e-15 | 129.048 × 149.048 | 6.0e-13 | 0 | 6 / 6 |
| ✅ | `box_thin` | box 120 × 70 × 25, t 1.5, R 1, gap 1, reliefs 3 | 1.5 | 4 | 0.44 | 18732.4282 | 18732.4282 | 2.7e-15 | 125.215 × 175.215 | 3.7e-13 | 0 | 6 / 6 |
| ✅ | `box_window` | `box_relief` with a 20 wide window across a long bend, 10 into the wall, 15 into the base | 2 | 4 | 0.44 | 13460.0211 | 13460.0211 | 3.9e-15 | 129.048 × 149.048 | 6.0e-13 | 1 | 6 / 6 |
| ✅ | `cbore_Z` | Z, t 5: counterbores Ø6.6 / Ø11, 3 deep, on the outer flats on opposite sides; a countersink in the web | 5 | 2 | 0.44 | 7487.9947 | 7487.9947 | 8.6e-15 | 60.000 × 126.336 | 1.1e-12 | 3 | 9 / 9 |
| ✅ | `csk_U` | U channel, t 3: 90° countersinks on the base and both flanges, both sides | 3 | 2 | 0.44 | 10630.1675 | 10630.1675 | 6.8e-15 | 80.000 × 133.572 | 9.1e-13 | 3 | 9 / 9 |
| ✅ | `csk_angles` | L, t 3: 82°, 100° and 120° countersinks | 3 | 1 | 0.44 | 7558.2986 | 7558.2986 | 3.6e-15 | 80.000 × 95.215 | 3.6e-13 | 3 | 9 / 9 |
| ✅ | `hb_notch` | the same L, a 12 deep notch cut in from the edge across the bend | 2 | 1 | 0.44 | 6270.4389 | 6270.4389 | 4.2e-15 | 80.000 × 81.095 | 4.0e-13 | 0 | 6 / 6 |
| ✅ | `hb_radial` | L, t 3, R 8: a Ø10 hole drilled radially through the middle of the bend | 3 | 1 | 0.44 | 4981.5915 | 4981.5915 | 8.4e-12 | 60.000 × 84.640 | 9.7e-13 | 1 | 5 / 5 |
| ✅ | `hb_refold_mid` | L, t 2, R 10: a Ø8 hole cut in the flat at mid-bend, folded at K 0.44 (radial wall) | 2 | 1 | 0.44 | 5175.1504 | 5175.1504 | 2.5e-12 | 60.000 × 87.090 | 1.1e-12 | 1 | 7 / 7 |
| ✅ | `hb_refold_tangent` | L, t 2, R 6: a Ø6 hole cut in the flat on the bend's tangent line, folded at K 0.44 | 2 | 1 | 0.44 | 4820.1504 | 4820.1504 | 1.4e-12 | 60.000 × 80.807 | 7.2e-13 | 1 | 7 / 7 |
| ✅ | `hb_slot` | L, t 2, R 2: an 8 wide slot across the bend, round ends 6 into each flange | 2 | 1 | 0.44 | 3688.9770 | 3688.9770 | 4.1e-15 | 60.000 × 64.524 | 3.0e-13 | 1 | 6 / 6 |
| ✅ | `hb_window` | L, t 2, R 3: a 20 wide window across the bend, 8 into one flange, 10 into the other | 2 | 1 | 0.44 | 6005.6814 | 6005.6814 | 3.9e-15 | 80.000 × 81.095 | 4.0e-13 | 1 | 6 / 6 |
| ✅ | `hem_closed_both` | plate, t 1, closed hems (R 0.1, 0.2 mm gap) at both ends on opposite sides | 1 | 2 | 0.44 | 5003.5752 | 5003.5752 | 2.5e-15 | 60.000 × 83.393 | 2.1e-13 | 0 | 7 / 7 |
| ✅ | `hem_hole_through_both` | L with a hem (R 0.5); a Ø5 hole drilled through the hem leg and the flange behind it | 2 | 2 | 0.44 | 4992.2876 | 4992.2876 | 7.1e-15 | 60.000 × 83.859 | 6.0e-13 | 2 | 10 / 10 |
| ✅ | `hem_open_U` | U channel, t 1.5, an open hem (R 1, 2 mm gap) folded into the channel | 1.5 | 3 | 0.44 | 8800.0707 | 8800.0707 | 7.0e-15 | 80.000 × 110.001 | 7.8e-13 | 0 | 7 / 7 |
| ✅ | `hem_rolled_270` | a rolled edge: one 270° bend, R 3, a 2 mm leg | 1.5 | 1 | 0.44 | 3554.8406 | 3554.8406 | 1.9e-14 | 59.247 × 60.000 | 1.1e-12 | 0 | 7 / 7 |
| ✅ | `joggle_offset` | 45° joggle, 5 mm between its bends, a Ø6 hole either side | 2 | 2 | 0.44 | 5314.8849 | 5314.8849 | 2.5e-14 | 60.000 × 89.524 | 2.2e-12 | 2 | 12 / 12 |
| ✅ | `joggle_on_flange` | L whose flange is joggled: two tangent 20° bends, R 4, no flat between | 2 | 3 | 0.44 | 4975.8466 | 4975.8466 | 1.9e-14 | 60.000 × 82.931 | 1.6e-12 | 0 | 9 / 9 |
| ✅ | `joggle_z_tight` | Z joggle: two 90° bends, R 1, 2 mm apart, t 1.5 | 1.5 | 2 | 0.44 | 3360.7522 | 3360.7522 | 5.3e-15 | 50.000 × 67.215 | 3.6e-13 | 0 | 9 / 9 |
| ✅ | `panel_hem_joggle_csk` | t 1.5: a hem at one end, a joggle at the other, two countersinks between | 1.5 | 3 | 0.44 | 8478.2566 | 8478.2566 | 7.3e-14 | 70.000 × 121.477 | 8.8e-12 | 2 | 10 / 10 |
