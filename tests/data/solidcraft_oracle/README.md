# SolidCraft sheet metal reference parts

The nine sheet metal cases of SolidCraft's analytic reference tests
([tangency-labs/solidcraft](https://github.com/tangency-labs/solidcraft) `tests/oracle/67`
to `76`, MIT OR Apache-2.0), rebuilt in OCP by `tests/oracle_parts.py`. Each part's
volume and area match SolidCraft's hand-worked `expected.json` to 1e-9. Each `<name>.json`
gives the flat pattern worked out by hand: thickness, K, the bends with their allowance
θ(R + K·t), the flat area, size and number of inner loops.

| File | Part | Flat by hand (BA = θ·(R + K·t)) |
|---|---|---|
| `sc67_base_flange` | plate 100 × 60 × 2.5, no bend | 100 × 60 |
| `sc68_edge_flange` | one 90° flange, R 2.5, t 2.5 | 100 × (55 + BA + 15) |
| `sc70_flange_45` | one 45° flange | 100 × (57.929 + BA + 12.929), setback (R + t)·tan 22.5° = 2.071 |
| `sc71_custom_rule` | U channel, t 1.5, R 3, K 0.33 | 80 × (15.5 + 41 + 15.5 + 2 BA) |
| `sc72_contour_flange` | Z profile, two bends turning opposite ways | 50 × (35 + 22.5 + 27.5 + 2 BA) |
| `sc73_hem` | flange with a flat hem: 180°, inner radius 0.005, leg 0.01 mm from the flange | 100 × (55 + 15 + 5.495 + BA₉₀ + BA_hem), BA_hem = π·(0.005 + K·t) |
| `sc74_unfold_refold` | U channel with a 10 × 22 slot cut across a bend in the flat | 80 × (40 + 2·15 + 2 BA) − 10 × (5 + BA + 11.345); the slot is 22 long at K 0.44 |
| `sc75_bracket_holes` | two flanges meeting at a mitred corner, two Ø6 holes | 75 × 35 − 2·π·3² + BA·(75 + 35) + 25 × 75.732 + 10 × 35.732 |
| `sc76_thickness_variant` | one flange, t 3, R 3 | 100 × (54 + BA + 14) |

Regenerate with `python tests/oracle_parts.py tests/data/solidcraft_oracle`;
`tests/test_oracle.py` checks both the generated parts and these files.
