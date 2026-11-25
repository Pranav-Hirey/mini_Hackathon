# RELATIV-GPS: GPS with and without Einstein

**C Programming Mini Hackathon** · Own problem statement (approved by faculty)

| | |
|---|---|
| **Team** | Antifragile |
| **Batch** | 56 |
| **Members** | Pranav Hirey, Parth Chhibber |
| **Repository** | `Batch56_Antifragile_RelativGPS` |

---

## 1. Problem statement

GPS finds your position by measuring **how long** radio signals take to arrive from satellites:

```
distance = speed of light × (time received − time sent)
```

So GPS is really a **clock system**. A clock error of just **1 microsecond** puts you **300 m** off.

GPS satellite clocks are affected by Einstein's relativity:

| Effect | Cause | Change per day |
|---|---|---|
| General relativity | The satellite is higher in Earth's gravity, so its clock runs **faster** | **+45.72 µs** |
| Special relativity | The satellite moves at 3.87 km/s, so its clock runs **slower** | **−7.21 µs** |
| **Net** | | **+38.51 µs** (11.54 km of range error per day) |

Real GPS satellites correct this before launch by setting their clocks slightly slow. **Our program derives that correction from first principles, then switches it off in a simulated GPS world to measure what goes wrong.**

## 2. Objectives

1. Calculate the relativistic clock drift using only physical constants (G·M, c, Earth's radius, orbit radius).
2. Simulate how a GPS receiver finds its position from satellite distances.
3. Compare the position **with** and **without** the relativity correction.
4. Show that a computer can lose the effect entirely if precision is handled badly (`float` vs `double`).
5. Accept user input with full validation, and verify the results automatically.

## 3. Features

The program is **menu-driven**. At start-up it asks for the receiver's location; every option then asks for its own input (pressing **Enter** keeps the default shown in brackets).

| # | Option | User input (range, default) | What it shows |
|---|---|---|---|
| | Start-up | Location: 1–4 preset, 5 custom latitude/longitude (default 1) | Dehradun Bidholi, Nashik, New Delhi, Shimla, or any point on Earth |
| 1 | Clock physics | Satellite altitude, 100–100,000 km (default 20,190 km = GPS) | Gravity (+45.72) and speed (−7.21) effects, net +38.51 µs/day, orbital speed and period, and the factory clock setting (10.22999999544 MHz vs the real 10.22999999543 MHz) |
| 2 | Precision lab | none | In `float`, the effect is computed as **0** (catastrophic cancellation); `double` and a rewritten formula fix it |
| 3 | Mode A | Hours since clock sync, 0–8,760 (default 24) | Receiver that trusts its own clock: **15.653 km** off after 24 h without correction |
| 4 | Mode B | Hours (24); receiver clock error, −1,000 to 1,000 µs (75) | Realistic receiver that also solves for its clock: position **0 m** off, but its time is wrong by 38.51 µs |
| 5 | Mode C | Orbit shape scale, 0–5 (default 1 = real GPS) | Elliptical orbits give each satellite a different clock error: **20.472 m** off, the part real phones correct in software |
| 6 | Error growth | Number of days, 1–365 (default 30) | Table and bar chart: error grows linearly (about 471.5 km after 30 days) |
| 7 | Altitude explorer | Any altitude in km (0 to return) | Drift from the ISS to geostationary orbit; the two effects **cancel at 3,186 km** |
| 8 | Verification panel | none | 9 automatic self-tests with PASS/FAIL |
| 9 | Settings | Location, hours, receiver clock error | Change any setting |

## 4. The physics (summary)

**Clock rate (weak-field limit of general relativity):**
```
dτ/dt ≈ 1 + Φ/c² − v²/(2c²)        Φ = −GM/r,   v = √(GM/r)

Gravity term = (GM/c²)(1/R − 1/r) = 5.29 × 10⁻¹⁰  →  +45.72 µs/day
Speed term   =  GM/(2rc²)          = 8.35 × 10⁻¹¹  →   −7.21 µs/day
Net                                = 4.46 × 10⁻¹⁰  →  +38.51 µs/day
```

**Measured distance (pseudorange):**
```
pseudorange = true distance + c × (receiver clock error) − c × (satellite clock error)
```

**Unknowns:** position (x, y, z) and the receiver's clock error, so at least **4 satellites** are needed.

**Elliptical-orbit offset (Mode C):** `Δt = F · e · √a · sin E`, where `F = −2√(GM)/c²`

**Cancel altitude:** setting gravity term = speed term gives `r = 3R/2`, an altitude of about **3,186 km**.

## 5. Algorithm

**Main program**
1. Start. Set defaults: 24 hours, receiver clock error 75 µs.
2. Ask for the receiver's location (1–5); reject invalid input and ask again.
3. Display the menu (options 0–9).
4. Read the choice; if it is not a whole number from 0 to 9, show an error and ask again.
5. If the choice is 0, display a goodbye message and stop.
6. Otherwise run the chosen option: ask for its input, validate it, calculate, and display the results.
7. Wait for Enter and go back to step 3.

**Simulation (Modes A, B, C)**
1. Set each satellite's clock error: 0 if the correction is ON; the relativistic drift × hours (A, B) or the eccentricity offset (C) if OFF.
2. Place the receiver at its true position (latitude/longitude → x, y, z).
3. Place 4 satellites on the GPS orbit using their azimuth and elevation (solving a quadratic for each distance).
4. Compute each pseudorange: true distance + receiver clock error − c × satellite clock error.
5. Give the solver only the satellite positions and pseudoranges (not the true position).
6. **Solver (Gauss-Newton least squares):** start at Earth's centre; repeat: predict the distances, compare with the pseudoranges, solve the linear equations for a correction (Gaussian elimination with partial pivoting), move the guess; stop when the move is below 1 micrometre.
7. Position error = distance between the solver's answer and the true position. Display the ON and OFF results side by side.

**Flowchart:** see [`flowchart.png`](flowchart.png).

## 6. C programming concepts used

| Concept | Where it is used |
|---|---|
| Variables and data types | `int` for counters and menu choices; `double` for all physics; `float` only in the precision lab; `char` arrays for text |
| Operators | arithmetic, relational, logical (`&&`, `!`), compound assignment (`+=`, `-=`), ternary (`? :`) |
| Conditionals | `if`/`else`, `switch` for the main menu |
| Loops | `for` (satellites, solver, tables), `while` (menu loop, input loop), nested loops (matrix building, elimination) |
| Arrays | 1D (`X[4]`, `pr[NSAT]`), 2D (`H[4][4]`, `N[4][5]`), arrays of structures (`presets[]`, `sky[]`) |
| Strings | `fgets`, `sscanf`, `sprintf`, `strcpy`, `strcspn` |
| Functions | 39 functions (including `main`), each with one job: physics, geometry, solver, input, output, menus |
| Structures | `Vec3`, `Satellite`, `City`, `Drift`, `Solution`, `Settings` (with nested structures) |
| Pointers (basic) | `Settings *` so options can update the settings; output parameters in `local_axes()` and `simulate()`; an optional default (`const double *`, or `NULL`) in `read_number()` |
| Input validation | `read_number()` rejects letters, extra text, out-of-range values, and decimals where a whole number is needed |
| `math.h` | `sqrt`, `sin`, `cos`, `asin`, `atan2`, `fabs`, `floor` |
| Preprocessor | `#include`, `#define` constants |

## 7. How to compile and run

```bash
gcc relativ_gps.c -o relativ_gps -lm
./relativ_gps
```
On Windows: `gcc relativ_gps.c -o relativ_gps.exe`, then `relativ_gps.exe`. (`-lm` links the math library on Linux/macOS.)

## 8. Sample test cases

All results below were produced by running the program (receiver: Dehradun Bidholi).

| # | Type | Test | Input | Expected | Actual |
|---|---|---|---|---|---|
| 1 | Normal | GPS clock drift | Menu 1, Enter | Net about +38.5 µs/day | +38.51 µs/day; 10.22999999544 MHz |
| 2 | Normal | Mode A, one day | Menu 3, `24` | 0 km ON, several km OFF | 0.000 km / 15.653 km |
| 3 | Normal | Mode A, two days | Menu 3, `48` | Error doubles | 31.310 km |
| 4 | Normal | Mode B | Menu 4, `24`, `75` | Position 0 m; bias 75 recovered | 0.000 m; 75.00 µs ON, 36.49 µs OFF |
| 5 | Normal | Mode B, different clock | Menu 4, `48`, `200` | Bias = 200 − 77.02 | 122.99 µs |
| 6 | Normal | Mode C | Menu 5, `1` | Metres of error OFF | 20.472 m |
| 7 | Normal | Error growth | Menu 6, `90` | Linear growth | 1,426.7 km at 90 days |
| 8 | Normal | Low orbit | Menu 1, `1000` | Clock runs slow | −17.83 µs/day |
| 9 | Boundary | No time elapsed | Menu 3, `0` | No error | 0.000 km |
| 10 | Boundary | Maximum time (1 year) | Menu 3, `8760` | Large, finite error | 5,985.278 km |
| 11 | Boundary | Circular orbits | Menu 5, `0` | No error | 0.000 m |
| 12 | Boundary | One day only | Menu 6, `1` | Single row | 15.7 km |
| 13 | Boundary | Cancel altitude | Menu 7 | Net = 0 at 3,186 km | +0.00 µs/day |
| 14 | Boundary | North Pole | Location `5`, `90`, `0` | Solver still converges | 15.653 km; 9/9 PASS |
| 15 | Invalid | Letters in menu | `abc` | Rejected, asks again | `[!] Invalid input` |
| 16 | Invalid | Menu out of range / decimal | `42`, `3.5` | Rejected | `[!] Invalid input` |
| 17 | Invalid | Location out of range | `7` | Rejected | `[!] Invalid input` |
| 18 | Invalid | Latitude / longitude | `95`, `200` | Rejected | `[!] Invalid input` |
| 19 | Invalid | Hours | `-5`, `ten`, `9000` | Rejected | `[!] Invalid input` |
| 20 | Invalid | Receiver clock error | `5000` | Rejected | `[!] Invalid input` |
| 21 | Invalid | Days | `0`, `400`, `2.5` | Rejected | `[!] Invalid input` |
| 22 | Invalid | Altitude too low | `50` | Rejected | `[!] Invalid input` |

Full console output: [`sample_output.txt`](sample_output.txt) · Screenshots: [`screenshots/`](screenshots/)

## 9. Assumptions and limits

- Earth is modelled as a perfect sphere that does not rotate (this explains the tiny gap from the real factory frequency).
- One snapshot of 4 satellites; no atmospheric delay or signal noise.
- Real receivers use 6–12 satellites, live orbit data and further corrections.

## 10. Conclusion

Nanosecond-scale clock physics becomes kilometre-scale position error. A simple receiver would be 15.65 km off after one day; a realistic receiver keeps its position but loses its time; and elliptical orbits still cause about 20 m of error, which is why every phone applies a relativity formula in software. Einstein's relativity is built into the hardware of every GPS satellite and the software of every phone.

## 11. Tools used

An AI assistant was used to help design and explain the program. The team studied, modified and tested all of the code.

## Repository structure

```
Batch56_Antifragile_RelativGPS/
├── README.md
├── relativ_gps.c
├── presentation.pptx
├── flowchart.png
├── sample_output.txt
└── screenshots/
    ├── menu.png  clock.png  precision.png  modea.png  modeb.png
    └── modec.png  growth.png  altitude.png  verify.png  invalid.png
```
