//
//RELATIV-GPS : GPS positioning with and without Einstein
// C Programming Mini Hackathon
// Derives the relativistic drift of GPS satellite clocks from
// dtau/dt = 1 + Phi/c^2 - v^2/(2c^2)
//and feeds it into a GPS position solver, correction ON vs OFF

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Physical constants (SI units) 
#define PI          3.14159265358979323846
#define C_LIGHT     299792458.0      //speed of light, m/s        
#define GM_EARTH    3.986004418e14   // G x Earth mass, m^3/s^2    
#define R_EARTH     6371000.0        // mean Earth radius, m       
#define A_GPS       26561750.0       //GPS orbit radius, m        
#define SEC_PER_DAY 86400.0
#define F0_GPS      10.23e6          // nominal GPS clock freq, Hz 

#define NSAT  4                      // satellites in the solver   
#define NCITY 4                      // preset receiver locations  

/* ---------------- Structures ---------------- */
typedef struct { double x, y, z; } Vec3;       // Earth-centred, m 

typedef struct {
    int    prn;                    // satellite ID                 
    double az_deg, el_deg;         //position in the sky          
    double ecc, ecc_anom_deg;      //orbit eccentricity, angle E  
    Vec3   pos;
} Satellite;

typedef struct { char name[32]; double lat_deg, lon_deg; } City;

typedef struct { double gr, sr, net; } Drift;  //fractional rate change 

typedef struct {
    Vec3   pos;                    //estimated receiver position  
    double bias_m;                 // estimated clock bias x c 
    double residual_m;             //RMS misfit of the ranges     
    int    iterations, ok;
} Solution;

typedef struct {
    City   city;
    double hours;                  //time since last clock sync   
    double bias_us;                //receiver's own clock error   
} Settings;

City presets[NCITY] = {
    {"Dehradun Bidholi",   30.4155, 77.9683},
    {"Nashik",             19.9975, 73.7898},
    {"New Delhi",          28.6139, 77.2090},
    {"Shimla",             31.1046, 77.1734}
};

/* PRN, azimuth, elevation, eccentricity, eccentric anomaly */
Satellite sky[NSAT] = {
    { 5,   0.0, 72.0, 0.010,  30.0, {0, 0, 0}},
    {12,  95.0, 38.0, 0.004, 250.0, {0, 0, 0}},
    {19, 215.0, 44.0, 0.018, 110.0, {0, 0, 0}},
    {27, 318.0, 33.5, 0.007, 320.0, {0, 0, 0}}
};

/* ---------------- Vector helpers ---------------- */
Vec3   vec(double x, double y, double z) { Vec3 v = {x, y, z}; return v; }
Vec3   vadd(Vec3 a, Vec3 b)     { return vec(a.x + b.x, a.y + b.y, a.z + b.z); }
Vec3   vsub(Vec3 a, Vec3 b)     { return vec(a.x - b.x, a.y - b.y, a.z - b.z); }
Vec3   vscale(Vec3 a, double k) { return vec(a.x * k, a.y * k, a.z * k); }
double vdot(Vec3 a, Vec3 b)     { return a.x * b.x + a.y * b.y + a.z * b.z; }
double vnorm(Vec3 a)            { return sqrt(vdot(a, a)); }
double rad(double deg)          { return deg * PI / 180.0; }
double us_per_day(double frac)  { return frac * SEC_PER_DAY * 1e6; }
double km_of_us(double us)      { return us * 1e-6 * C_LIGHT / 1000.0; }

Vec3 geo_to_xyz(double lat_deg, double lon_deg)
{
    double la = rad(lat_deg), lo = rad(lon_deg);
    return vec(R_EARTH * cos(la) * cos(lo),
               R_EARTH * cos(la) * sin(lo),
               R_EARTH * sin(la));
}

double height_of(Vec3 p)
{
    double h = vnorm(p) - R_EARTH;
    return fabs(h) < 1e-6 ? 0.0 : h;      /* hide round-off like -0.000 */
}

/* ---------------- Physics ---------------- */
/* Satellite clock rate relative to a clock on the ground:
 *   gravity term (GR) = (GM/c^2)(1/R - 1/r)   -> runs faster
 *   speed term   (SR) = v^2/(2c^2), v^2=GM/r   -> runs slower   */
Drift relativity_drift(double r_orbit)
{
    Drift d;
    double c2 = C_LIGHT * C_LIGHT;
    d.gr  = (GM_EARTH / c2) * (1.0 / R_EARTH - 1.0 / r_orbit);
    d.sr  = (GM_EARTH / r_orbit) / (2.0 * c2);
    d.net = d.gr - d.sr;
    return d;
}

/* Extra clock offset from an elliptical orbit (seconds):
 *   dt = F e sqrt(a) sin(E),  F = -2 sqrt(GM)/c^2               */
double eccentricity_offset(double ecc, double ecc_anom_deg)
{
    double F = -2.0 * sqrt(GM_EARTH) / (C_LIGHT * C_LIGHT);
    return F * ecc * sqrt(A_GPS) * sin(rad(ecc_anom_deg));
}

/* ---------------- Satellite geometry ---------------- */
/* East, north and up unit vectors at a point on Earth */
void local_axes(City c, Vec3 *east, Vec3 *north, Vec3 *up)
{
    double la = rad(c.lat_deg), lo = rad(c.lon_deg);
    *east  = vec(-sin(lo), cos(lo), 0.0);
    *north = vec(-sin(la) * cos(lo), -sin(la) * sin(lo), cos(la));
    *up    = vec(cos(la) * cos(lo), cos(la) * sin(lo), sin(la));
}

/* Put each satellite on the GPS orbit sphere, in the direction
 * (azimuth, elevation) seen from the receiver.                  */
void build_constellation(Satellite sats[], City c)
{
    Vec3 e, n, u, rx = geo_to_xyz(c.lat_deg, c.lon_deg);
    int i;
    local_axes(c, &e, &n, &u);

    for (i = 0; i < NSAT; i++) {
        double az = rad(sats[i].az_deg), el = rad(sats[i].el_deg);
        Vec3 dir = vadd(vadd(vscale(e, cos(el) * sin(az)),
                             vscale(n, cos(el) * cos(az))),
                        vscale(u, sin(el)));
        /* solve |rx + t*dir| = A_GPS for the distance t */
        double b = vdot(rx, dir);
        double t = -b + sqrt(b * b - vdot(rx, rx) + A_GPS * A_GPS);
        sats[i].pos = vadd(rx, vscale(dir, t));
    }
}

/* ---------------- Numerics ---------------- */
/* Gaussian elimination with partial pivoting.
 * A is an n x (n+1) augmented matrix, n <= 4. Returns 0 if singular. */
int solve_linear(double A[4][5], int n, double x[])
{
    int i, j, k, p;
    for (k = 0; k < n; k++) {
        p = k;
        for (i = k + 1; i < n; i++)
            if (fabs(A[i][k]) > fabs(A[p][k])) p = i;
        if (fabs(A[p][k]) < 1e-12) return 0;
        for (j = 0; j <= n; j++) {
            double t = A[k][j]; A[k][j] = A[p][j]; A[p][j] = t;
        }
        for (i = k + 1; i < n; i++) {
            double f = A[i][k] / A[k][k];
            for (j = k; j <= n; j++) A[i][j] -= f * A[k][j];
        }
    }
    for (i = n - 1; i >= 0; i--) {
        double s = A[i][n];
        for (j = i + 1; j < n; j++) s -= A[i][j] * x[j];
        x[i] = s / A[i][i];
    }
    return 1;
}

/* GPS solver (Gauss-Newton least squares).
 * Unknowns: x, y, z (+ clock bias if solve_bias). Start at Earth's
 * centre, compare predicted and measured ranges, correct, repeat. */
Solution solve_position(Satellite sats[], double pr[], int solve_bias)
{
    Solution s;
    double X[4] = {0, 0, 0, 0}, sum = 0;
    int n = solve_bias ? 4 : 3, it, i, j, k;

    s.ok = 0;
    s.iterations = 0;
    for (it = 1; it <= 30; it++) {
        double H[NSAT][4], res[NSAT], N[4][5], dx[4], step = 0;
        Vec3 p = vec(X[0], X[1], X[2]);

        for (i = 0; i < NSAT; i++) {
            Vec3 d = vsub(sats[i].pos, p);
            double r = vnorm(d);
            res[i] = pr[i] - (r + X[3]);
            H[i][0] = -d.x / r;  H[i][1] = -d.y / r;  H[i][2] = -d.z / r;
            H[i][3] = 1.0;
        }
        /* normal equations: (H^T H) dx = H^T res */
        for (j = 0; j < n; j++) {
            for (k = 0; k <= n; k++) {
                N[j][k] = 0;
                for (i = 0; i < NSAT; i++)
                    N[j][k] += H[i][j] * (k < n ? H[i][k] : res[i]);
            }
        }
        if (!solve_linear(N, n, dx)) break;

        for (j = 0; j < n; j++) { X[j] += dx[j]; step += dx[j] * dx[j]; }
        s.iterations = it;
        if (sqrt(step) < 1e-6) { s.ok = 1; break; }
    }

    s.pos = vec(X[0], X[1], X[2]);
    s.bias_m = X[3];
    for (i = 0; i < NSAT; i++) {
        double r = vnorm(vsub(sats[i].pos, s.pos)) + s.bias_m;
        sum += (pr[i] - r) * (pr[i] - r);
    }
    s.residual_m = sqrt(sum / NSAT);
    return s;
}

/* Simulate the measurements and solve.
 * A satellite clock running AHEAD by dt stamps a later send time,
 * so its signal looks shorter by c*dt.                           */
Solution simulate(Settings *st, double sat_clock_err_s[], double rx_bias_m,
                  int solve_bias, double *error_m)
{
    Satellite sats[NSAT];
    double pr[NSAT];
    Vec3 truth = geo_to_xyz(st->city.lat_deg, st->city.lon_deg);
    Solution s;
    int i;

    memcpy(sats, sky, sizeof sats);
    build_constellation(sats, st->city);
    for (i = 0; i < NSAT; i++)
        pr[i] = vnorm(vsub(sats[i].pos, truth)) + rx_bias_m
                - C_LIGHT * sat_clock_err_s[i];

    s = solve_position(sats, pr, solve_bias);
    *error_m = vnorm(vsub(s.pos, truth));
    return s;
}

void fill(double arr[], double value)
{
    int i;
    for (i = 0; i < NSAT; i++) arr[i] = value;
}

/* ---------------- Input and validation ---------------- */
int valid_range(double v, double lo, double hi) { return v >= lo && v <= hi; }

/* Reads one number. Rejects letters, extra text and out-of-range
 * values; whole_only also rejects decimals (used for menus).
 * If def is not NULL, pressing just Enter returns *def.           */
double read_number(const char *prompt, double lo, double hi, int whole_only,
                   const double *def)
{
    char buf[128], extra;
    double v;
    while (1) {
        printf("%s", prompt);
        fflush(stdout);
        if (fgets(buf, sizeof buf, stdin) == NULL) { printf("\n"); exit(0); }
        buf[strcspn(buf, "\r\n")] = '\0';          /* remove the Enter key */
        if (buf[0] == '\0' && def != NULL)            /* only Enter pressed  */
            return *def;
        if (sscanf(buf, "%lf %c", &v, &extra) == 1 && valid_range(v, lo, hi)
            && (!whole_only || v == floor(v)))
            return v;
        if (whole_only)
            printf("  [!] Invalid input. Enter a whole number from %g to %g.\n", lo, hi);
        else
            printf("  [!] Invalid input. Enter a number from %g to %g.\n", lo, hi);
    }
}

/* Asks for a value, showing its range and default (Enter keeps it) */
double ask(const char *label, double lo, double hi, int whole_only, double def)
{
    char prompt[160];
    sprintf(prompt, "  %s (%g to %g) [Enter = %g]: ", label, lo, hi, def);
    return read_number(prompt, lo, hi, whole_only, &def);
}

void pause_screen(void)
{
    char buf[16];
    printf("\n  Press Enter to return to the menu...");
    fflush(stdout);
    if (fgets(buf, sizeof buf, stdin) == NULL) exit(0);
}

/* ---------------- Console output ---------------- */
void rule(char ch)
{
    int i;
    printf(" ");
    for (i = 0; i < 62; i++) putchar(ch);
    printf("\n");
}

void header(const char *title, const char *sub)
{
    printf("\n");
    rule('=');
    printf("   %s\n", title);
    if (sub[0]) printf("   %s\n", sub);
    rule('=');
}

void print_latlon(Vec3 p)
{
    double lat = asin(p.z / vnorm(p)) * 180.0 / PI;
    double lon = atan2(p.y, p.x) * 180.0 / PI;
    printf("%7.4f%c %8.4f%c", fabs(lat), lat >= 0 ? 'N' : 'S',
           fabs(lon), lon >= 0 ? 'E' : 'W');
}

/* ---------------- 1. Clock physics ---------------- */
void menu_clock_physics(void)
{
    double alt_km = ask("Satellite altitude in km", 100, 100000, 0,
                        (A_GPS - R_EARTH) / 1000);
    double r = R_EARTH + alt_km * 1000;
    Drift d = relativity_drift(r);
    double v = sqrt(GM_EARTH / r);
    double T = 2.0 * PI * sqrt(r * r * r / GM_EARTH);

    header("CLOCK PHYSICS | derived from first principles",
           "dtau/dt = 1 + Phi/c^2 - v^2/(2c^2)   (weak-field metric)");
    printf(" CONSTANTS + YOUR ALTITUDE (every result below is derived)\n");
    printf("   GM = %.9e m^3/s^2    c = %.0f m/s\n", GM_EARTH, C_LIGHT);
    printf("   Earth radius = %.0f km    Orbit radius = %.2f km (altitude %.0f km)\n",
           R_EARTH / 1000, r / 1000, alt_km);
    rule('-');
    printf(" DERIVED ORBIT\n");
    printf("   Orbital speed v = sqrt(GM/r)        %8.1f m/s\n", v);
    printf("   Orbital period  = 2*pi*sqrt(r^3/GM) %8.3f h\n", T / 3600);
    rule('-');
    printf(" EINSTEIN'S TWO EFFECTS\n");
    printf("   Gravity (GR)  dPhi/c^2   = %+.4e -> %+7.2f us/day (faster)\n",
           d.gr, us_per_day(d.gr));
    printf("   Speed   (SR)  v^2/2c^2   = %+.4e -> %+7.2f us/day (slower)\n",
           -d.sr, -us_per_day(d.sr));
    rule('-');
    printf("   NET                      = %+.4e -> %+7.2f us/day\n",
           d.net, us_per_day(d.net));
    printf("   x speed of light         = %.2f km of range error per day\n",
           km_of_us(us_per_day(d.net)));
    rule('-');
    printf(" THE FIX: set the satellite clock's rate before launch\n");
    printf("   Nominal frequency          %.11f MHz\n", F0_GPS / 1e6);
    printf("   Factory setting needed     %.11f MHz\n", F0_GPS * (1.0 - d.net) / 1e6);
    printf("   Real GPS factory setting   10.22999999543 MHz (at 20,190 km)\n");
    printf("   (tiny gap: we ignore Earth's rotation and flattening)\n");
    pause_screen();
}

/* ---------------- 2. Precision lab ---------------- */
void menu_precision_lab(void)
{
    double b2 = (GM_EARTH / A_GPS) / (C_LIGHT * C_LIGHT);   /* v^2/c^2 */
    float  b2f = (float)b2;
    float  f_naive  = 1.0f - sqrtf(1.0f - b2f);
    float  f_series = b2f / 2.0f + b2f * b2f / 8.0f;
    double d_naive  = 1.0 - sqrt(1.0 - b2);
    double d_series = b2 / 2.0 + b2 * b2 / 8.0;              /* reference */

    header("PRECISION LAB | can a computer even see Einstein?",
           "Speed effect = 1 - sqrt(1 - v^2/c^2),  v^2/c^2 = 1.67e-10");
    printf("                        RESULT           us/day    ERROR\n");
    rule('-');
    printf(" float,  naive formula  %.6e   %7.3f   %s\n", f_naive,
           us_per_day(f_naive), f_naive == 0.0f ? "100%  EINSTEIN VANISHED" : "");
    printf(" double, naive formula  %.6e   %7.3f   %.1e\n", d_naive,
           us_per_day(d_naive), fabs(d_naive - d_series) / d_series);
    printf(" float,  series formula %.6e   %7.3f   %.1e\n", f_series,
           us_per_day(f_series), fabs(f_series - d_series) / d_series);
    printf(" double, series formula %.6e   %7.3f   reference\n", d_series,
           us_per_day(d_series));
    rule('-');
    printf("  A float keeps ~7 digits: 1 - 0.000000000167 rounds to exactly 1.\n");
    printf("  This is CATASTROPHIC CANCELLATION. Fix 1: use double (~16 digits).\n");
    printf("  Fix 2: 1 - sqrt(1-x) = x/2 + x^2/8 has no subtraction at all.\n");
    pause_screen();
}

/* ---------------- 3. Mode A: trusted receiver clock ---------------- */
void menu_mode_a(Settings *st)
{
    Drift d = relativity_drift(A_GPS);
    double clk[NSAT], err_on, err_off, clk_us;
    Solution on, off;
    char sub[96];

    st->hours = ask("Hours since last clock sync", 0, 8760, 0, st->hours);
    clk_us = us_per_day(d.net) * st->hours / 24.0;

    fill(clk, 0.0);
    on = simulate(st, clk, 0.0, 0, &err_on);
    fill(clk, clk_us * 1e-6);
    off = simulate(st, clk, 0.0, 0, &err_off);

    sprintf(sub, "Receiver: %s | Mode A | %.0f h since sync", st->city.name, st->hours);
    header("RELATIV-GPS | trusted receiver clock (x, y, z)", sub);
    printf(" Satellite clocks gain %.2f us/day (GR %+.2f, SR %+.2f)\n\n",
           us_per_day(d.net), us_per_day(d.gr), -us_per_day(d.sr));
    printf("                       CORRECTION ON       CORRECTION OFF\n");
    printf(" Clock error           %8.2f us          %8.2f us\n", 0.0, clk_us);
    printf(" Range error           %8.3f km          %8.3f km\n", 0.0, km_of_us(clk_us));
    printf(" Position            ");
    print_latlon(on.pos);
    printf("  ");
    print_latlon(off.pos);
    printf("\n");
    printf(" Height error          %8.3f km          %8.3f km\n",
           height_of(on.pos) / 1000, height_of(off.pos) / 1000);
    printf(" POSITION ERROR        %8.3f km          %8.3f km\n",
           err_on / 1000, err_off / 1000);
    printf(" Solver                %d iter, %.0e m    %d iter, fit %.1f km\n",
           on.iterations, on.residual_m, off.iterations, off.residual_m / 1000);
    printf(" STATUS                [ LOCKED ]          [ OFF TARGET ]\n");
    rule('-');
    printf(" Satellite geometry turns the range error into a bigger, mostly\n");
    printf(" vertical, position error.\n");
    pause_screen();
}

/* ---------------- 4. Mode B: realistic receiver ---------------- */
void menu_mode_b(Settings *st)
{
    Drift d = relativity_drift(A_GPS);
    double clk[NSAT], err_on, err_off, drift_us, rx_bias_m;
    Solution on, off;
    char sub[96];

    st->hours   = ask("Hours since last clock sync", 0, 8760, 0, st->hours);
    st->bias_us = ask("Receiver clock error in us", -1000, 1000, 0, st->bias_us);
    drift_us  = us_per_day(d.net) * st->hours / 24.0;
    rx_bias_m = st->bias_us * 1e-6 * C_LIGHT;

    fill(clk, 0.0);
    on = simulate(st, clk, rx_bias_m, 1, &err_on);
    fill(clk, drift_us * 1e-6);
    off = simulate(st, clk, rx_bias_m, 1, &err_off);

    sprintf(sub, "Receiver: %s | Mode B | %.0f h | true rx bias %.2f us",
            st->city.name, st->hours, st->bias_us);
    header("RELATIV-GPS | realistic receiver (x, y, z, clock bias)", sub);
    printf("                        CORRECTION ON     CORRECTION OFF\n");
    printf(" Satellite drift        %8.2f us        %8.2f us\n", 0.0, drift_us);
    printf(" Position error         %8.3f m         %8.3f m\n", err_on, err_off);
    printf(" Estimated clock bias   %8.2f us        %8.2f us\n",
           on.bias_m / C_LIGHT * 1e6, off.bias_m / C_LIGHT * 1e6);
    printf(" Bias change                --          %+8.2f us (%.2f km)\n",
           (off.bias_m - on.bias_m) / C_LIGHT * 1e6, fabs(off.bias_m - on.bias_m) / 1000);
    rule('-');
    printf("  A drift shared by all satellites looks like receiver clock error,\n");
    printf("  so the 4th unknown absorbs it: position right, TIME wrong by\n");
    printf("  %.2f us. Only DIFFERENCES between satellites move you (Mode C).\n", drift_us);
    pause_screen();
}

/* ---------------- 5. Mode C: eccentric orbits ---------------- */
void menu_mode_c(Settings *st)
{
    double clk[NSAT], err_on, err_off, scale, ecc;
    double rx_bias_m = st->bias_us * 1e-6 * C_LIGHT;
    char sub[96];
    int i;

    scale = ask("Orbit shape scale, 0 = circular, 1 = real GPS", 0, 5, 0, 1);

    sprintf(sub, "Receiver: %s | Mode C | eccentricity x %.2f", st->city.name, scale);
    header("RELATIV-GPS | eccentric orbits (the effect that survives)", sub);
    printf("  dt = F e sqrt(a) sin E\n");
    printf("  PRN   eccentricity  E (deg)   clock offset   range effect\n");
    for (i = 0; i < NSAT; i++) {
        ecc = sky[i].ecc * scale;
        clk[i] = eccentricity_offset(ecc, sky[i].ecc_anom_deg);
        printf("  %3d      %.3f      %5.0f      %+7.2f ns      %+6.2f m\n",
               sky[i].prn, ecc, sky[i].ecc_anom_deg,
               clk[i] * 1e9, -clk[i] * C_LIGHT);
    }
    simulate(st, clk, rx_bias_m, 1, &err_off);
    fill(clk, 0.0);
    simulate(st, clk, rx_bias_m, 1, &err_on);
    rule('-');
    printf("                        CORRECTION ON     CORRECTION OFF\n");
    printf(" Position error         %8.3f m         %8.3f m\n", err_on, err_off);
    printf(" STATUS                 [ LOCKED ]        [ DRIFTING ]\n");
    rule('-');
    printf("  Each clock is off by a DIFFERENT amount, so the clock bias cannot\n");
    printf("  absorb it. Real receivers apply this formula in software.\n");
    pause_screen();
}

/* ---------------- 6. Error growth ---------------- */
void menu_growth(Settings *st)
{
    int steps[10] = {1, 2, 3, 7, 14, 30, 60, 90, 180, 365};
    int days[11], n = 0, max_days, i, j, bars;
    double err_km[11], clk[NSAT], clk_us, err;
    double net_us = us_per_day(relativity_drift(A_GPS).net);
    char sub[96];

    max_days = (int)ask("Show error growth up to how many days", 1, 365, 1, 30);
    for (i = 0; i < 10 && steps[i] < max_days; i++)   /* rows below the limit */
        days[n++] = steps[i];
    days[n++] = max_days;                             /* last row = user's value */

    for (i = 0; i < n; i++) {
        fill(clk, net_us * days[i] * 1e-6);
        simulate(st, clk, 0.0, 0, &err);
        err_km[i] = err / 1000;
    }
    sprintf(sub, "Receiver: %s | Mode A | correction OFF", st->city.name);
    header("ERROR GROWTH | what if nobody corrected for Einstein?", sub);
    printf("  days  clock err   range err  position err\n");
    for (i = 0; i < n; i++) {
        clk_us = net_us * days[i];
        bars = (int)(err_km[i] / err_km[n - 1] * 24 + 0.5);
        printf("  %3d   %7.1f us  %7.1f km  |", days[i], clk_us, km_of_us(clk_us));
        for (j = 0; j < 24; j++) putchar(j < bars ? '#' : ' ');
        printf("| %6.1f km\n", err_km[i]);
    }
    rule('-');
    printf("  Grows LINEARLY: the clocks keep ticking at the wrong rate.\n");
    pause_screen();
}

/* ---------------- 7. Altitude explorer ---------------- */
void print_alt_row(const char *label, double alt_km, double max_us)
{
    Drift d = relativity_drift(R_EARTH + alt_km * 1000);
    double net = us_per_day(d.net);
    int i, n = (int)(fabs(net) / max_us * 14 + 0.5);

    printf(" %-10s %7.0f %7.2f %7.2f %+7.2f  ", label, alt_km,
           us_per_day(d.gr), -us_per_day(d.sr), net);
    for (i = 14; i > 0; i--) putchar(net < 0 && i <= n ? '#' : ' ');
    putchar('|');
    for (i = 1; i <= 14; i++) putchar(net > 0 && i <= n ? '#' : ' ');
    printf("\n");
}

void menu_altitude(void)
{
    /* GR = SR  ->  1/R - 1/r = 1/(2r)  ->  r = 3R/2 */
    double alt_cancel = 0.5 * R_EARTH / 1000;
    double max_us = us_per_day(relativity_drift(R_EARTH + 35786e3).net);
    double a;

    header("ALTITUDE EXPLORER | where do Einstein's effects cancel?",
           "net = GM/c^2 (1/R - 1/r) - GM/(2rc^2)   [us/day]");
    printf(" orbit        alt km    GR      SR      NET   slower <-|-> faster\n");
    rule('-');
    print_alt_row("ISS", 408, max_us);
    print_alt_row("Low orbit", 1000, max_us);
    print_alt_row("", 2000, max_us);
    print_alt_row("CANCEL", alt_cancel, max_us);
    print_alt_row("", 5000, max_us);
    print_alt_row("Medium", 10000, max_us);
    print_alt_row("GPS", (A_GPS - R_EARTH) / 1000, max_us);
    print_alt_row("Geostat.", 35786, max_us);
    rule('-');
    printf("  GR = SR gives r = 3R/2, so at %.0f km the clock keeps ground\n", alt_cancel);
    printf("  time. Below it astronauts AGE SLOWER than us; above it, faster.\n");
    rule('-');
    while ((a = read_number("  Try your own altitude in km (0 to go back): ", 0, 1e6, 0, NULL)) != 0)
        print_alt_row("Yours", a, max_us);
}

/* ---------------- 8. Verification panel ---------------- */
int check(const char *name, int pass, const char *fmt, double value)
{
    printf(" [%s] %-30s -> ", pass ? "PASS" : "FAIL", name);
    printf(fmt, value);
    printf("\n");
    return pass;
}

void menu_verify(Settings *st)
{
    Drift d = relativity_drift(A_GPS);
    double clk[NSAT], err, bias_us;
    double gr = us_per_day(d.gr), sr = us_per_day(d.sr);
    double cancel = us_per_day(relativity_drift(1.5 * R_EARTH).net);
    Solution s;
    int passed = 0, total = 9;

    header("VERIFICATION PANEL | the program tests itself", st->city.name);
    passed += check("1 us clock error = range", fabs(1e-6 * C_LIGHT - 299.79) < 0.01,
                    "%.1f m", 1e-6 * C_LIGHT);
    passed += check("GR effect in known range", valid_range(gr, 45.5, 45.9), "%.2f us/day", gr);
    passed += check("SR effect in known range", valid_range(sr, 7.0, 7.4), "%.2f us/day", sr);
    passed += check("Effects cancel at r = 3R/2", fabs(cancel) < 1e-9, "net %.1e us/day", cancel);

    fill(clk, 0.0);
    s = simulate(st, clk, 0.0, 0, &err);
    passed += check("Mode A, correction ON", err < 1e-3, "error %.1e m", err);
    passed += check("Solver residual < 1 mm", s.ok && s.residual_m < 1e-3,
                    "residual %.1e m", s.residual_m);

    s = simulate(st, clk, 75e-6 * C_LIGHT, 1, &err);
    bias_us = s.bias_m / C_LIGHT * 1e6;
    passed += check("Mode B recovers 75 us bias", fabs(bias_us - 75) < 1e-6, "%.4f us", bias_us);

    fill(clk, d.net * SEC_PER_DAY);
    simulate(st, clk, 75e-6 * C_LIGHT, 1, &err);
    passed += check("Mode B absorbs common drift", err < 1e-3, "error %.1e m", err);

    passed += check("Rejects bad input (lat 95)", !valid_range(95, -90, 90), "rejected", 0);

    rule('-');
    printf("  %d / %d checks passed  %s\n", passed, total,
           passed == total ? "[ ALL GOOD ]" : "[ CHECK FAILURES ]");
    pause_screen();
}

/* ---------------- 9. Settings ---------------- */
void choose_city(Settings *st)
{
    int i, pick;
    printf("\n  Where is the GPS receiver?\n");
    for (i = 0; i < NCITY; i++) printf("    %d. %s\n", i + 1, presets[i].name);
    printf("    %d. Custom latitude/longitude\n", NCITY + 1);
    pick = (int)ask("Pick a location", 1, NCITY + 1, 1, 1);

    if (pick <= NCITY) {
        st->city = presets[pick - 1];
    } else {
        strcpy(st->city.name, "Custom location");
        st->city.lat_deg = ask("Latitude in degrees", -90, 90, 0, st->city.lat_deg);
        st->city.lon_deg = ask("Longitude in degrees", -180, 180, 0, st->city.lon_deg);
    }
}

void menu_settings(Settings *st)
{
    int ch;
    while (1) {
        header("SETTINGS", "");
        printf("  Receiver : %s (", st->city.name);
        print_latlon(geo_to_xyz(st->city.lat_deg, st->city.lon_deg));
        printf(")\n");
        printf("  Duration : %.1f hours since last clock sync\n", st->hours);
        printf("  Rx bias  : %.2f us (receiver's own clock error)\n", st->bias_us);
        rule('-');
        printf("  1. Change receiver location\n");
        printf("  2. Change duration\n");
        printf("  3. Change receiver clock error\n");
        printf("  0. Back\n");
        ch = (int)read_number("  Choice: ", 0, 3, 1, NULL);

        if (ch == 0) return;
        if (ch == 1) choose_city(st);
        if (ch == 2) st->hours = ask("Hours since last clock sync", 0, 8760, 0, st->hours);
        if (ch == 3) st->bias_us = ask("Receiver clock error in us", -1000, 1000, 0, st->bias_us);
    }
}

/* ---------------- Main menu ---------------- */
int main(void)
{
    Settings st;
    int choice;

    st.city = presets[0];
    st.hours = 24.0;
    st.bias_us = 75.0;

    header("RELATIV-GPS | GPS with and without Einstein",
           "A C simulation of relativity in satellite navigation");
    choose_city(&st);

    while (1) {
        header("RELATIV-GPS | GPS with and without Einstein", "");
        printf("  Receiver: %s  |  %.0f h since sync\n", st.city.name, st.hours);
        rule('-');
        printf("  1. Clock physics    : derive the drift from first principles\n");
        printf("  2. Precision lab    : float vs double (where Einstein vanishes)\n");
        printf("  3. Mode A           : trusted receiver clock (3 unknowns)\n");
        printf("  4. Mode B           : realistic receiver (4 unknowns)\n");
        printf("  5. Mode C           : eccentric orbits (what still breaks)\n");
        printf("  6. Error growth     : days without correction\n");
        printf("  7. Altitude explorer: where GR and SR cancel\n");
        printf("  8. Verification panel\n");
        printf("  9. Settings\n");
        printf("  0. Exit\n");
        rule('-');
        choice = (int)read_number("  Enter choice: ", 0, 9, 1, NULL);

        switch (choice) {
            case 1: menu_clock_physics(); break;
            case 2: menu_precision_lab(); break;
            case 3: menu_mode_a(&st);     break;
            case 4: menu_mode_b(&st);     break;
            case 5: menu_mode_c(&st);     break;
            case 6: menu_growth(&st);     break;
            case 7: menu_altitude();      break;
            case 8: menu_verify(&st);     break;
            case 9: menu_settings(&st);   break;
            case 0:
                printf("\n  Nanoseconds of clock physics -> kilometres of position.\n");
                printf("  Thank you, Einstein.\n\n");
                return 0;
        }
    }
}
