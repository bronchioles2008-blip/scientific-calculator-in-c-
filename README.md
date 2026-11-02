# scientific-calculator-in-c-
/*
 * scicalc.c - Casio-style scientific calculator (pure C, console)
 *
 * Build:  gcc -O2 -o scicalc scicalc.c -lm
 * Run:    ./scicalc
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <setjmp.h>
#include <time.h>

#define PI 3.14159265358979323846
#define E_CONST 2.71828182845904523536

typedef enum { DEG, RAD, GRA } Angle;
typedef enum { NORM, FIX, SCI } Disp;

static Angle  amode = DEG;
static Disp   dmode = NORM;
static int    ddigits = 10;
static double ans = 0.0;
static double vars[26];            /* A-Z, M is memory */
static const char *p;              /* parser cursor */
static jmp_buf jb;
static const char *errmsg;

static double parse_expr(void);
static double parse_unary(void);

static void fail(const char *m) { errmsg = m; longjmp(jb, 1); }

static double check(double r) {
    if (isnan(r) || isinf(r)) fail("Math ERROR");
    return r;
}
static void skip(void) { while (*p == ' ' || *p == '\t') p++; }

/* ---------- angle helpers ---------- */
static double to_rad(double x) {
    return amode == DEG ? x * PI / 180.0 : amode == GRA ? x * PI / 200.0 : x;
}
static double from_rad(double x) {
    return amode == DEG ? x * 180.0 / PI : amode == GRA ? x * 200.0 / PI : x;
}
static double snap(double r) { return fabs(r) < 1e-14 ? 0.0 : r; }

/* ---------- number theory ---------- */
static double is_int(double x) { return x == floor(x) && fabs(x) < 9e15; }

static double fact(double x) {
    if (x < 0 || !is_int(x) || x > 170) fail("Math ERROR");
    double r = 1;
    for (int i = 2; i <= (int)x; i++) r *= i;
    return r;
}
static double npr(double n, double r) {
    if (!is_int(n) || !is_int(r) || r < 0 || r > n) fail("Math ERROR");
    double res = 1;
    for (int i = 0; i < (int)r; i++) res *= (n - i);
    return check(res);
}
static double ncr(double n, double r) {
    if (!is_int(n) || !is_int(r) || r < 0 || r > n) fail("Math ERROR");
    double k = r < n - r ? r : n - r, res = 1;
    for (int i = 1; i <= (int)k; i++) res = res * (n - k + i) / i;
    return check(round(res));
}
static double gcd_(double a, double b) {
    if (!is_int(a) || !is_int(b)) fail("Math ERROR");
    long long x = llabs((long long)a), y = llabs((long long)b);
    while (y) { long long t = x % y; x = y; y = t; }
    return (double)x;
}

/* ---------- functions ---------- */
static double call1(const char *f, double x) {
    if (!strcmp(f, "sin"))   return snap(sin(to_rad(x)));
    if (!strcmp(f, "cos"))   return snap(cos(to_rad(x)));
    if (!strcmp(f, "tan")) {
        double r = to_rad(x);
        if (fabs(cos(r)) < 1e-12) fail("Math ERROR");
        return snap(tan(r));
    }
    if (!strcmp(f, "asin"))  return check(snap(from_rad(asin(x))));
    if (!strcmp(f, "acos"))  return check(snap(from_rad(acos(x))));
    if (!strcmp(f, "atan"))  return snap(from_rad(atan(x)));
    if (!strcmp(f, "sinh"))  return check(sinh(x));
    if (!strcmp(f, "cosh"))  return check(cosh(x));
    if (!strcmp(f, "tanh"))  return tanh(x);
    if (!strcmp(f, "asinh")) return asinh(x);
    if (!strcmp(f, "acosh")) return check(acosh(x));
    if (!strcmp(f, "atanh")) return check(atanh(x));
    if (!strcmp(f, "log"))   return check(log10(x));
    if (!strcmp(f, "ln"))    return check(log(x));
    if (!strcmp(f, "exp"))   return check(exp(x));
    if (!strcmp(f, "sqrt"))  return check(sqrt(x));
    if (!strcmp(f, "cbrt"))  return cbrt(x);
    if (!strcmp(f, "sq"))    return check(x * x);
    if (!strcmp(f, "rec"))   { if (x == 0) fail("Math ERROR"); return 1.0 / x; }
    if (!strcmp(f, "abs"))   return fabs(x);
    if (!strcmp(f, "fact"))  return fact(x);
    if (!strcmp(f, "int"))   return trunc(x);
    if (!strcmp(f, "intg"))  return floor(x);
    if (!strcmp(f, "frac"))  return x - trunc(x);
    if (!strcmp(f, "round")) return round(x);
    if (!strcmp(f, "ceil"))  return ceil(x);
    if (!strcmp(f, "floor")) return floor(x);
    if (!strcmp(f, "deg"))   return x * 180.0 / PI;   /* rad -> deg */
    if (!strcmp(f, "rad"))   return x * PI / 180.0;   /* deg -> rad */
    fail("Syntax ERROR");
    return 0;
}

static double call2(const char *f, double a, double b) {
    if (!strcmp(f, "npr"))   return npr(a, b);
    if (!strcmp(f, "ncr"))   return ncr(a, b);
    if (!strcmp(f, "pow"))   { if (a == 0 && b == 0) fail("Math ERROR"); return check(pow(a, b)); }
    if (!strcmp(f, "root"))  { if (a == 0) fail("Math ERROR");
                               if (a == 3) return cbrt(b);
                               return check(pow(b, 1.0 / a)); }
    if (!strcmp(f, "logb"))  { if (a <= 0 || a == 1 || b <= 0) fail("Math ERROR");
                               return check(log(b) / log(a)); }
    if (!strcmp(f, "mod"))   { if (b == 0) fail("Math ERROR"); return a - b * floor(a / b); }
    if (!strcmp(f, "atan2")) return snap(from_rad(atan2(a, b)));
    if (!strcmp(f, "gcd"))   return gcd_(a, b);
    if (!strcmp(f, "lcm"))   { double g = gcd_(a, b); if (g == 0) return 0; return check(fabs(a / g * b)); }
    if (!strcmp(f, "hypot")) return check(hypot(a, b));
    fail("Syntax ERROR");
    return 0;
}

/* ---------- parser ---------- */
static double parse_number(void) {
    char buf[64]; int n = 0;
    while (isdigit((unsigned char)*p) && n < 60) buf[n++] = *p++;
    if (*p == '.') { buf[n++] = *p++; while (isdigit((unsigned char)*p) && n < 60) buf[n++] = *p++; }
    if ((*p == 'e' || *p == 'E') &&
        (isdigit((unsigned char)p[1]) ||
         ((p[1] == '+' || p[1] == '-') && isdigit((unsigned char)p[2])))) {
        buf[n++] = *p++;
        if (*p == '+' || *p == '-') buf[n++] = *p++;
        while (isdigit((unsigned char)*p) && n < 60) buf[n++] = *p++;
    }
    buf[n] = 0;
    if (n == 0 || (n == 1 && buf[0] == '.')) fail("Syntax ERROR");
    return strtod(buf, NULL);
}

static double parse_primary(void) {
    skip();
    if (isdigit((unsigned char)*p) || *p == '.') return parse_number();

    if (*p == '(') {
        p++;
        double v = parse_expr();
        skip();
        if (*p == ')') p++;
        else if (*p) fail("Syntax ERROR");     /* missing ')' at end is allowed */
        return v;
    }

    if (isalpha((unsigned char)*p)) {
        char id[24]; int n = 0;
        while (isalpha((unsigned char)*p) && n < 20) id[n++] = *p++;
        id[n] = 0;
        skip();

        if (n == 1 && isupper((unsigned char)id[0]) && *p != '(')
            return vars[id[0] - 'A'];

        for (char *c = id; *c; c++) *c = (char)tolower((unsigned char)*c);

        if (*p == '(') {                           /* function call */
            p++;
            double a[3]; int cnt = 0;
            skip();
            if (*p == ')') p++;
            else {
                for (;;) {
                    if (cnt >= 3) fail("Syntax ERROR");
                    a[cnt++] = parse_expr();
                    skip();
                    if (*p == ',') { p++; continue; }
                    if (*p == ')') { p++; break; }
                    if (*p == 0) break;            /* auto-close */
                    fail("Syntax ERROR");
                }
            }
            if (!strcmp(id, "rand") && cnt == 0) return (double)rand() / ((double)RAND_MAX + 1.0);
            if (cnt == 1) return call1(id, a[0]);
            if (cnt == 2) return call2(id, a[0], a[1]);
            fail("Syntax ERROR");
        }
        if (!strcmp(id, "pi"))  return PI;
        if (!strcmp(id, "e"))   return E_CONST;
        if (!strcmp(id, "ans")) return ans;
        fail("Syntax ERROR");
    }
    fail("Syntax ERROR");
    return 0;
}

static double parse_postfix(void) {
    double v = parse_primary();
    for (;;) {
        skip();
        if (*p == '!')      { p++; v = fact(v); }
        else if (*p == '%') { p++; v = v / 100.0; }
        else break;
    }
    return v;
}

static double parse_power(void) {
    double b = parse_postfix();
    skip();
    if (*p == '^') {
        p++;
        double e = parse_unary();                  /* right-assoc, allows 2^-3 */
        if (b == 0 && e == 0) fail("Math ERROR");
        return check(pow(b, e));
    }
    return b;
}

static double parse_unary(void) {
    skip();
    if (*p == '-') { p++; return -parse_unary(); }
    if (*p == '+') { p++; return  parse_unary(); }
    return parse_power();
}

static double parse_term(void) {
    double v = parse_unary();
    for (;;) {
        skip();
        if (*p == '*') { p++; v = check(v * parse_unary()); }
        else if (*p == '/') {
            p++;
            double d = parse_unary();
            if (d == 0) fail("Math ERROR");
            v = check(v / d);
        }
        else if (*p == '(' || isalpha((unsigned char)*p) ||
                 isdigit((unsigned char)*p) || *p == '.')
            v = check(v * parse_unary());          /* implicit multiplication */
        else break;
    }
    return v;
}

static double parse_expr(void) {
    double v = parse_term();
    for (;;) {
        skip();
        if (*p == '+')      { p++; v = check(v + parse_term()); }
        else if (*p == '-') { p++; v = check(v - parse_term()); }
        else break;
    }
    return v;
}

static int evaluate(const char *s, double *out) {
    p = s;
    errmsg = NULL;
    if (setjmp(jb)) return 0;
    double v = parse_expr();
    skip();
    if (*p) fail("Syntax ERROR");
    *out = check(v);
    return 1;
}

/* ---------- display ---------- */
static void format(double v, char *buf, size_t n) {
    if (v == 0) v = 0.0;                           /* kill -0 */
    switch (dmode) {
        case FIX: snprintf(buf, n, "%.*f", ddigits, v); break;
        case SCI: snprintf(buf, n, "%.*E", ddigits > 0 ? ddigits - 1 : 0, v); break;
        default:  snprintf(buf, n, "%.10G", v); break;
    }
}

static const char *mode_name(void) {
    return amode == DEG ? "DEG" : amode == RAD ? "RAD" : "GRA";
}

static void help(void) {
    puts(
"Operators : + - * / ^ ! %   ( )   implicit mult: 2pi, 3(4+1), 2sin(30)\n"
"Constants : pi  e  ans  rand()     Variables: A-Z (M = memory)\n"
"Trig      : sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh\n"
"Other     : log ln exp sqrt cbrt sq rec abs fact int intg frac round\n"
"2 args    : npr(n,r) ncr(n,r) pow(a,b) root(n,x) logb(base,x) mod(a,b)\n"
"            atan2(y,x) gcd(a,b) lcm(a,b) hypot(a,b)\n"
"Store     : 5*3 -> A   or   A=5*3        Recall: just use A in an expression\n"
"Memory    : m+ [expr]  m- [expr]  (uses ans if no expr)   mc   vars\n"
"Modes     : deg rad gra | fix N | sci N | norm\n"
"Continue  : a line starting with * / ^ + continues from ans (e.g. '*2')\n"
"Other cmds: ac (clear all)  help  quit");
}

static void trim(char *s) {
    size_t n = strlen(s);
    while (n && isspace((unsigned char)s[n - 1])) s[--n] = 0;
}

static void show_result(double v) {
    char b[64];
    ans = v;
    format(v, b, sizeof b);
    printf("  = %s\n", b);
}

int main(void) {
    char line[512], low[512], expr[540];
    srand((unsigned)time(NULL));

    puts("=== Scientific Calculator (type 'help' for commands, 'quit' to exit) ===");

    for (;;) {
        printf("[%s] > ", mode_name());
        fflush(stdout);
        if (!fgets(line, sizeof line, stdin)) break;
        trim(line);
        char *s = line;
        while (isspace((unsigned char)*s)) s++;
        if (!*s) continue;

        for (size_t i = 0; s[i]; i++) low[i] = (char)tolower((unsigned char)s[i]);
        low[strlen(s)] = 0;

        /* ---- commands ---- */
        if (!strcmp(low, "quit") || !strcmp(low, "exit") || !strcmp(low, "q")) break;
        if (!strcmp(low, "help")) { help(); continue; }
        if (!strcmp(low, "deg"))  { amode = DEG; continue; }
        if (!strcmp(low, "rad"))  { amode = RAD; continue; }
        if (!strcmp(low, "gra"))  { amode = GRA; continue; }
        if (!strcmp(low, "norm")) { dmode = NORM; continue; }
        if (!strcmp(low, "mc"))   { vars['M' - 'A'] = 0; puts("  M cleared"); continue; }
        if (!strcmp(low, "ac"))   { memset(vars, 0, sizeof vars); ans = 0; puts("  cleared"); continue; }
        if (!strcmp(low, "vars")) {
            int any = 0;
            for (int i = 0; i < 26; i++)
                if (vars[i] != 0) { char b[64]; format(vars[i], b, sizeof b); printf("  %c = %s\n", 'A' + i, b); any = 1; }
            if (!any) puts("  (all variables are 0)");
            continue;
        }
        int nd;
        if (sscanf(low, "fix %d", &nd) == 1) { if (nd < 0 || nd > 9) puts("  Argument ERROR"); else { dmode = FIX; ddigits = nd; } continue; }
        if (sscanf(low, "sci %d", &nd) == 1) { if (nd < 1 || nd > 10) puts("  Argument ERROR"); else { dmode = SCI; ddigits = nd; } continue; }

        double v;

        if (low[0] == 'm' && (low[1] == '+' || low[1] == '-')) {
            const char *rest = s + 2;
            while (isspace((unsigned char)*rest)) rest++;
            if (*rest) { if (!evaluate(rest, &v)) { printf("  %s\n", errmsg); continue; } }
            else v = ans;
            vars['M' - 'A'] += (low[1] == '+') ? v : -v;
            char b[64]; format(vars['M' - 'A'], b, sizeof b);
            printf("  M = %s\n", b);
            continue;
        }

        /* ---- store: "expr -> X" ---- */
        char *arrow = strstr(s, "->");
        if (arrow) {
            char *t = arrow + 2;
            while (isspace((unsigned char)*t)) t++;
            if (!isupper((unsigned char)*t) || t[1]) { puts("  Syntax ERROR"); continue; }
            *arrow = 0;
            if (!evaluate(s, &v)) { printf("  %s\n", errmsg); continue; }
            vars[*t - 'A'] = v;
            printf("  %c ", *t); show_result(v);
            continue;
        }

        /* ---- store: "X=expr" ---- */
        if (isupper((unsigned char)s[0])) {
            const char *q = s + 1;
            while (isspace((unsigned char)*q)) q++;
            if (*q == '=') {
                if (!evaluate(q + 1, &v)) { printf("  %s\n", errmsg); continue; }
                vars[s[0] - 'A'] = v;
                printf("  %c ", s[0]); show_result(v);
                continue;
            }
        }

        /* ---- plain expression (with ans-continuation) ---- */
        if (*s == '*' || *s == '/' || *s == '^' || *s == '+') {
            snprintf(expr, sizeof expr, "ans%s", s);
            s = expr;
        }
        if (evaluate(s, &v)) show_result(v);
        else printf("  %s\n", errmsg);
    }
    puts("Bye!");
    return 0;
}

