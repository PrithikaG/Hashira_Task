import java.io.*;
import java.math.BigInteger;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Main {

    // ---- Early-invalid shares (bad digits/base) ----
    static final Set<Integer> EARLY_BAD = new HashSet<>();

    // ---- Exact Rational arithmetic with BigInteger ----
    static final class Rational {
        BigInteger num, den; // den > 0
        Rational(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new IllegalArgumentException("Zero denominator");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            if (!g.equals(BigInteger.ONE)) { n = n.divide(g); d = d.divide(g); }
            this.num = n; this.den = d;
        }
        static Rational of(BigInteger n) { return new Rational(n, BigInteger.ONE); }
        static Rational of(long n) { return of(BigInteger.valueOf(n)); }
        Rational add(Rational o) {
            BigInteger g = den.gcd(o.den);
            BigInteger b1 = den.divide(g);
            BigInteger d2 = o.den.divide(g);
            BigInteger nsum = num.multiply(d2).add(o.num.multiply(b1));
            BigInteger dsum = b1.multiply(o.den);
            return new Rational(nsum, dsum);
        }
        Rational sub(Rational o) { return add(new Rational(o.num.negate(), o.den)); }
        Rational mul(Rational o) {
            BigInteger g1 = num.gcd(o.den);
            BigInteger g2 = o.num.gcd(den);
            BigInteger n1 = num.divide(g1);
            BigInteger d1 = den.divide(g2);
            BigInteger n2 = o.num.divide(g2);
            BigInteger d2 = o.den.divide(g1);
            return new Rational(n1.multiply(n2), d1.multiply(d2));
        }
        Rational div(Rational o) { return new Rational(num.multiply(o.den), den.multiply(o.num)); }
        boolean isZero() { return num.signum() == 0; }
        boolean isInteger() { return den.equals(BigInteger.ONE); }
        BigInteger asInteger() {
            if (!isInteger()) throw new IllegalStateException("Not an integer");
            return num;
        }
        @Override public boolean equals(Object x) {
            if (!(x instanceof Rational)) return false;
            Rational r = (Rational)x;
            return num.equals(r.num) && den.equals(r.den);
        }
        @Override public int hashCode() { return num.hashCode()*31 + den.hashCode(); }
        @Override public String toString() { return den.equals(BigInteger.ONE) ? num.toString() : (num + "/" + den); }
    }

    static final class Share {
        int x;               // x coordinate (from JSON key)
        int base;            // radix declared
        String raw;          // value digits in 'base'
        BigInteger y;        // parsed y
    }

    public static void main(String[] args) throws Exception {
        String json = readAll();
        int n = extractInt(json, "\"n\"\\s*:\\s*(\\d+)");
        int k = extractInt(json, "\"k\"\\s*:\\s*(\\d+)");

        List<Share> shares = parseShares(json);

        // If too few valid shares, report gracefully and exit
        if (shares.size() < k) {
            List<Integer> wrong = new ArrayList<>(EARLY_BAD);
            Collections.sort(wrong);
            int wrongCount = wrong.size();
            int accuracy100 = (int)Math.round(100.0 * (n - wrongCount) / Math.max(1, n));
            System.out.println("Secret: (undetermined)");
            System.out.println("Reason: Not enough valid shares to reconstruct (have " + shares.size() + ", need " + k + ").");
            System.out.println("Wrong shares (x-ids): " + wrong);
            System.out.println("Accuracy (shares consistent) : " + accuracy100 + "/100");
            System.out.println("Consensus confidence (k-subset agreement) : 0/100");
            return;
        }

        // 1) Compute c0 (constant) over all size-k subsets via Gaussian elimination on Vandermonde
        Map<Rational, Integer> freq = new HashMap<>();
        List<int[]> kCombos = combinations(shares.size(), k);
        for (int[] idx : kCombos) {
            List<Share> sub = subset(shares, idx);
            try {
                List<Rational> coeffs = solveVandermonde(sub); // c0..c_{k-1}
                Rational c0 = coeffs.get(0);
                freq.merge(c0, 1, Integer::sum);
            } catch (RuntimeException e) {
                // singular/degenerate (shouldn’t happen with distinct x), skip
            }
        }
        Map.Entry<Rational,Integer> modalEntry = freq.entrySet().stream()
                .max(Comparator.comparingInt(Map.Entry::getValue)).orElseThrow();
        Rational modal = modalEntry.getKey();
        int modalCount = modalEntry.getValue();
        int totalCombos = kCombos.size();

        if (!modal.isInteger()) {
            // Extremely unlikely with integer shares; handle gracefully
            System.out.println("Secret: (undetermined)");
            System.out.println("Reason: Modal secret is not an integer: " + modal);
            List<Integer> wrong = new ArrayList<>(EARLY_BAD);
            Collections.sort(wrong);
            System.out.println("Wrong shares (x-ids): " + wrong);
            System.out.println("Accuracy (shares consistent) : " +
                    (int)Math.round(100.0 * (n - wrong.size()) / Math.max(1, n)) + "/100");
            System.out.println("Consensus confidence (k-subset agreement) : " +
                    (int)Math.round(100.0 * modalCount / Math.max(1, totalCombos)) + "/100");
            return;
        }
        BigInteger secret = modal.asInteger();

        // 2) Choose one good k-subset that yields modal secret; use it to evaluate all shares
        List<Rational> goodCoeffs = null;
        for (int[] idx : kCombos) {
            List<Share> sub = subset(shares, idx);
            try {
                List<Rational> coeffs = solveVandermonde(sub);
                if (coeffs.get(0).equals(modal)) { goodCoeffs = coeffs; break; }
            } catch (RuntimeException ignored) {}
        }
        if (goodCoeffs == null) {
            System.out.println("Secret: (undetermined)");
            System.out.println("Reason: Could not find a consistent k-subset matching modal secret.");
            List<Integer> wrong = new ArrayList<>(EARLY_BAD);
            Collections.sort(wrong);
            System.out.println("Wrong shares (x-ids): " + wrong);
            System.out.println("Accuracy (shares consistent) : " +
                    (int)Math.round(100.0 * (n - wrong.size()) / Math.max(1, n)) + "/100");
            System.out.println("Consensus confidence (k-subset agreement) : " +
                    (int)Math.round(100.0 * modalCount / Math.max(1, totalCombos)) + "/100");
            return;
        }

        // 3) Check each share against the reconstructed polynomial
        Set<Integer> wrongSet = new TreeSet<>();
        for (Share s : shares) {
            Rational pred = evalPoly(goodCoeffs, BigInteger.valueOf(s.x));
            if (!pred.isInteger() || !pred.asInteger().equals(s.y)) {
                wrongSet.add(s.x);
            }
        }
        wrongSet.addAll(EARLY_BAD);

        // 4) Scores
        int wrongCount = wrongSet.size();
        int accuracy100 = (int)Math.round(100.0 * (n - wrongCount) / Math.max(1, n));
        int consensus100 = (int)Math.round(100.0 * modalCount / Math.max(1, totalCombos));

        // 5) Output
        System.out.println("Secret: " + secret);
        if (wrongSet.isEmpty()) {
            System.out.println("Wrong shares: (none detected)");
        } else {
            System.out.println("Wrong shares (x-ids): " + new ArrayList<>(wrongSet));
        }
        System.out.println("Accuracy (shares consistent) : " + accuracy100 + "/100");
        System.out.println("Consensus confidence (k-subset agreement) : " + consensus100 + "/100");
    }

    // ===== Vandermonde system: build A and solve A*c = y with Gaussian Elimination =====
    // Given k points (x_i, y_i), degree = k-1, A[i][j] = x_i^j, unknowns c[j], y[i] known.
    static List<Rational> solveVandermonde(List<Share> pts) {
        final int k = pts.size();
        Rational[][] aug = new Rational[k][k + 1]; // augmented matrix [A | y]
        for (int i = 0; i < k; i++) {
            BigInteger xi = BigInteger.valueOf(pts.get(i).x);
            // Precompute powers of xi: xi^0 .. xi^{k-1}
            BigInteger pow = BigInteger.ONE;
            for (int j = 0; j < k; j++) {
                aug[i][j] = Rational.of(pow);
                pow = pow.multiply(xi);
            }
            aug[i][k] = new Rational(pts.get(i).y, BigInteger.ONE); // RHS y_i
        }
        // Gaussian elimination with partial pivoting (exact rationals)
        gaussEliminate(aug);
        // Back substitution
        Rational[] c = new Rational[k];
        for (int i = k - 1; i >= 0; i--) {
            Rational sum = Rational.of(BigInteger.ZERO);
            for (int j = i + 1; j < k; j++) {
                sum = sum.add(aug[i][j].mul(c[j]));
            }
            Rational rhs = aug[i][k].sub(sum);
            if (aug[i][i].isZero()) throw new RuntimeException("Singular matrix");
            c[i] = rhs.div(aug[i][i]);
        }
        return Arrays.asList(c);
    }

    static void gaussEliminate(Rational[][] a) {
        int n = a.length;
        int m = a[0].length; // k+1
        for (int col = 0, row = 0; col < n && row < n; col++, row++) {
            // Pivot: find non-zero in [row..n-1] at column col
            int piv = row;
            while (piv < n && a[piv][col].isZero()) piv++;
            if (piv == n) throw new RuntimeException("Singular matrix (no pivot)");
            if (piv != row) {
                Rational[] tmp = a[piv]; a[piv] = a[row]; a[row] = tmp;
            }
            // Normalize pivot row to make pivot = 1
            Rational pivVal = a[row][col];
            for (int j = col; j < m; j++) a[row][j] = a[row][j].div(pivVal);
            // Eliminate other rows
            for (int r = 0; r < n; r++) {
                if (r == row) continue;
                Rational factor = a[r][col];
                if (factor.isZero()) continue;
                for (int j = col; j < m; j++) {
                    a[r][j] = a[r][j].sub(factor.mul(a[row][j]));
                }
            }
        }
    }

    // Evaluate polynomial c0 + c1 x + ... + c_{k-1} x^{k-1} at X (exactly)
    static Rational evalPoly(List<Rational> coeffs, BigInteger X) {
        Rational sum = Rational.of(BigInteger.ZERO);
        Rational xpow = Rational.of(BigInteger.ONE);
        Rational rx = new Rational(X, BigInteger.ONE);
        for (Rational c : coeffs) {
            sum = sum.add(c.mul(xpow));
            xpow = xpow.mul(rx);
        }
        return sum;
    }

    // ---- Combinatorics ----
    static List<int[]> combinations(int n, int k) {
        List<int[]> res = new ArrayList<>();
        if (k < 0 || k > n) return res;
        int[] comb = new int[k];
        for (int i = 0; i < k; i++) comb[i] = i;
        while (true) {
            res.add(comb.clone());
            int i = k - 1;
            while (i >= 0 && comb[i] == i + n - k) i--;
            if (i < 0) break;
            comb[i]++;
            for (int j = i + 1; j < k; j++) comb[j] = comb[j - 1] + 1;
        }
        return res;
    }
    static List<Share> subset(List<Share> all, int[] idx) {
        List<Share> out = new ArrayList<>(idx.length);
        for (int v : idx) out.add(all.get(v));
        return out;
    }

    // ---- JSON parsing ----
    static String readAll() throws IOException {
        try (BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
            StringBuilder sb = new StringBuilder();
            String line; while ((line = br.readLine()) != null) sb.append(line).append('\n');
            return sb.toString();
        }
    }
    static int extractInt(String s, String regex) {
        Matcher m = Pattern.compile(regex, Pattern.DOTALL).matcher(s);
        if (!m.find()) throw new RuntimeException("Missing int for " + regex);
        return Integer.parseInt(m.group(1));
    }
    static boolean digitsValidForBase(String s, int base) {
        for (int i = 0; i < s.length(); i++) {
            char ch = Character.toLowerCase(s.charAt(i));
            int v = Character.digit(ch, base);
            if (v == -1) return false;
        }
        return true;
    }
    static List<Share> parseShares(String json) {
        Pattern p = Pattern.compile(
            "\"(\\d+)\"\\s*:\\s*\\{\\s*\"base\"\\s*:\\s*\"(\\d+)\"\\s*,\\s*\"value\"\\s*:\\s*\"([0-9a-zA-Z]+)\"\\s*\\}",
            Pattern.DOTALL
        );
        Matcher m = p.matcher(json);
        List<Share> out = new ArrayList<>();
        while (m.find()) {
            Share s = new Share();
            s.x = Integer.parseInt(m.group(1));
            s.base = Integer.parseInt(m.group(2));
            s.raw = m.group(3);

            if (s.base < Character.MIN_RADIX || s.base > Character.MAX_RADIX) {
                EARLY_BAD.add(s.x);
                continue;
            }
            if (!digitsValidForBase(s.raw, s.base)) {
                EARLY_BAD.add(s.x);
                continue;
            }
            s.y = new BigInteger(s.raw.toLowerCase(), s.base);
            out.add(s);
        }
        out.sort(Comparator.comparingInt(a -> a.x));
        return out;
    }
}
