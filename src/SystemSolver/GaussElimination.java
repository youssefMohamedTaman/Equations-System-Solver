package SystemSolver;

public class GaussElimination  {
	public static double[] gauss(int n, double[][] a, double[] b, int sf) {
        // create the array that holds the scaling factors
        double[] s = new double[n];
        for (int i=0; i<n; i++) {
            s[i] = Math.abs(a[i][1]);
            for (int j=1; j<n; j++) {
                double abs = Math.abs(a[i][j]);
                if (abs>s[i]) {
                    s[i] = abs;
                }
            }
        }

        forwardEleminate(n, a, b, s, sf);
        double[] x = new double[n];
        backwardSubstitute(n, a, b, x, sf);
        return x;
    } 
    private static void partialPivot(int n, double[][] a, double[] b, double[] s, int r, int sf) {
        int pivotRow = r;

        double pivot = Math.abs(a[r][r]/s[r]);
        pivot = roundToSignificantFigures(pivot, sf);

        for (int i=r+1; i<n; i++) {
            double num = Math.abs(a[i][r]/s[r]);
            num = roundToSignificantFigures(num, sf);
            if (num>pivot) {
                pivot = num;
                pivotRow = i;
            }
        }

        if (pivotRow!=r) {
            for (int j=r; j<n; j++) {
                double temp = a[pivotRow][j];
                a[pivotRow][j] = a[r][j];
                a[r][j] = temp;
            }

            double temp = b[pivotRow];
            b[pivotRow] = b[r];
            b[r] = temp;

            temp = s[pivotRow];
            s[pivotRow] = s[r];
            s[r] = temp;
        }
    }

    private static void forwardEleminate(int n, double[][] a, double[] b, double[] s, int sf) {
        for (int r=0; r<n-1; r++) {
            partialPivot(n, a, b, s, r, sf);

            for (int i=r+1; i<n; i++) {
                double factor = a[i][r]/a[r][r];
                factor = roundToSignificantFigures(factor, sf);
                for (int j=r+1; j<n; j++) {
                    a[i][j] = a[i][j] - roundToSignificantFigures((factor*a[r][j]), sf);
                    a[i][j] = roundToSignificantFigures(a[i][j], sf);
                }
                b[i] = b[i] - roundToSignificantFigures((factor*b[r]), sf);
                b[i] = roundToSignificantFigures(b[i], sf);
            }
        }
    }

    private static void backwardSubstitute(int n, double[][] a, double[] b, double[] x, int sf) {
        x[n-1] = b[n-1]/a[n-1][n-1];
        x[n-1] = roundToSignificantFigures(x[n-1], sf);
        for (int i=n-2; i>=0; i--) {
            double sum = 0;
            for (int j=i+1; j<n; j++) {
                sum += roundToSignificantFigures(a[i][j]*x[j], sf);
                sum = roundToSignificantFigures(sum, sf);
            }
            x[i] = roundToSignificantFigures((b[i]-sum), sf)/a[i][i];
            x[i] = roundToSignificantFigures(x[i], sf);
        }
    }
    private static double roundToSignificantFigures(double value, int significantFigures) {
        if (value == 0) return 0; // Avoid issues with log(0)

        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
        return Math.round(value * magnitude) / magnitude;
    }
}
