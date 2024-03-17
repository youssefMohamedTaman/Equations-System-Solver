package SystemSolver;

public class GuassJordanElimination {
	public static double[] gaussJordan(int n, double[][] a, double[] b, int precision) {
        fullEliminate(n, a, b, precision);
        return b;
    }
    private static void partialPivot2(int n, double[][] a, double[] b, int r) {
        int pivotRow = r;

        double pivot = a[r][r];
        for (int i=r+1; i<n; i++) {
            double num = Math.abs(a[i][r]);
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
        }
    }
    private static void fullEliminate(int n, double[][] a, double[] b, int sf) {
        for (int r=0; r<n; r++) {
            partialPivot2(n, a, b, r);

            //scaling the pivot row
            for (int j=r+1; j<n; j++) {
                a[r][j] = a[r][j] / a[r][r];
                a[r][j] = roundToSignificantFigures(a[r][j], sf);
            }
            b[r] = b[r]/a[r][r];
            b[r] = roundToSignificantFigures(b[r], sf);

            for (int i=0; i<n; i++) {
                if (i==r) continue;
                double factor = a[i][r];
                for (int j=r+1; j<n; j++) {
                    a[i][j] = a[i][j] - roundToSignificantFigures((factor*a[r][j]), sf);
                    a[i][j] = roundToSignificantFigures(a[i][j], sf);
                }
                b[i] = b[i] - roundToSignificantFigures((factor*b[r]), sf);
                b[i] = roundToSignificantFigures(b[i], sf);
            }
        }
    }
    private static double roundToSignificantFigures(double value, int significantFigures) {
        if (value == 0) return 0; // Avoid issues with log(0)

        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
        return Math.round(value * magnitude) / magnitude;
    }
}
