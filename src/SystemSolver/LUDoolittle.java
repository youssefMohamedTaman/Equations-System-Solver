package SystemSolver;

public class LUDoolittle{
	  //LU Decomposition Doolittle Form
    public static double[] solveLUdo(double[][] coefficients, double[] constants, int significantFigures) {
        int n = coefficients.length;

        // Apply partial pivoting
        for (int i = 0; i < n; i++) {
            int maxIndex = i;
            for (int j = i + 1; j < n; j++) {
                if (Math.abs(coefficients[j][i]) > Math.abs(coefficients[maxIndex][i])) {
                    maxIndex = j;
                }
            }

            // Swap rows in coefficients and constants
            double[] tempRow = coefficients[i];
            coefficients[i] = coefficients[maxIndex];
            coefficients[maxIndex] = tempRow;

            double tempConstant = constants[i];
            constants[i] = constants[maxIndex];
            constants[maxIndex] = tempConstant;
        }

        // Perform LU decomposition (Doolittle's method)
        double[][] lower = new double[n][n];
        double[][] upper = new double[n][n];

        for (int i = 0; i < n; i++) {
            // Upper triangular matrix
            for (int k = i; k < n; k++) {
                double sum = 0.0;
                for (int j = 0; j < i; j++) {
                    sum += roundToSignificantFigures(lower[i][j] * upper[j][k], significantFigures);
                }
                upper[i][k] = roundToSignificantFigures(coefficients[i][k] - sum, significantFigures);
            }

            // Lower triangular matrix
            lower[i][i] = 1; // Diagonal elements are 1
            for (int k = i + 1; k < n; k++) {
                double sum = 0.0;
                for (int j = 0; j < i; j++) {
                    sum += lower[k][j] * upper[j][i];
                }
                lower[k][i] = roundToSignificantFigures((coefficients[k][i] - sum) / upper[i][i], significantFigures);
            }
        }

        // Solve Ly = b using forward substitution
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += roundToSignificantFigures(lower[i][j] * y[j], significantFigures);
            }
            y[i] = roundToSignificantFigures(constants[i] - sum, significantFigures);
        }

        // Solve Ux = y using backward substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += roundToSignificantFigures(upper[i][j] * x[j], significantFigures);
            }
            x[i] = roundToSignificantFigures((y[i] - sum) / upper[i][i], significantFigures);
        }
        for (int i = 0; i < n; i++) {
            x[i] = roundToSignificantFigures(x[i], significantFigures);
        }
		return x;
    }
    private static double roundToSignificantFigures(double value, int significantFigures) {
        if (value == 0) return 0; // Avoid issues with log(0)

        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
        return Math.round(value * magnitude) / magnitude;
    }
}
