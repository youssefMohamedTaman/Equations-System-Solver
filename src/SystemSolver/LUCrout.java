package SystemSolver;

public class LUCrout {
	public static double[] solveLUcr(double[][] coefficients, double[] constants, int significantFigures) {
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

        // Perform LU decomposition (Crout's method)
        double[][] lower = new double[n][n];
        double[][] upper = new double[n][n];

        for(int i = 0; i < n; i++) {
        	upper[i][i] = 1;
        }
        
       for(int j = 0; j < n; j++) {
        	for(int i = j; i < n ; i++) {
        		double sum = 0;
        		for(int k = 0; k < j ; k++) {
        			sum += roundToSignificantFigures(lower[i][k] * upper[k][j],significantFigures);
        		}
        		lower[i][j] = roundToSignificantFigures(coefficients[i][j] - sum,significantFigures);
        	}
        	
        	for (int i = j; i < n; i++) {
    			double sum = 0;
    			for(int k = 0; k < j; k++) {
    				sum += roundToSignificantFigures(lower[j][k] * upper[k][i],significantFigures);
    			}
    			if (lower[j][j] == 0) {
    				System.out.println("det(L) close to 0!\n Can't divide by 0...\n");
    				System.exit(0);
    			}
    			upper[j][i] = roundToSignificantFigures((coefficients[j][i] - sum) / lower[j][j],significantFigures);
    		}
        }
        

        // Solve Ly = b using forward substitution
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += roundToSignificantFigures(lower[i][j] * y[j],significantFigures);
            }
            y[i] = roundToSignificantFigures((constants[i] - sum) / lower[i][i],significantFigures);
        }

        // Solve Ux = y using backward substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += roundToSignificantFigures(upper[i][j] * x[j],significantFigures);
            }
            x[i] = roundToSignificantFigures((y[i] - sum) / upper[i][i],significantFigures);
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
