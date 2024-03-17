package SystemSolver;

public abstract class Solver {

	public static boolean check(double []solution) {
		int n = solution.length;
		for(int i=0;i<n;i++) {
			if(Double.isNaN(solution[i]))
				return false;
		}
		return true;
	} 
    public static int checkSolutionType(double[][] coefficients, double[] constants) {
        int numRows = coefficients.length;
        int numCols = coefficients[0].length;

        // Create an augmented matrix [A | B]
        double[][] augmentedMatrix = new double[numRows][numCols + 1];
        for (int i = 0; i < numRows; i++) {
            System.arraycopy(coefficients[i], 0, augmentedMatrix[i], 0, numCols);
            augmentedMatrix[i][numCols] = constants[i];
        }

        // Apply Gaussian elimination
        for (int i = 0; i < numRows; i++) {
            // Make the diagonal element 1
            double diagonalElement = augmentedMatrix[i][i];
            if (diagonalElement == 0) {
                // Check for a row of zeros (infinite solutions) or inconsistent system
                if (!hasNonZeroRow(augmentedMatrix, i, numCols)) {
                    return 2;
                } else {
                    return 0;
                }
            }
            for (int j = 0; j <= numCols; j++) {
                augmentedMatrix[i][j] /= diagonalElement;
            }

            // Eliminate other elements in the current column
            for (int k = 0; k < numRows; k++) {
                if (k != i) {
                    double factor = augmentedMatrix[k][i];
                    for (int j = 0; j <= numCols; j++) {
                        augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                    }
                }
            }
        }
        // Check for free variables (infinite solutions)
        for (int i = 0; i < numRows; i++) {
            boolean rowHasPivot = false;
            for (int j = 0; j < numCols; j++) {
                if (augmentedMatrix[i][j] != 0) {
                    rowHasPivot = true;
                    break;
                }
            }
            if (!rowHasPivot && augmentedMatrix[i][numCols] != 0) {
                return 2;
            }
        }

        return 1;
    }

    public static boolean hasNonZeroRow(double[][] matrix, int row, int numCols) {
        for (int j = 0; j < numCols; j++) {
            if (matrix[row][j] != 0) {
                return true;
            }
        }
        return false;
    }
    public static boolean forSymmetry(double[][] a){
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a.length; j++) {
                if (i!=j) {
                    if (a[i][j] != a[j][i]) {
                        return false;
                    }
                }
            }
        }
        return true;
     }
    public static boolean forDiagonalDominance(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            double sum  = 0;
            for (int j = 0; j < a.length; j++) {
                if (i!=j) {
                    sum += Math.abs(a[i][j]);
                }
            }
            if (Math.abs(a[i][i]) < sum) {
                return false;
            }
        }
        return true;
     }
}
//  //LU Decomposition Doolittle Form
//        public static double[] solveLUdo(double[][] coefficients, double[] constants, int significantFigures) {
//            int n = coefficients.length;
//
//            // Apply partial pivoting
//            for (int i = 0; i < n; i++) {
//                int maxIndex = i;
//                for (int j = i + 1; j < n; j++) {
//                    if (Math.abs(coefficients[j][i]) > Math.abs(coefficients[maxIndex][i])) {
//                        maxIndex = j;
//                    }
//                }
//
//                // Swap rows in coefficients and constants
//                double[] tempRow = coefficients[i];
//                coefficients[i] = coefficients[maxIndex];
//                coefficients[maxIndex] = tempRow;
//
//                double tempConstant = constants[i];
//                constants[i] = constants[maxIndex];
//                constants[maxIndex] = tempConstant;
//            }
//
//            // Perform LU decomposition (Doolittle's method)
//            double[][] lower = new double[n][n];
//            double[][] upper = new double[n][n];
//
//            for (int i = 0; i < n; i++) {
//                // Upper triangular matrix
//                for (int k = i; k < n; k++) {
//                    double sum = 0.0;
//                    for (int j = 0; j < i; j++) {
//                        sum += lower[i][j] * upper[j][k];
//                    }
//                    upper[i][k] = coefficients[i][k] - sum;
//                }
//
//                // Lower triangular matrix
//                lower[i][i] = 1; // Diagonal elements are 1
//                for (int k = i + 1; k < n; k++) {
//                    double sum = 0.0;
//                    for (int j = 0; j < i; j++) {
//                        sum += lower[k][j] * upper[j][i];
//                    }
//                    lower[k][i] = (coefficients[k][i] - sum) / upper[i][i];
//                }
//            }
//
//            // Solve Ly = b using forward substitution
//            double[] y = new double[n];
//            for (int i = 0; i < n; i++) {
//                double sum = 0.0;
//                for (int j = 0; j < i; j++) {
//                    sum += lower[i][j] * y[j];
//                }
//                y[i] = constants[i] - sum;
//            }
//
//            // Solve Ux = y using backward substitution
//            double[] x = new double[n];
//            for (int i = n - 1; i >= 0; i--) {
//                double sum = 0.0;
//                for (int j = i + 1; j < n; j++) {
//                    sum += upper[i][j] * x[j];
//                }
//                x[i] = (y[i] - sum) / upper[i][i];
//            }
//            for (int i = 0; i < n; i++) {
//                x[i] = roundToSignificantFigures(x[i], significantFigures);
//            }
//    		return x;
//        }
//    public static double[] solveLUcr(double[][] coefficients, double[] constants, int significantFigures) {
//            int n = coefficients.length;
//
//            // Apply partial pivoting
//            for (int i = 0; i < n; i++) {
//                int maxIndex = i;
//                for (int j = i + 1; j < n; j++) {
//                    if (Math.abs(coefficients[j][i]) > Math.abs(coefficients[maxIndex][i])) {
//                        maxIndex = j;
//                    }
//                }
//
//                // Swap rows in coefficients and constants
//                double[] tempRow = coefficients[i];
//                coefficients[i] = coefficients[maxIndex];
//                coefficients[maxIndex] = tempRow;
//
//                double tempConstant = constants[i];
//                constants[i] = constants[maxIndex];
//                constants[maxIndex] = tempConstant;
//            }
//
//            // Perform LU decomposition (Crout's method)
//            double[][] lower = new double[n][n];
//            double[][] upper = new double[n][n];
//
//            for(int i = 0; i < n; i++) {
//            	upper[i][i] = 1;
//            }
//            
//           for(int j = 0; j < n; j++) {
//            	for(int i = j; i < n ; i++) {
//            		double sum = 0;
//            		for(int k = 0; k < j ; k++) {
//            			sum += lower[i][k] * upper[k][j];
//            		}
//            		lower[i][j] = coefficients[i][j] - sum;
//            	}
//            	
//            	for (int i = j; i < n; i++) {
//        			double sum = 0;
//        			for(int k = 0; k < j; k++) {
//        				sum = sum + lower[j][k] * upper[k][i];
//        			}
//        			if (lower[j][j] == 0) {
//        				System.out.println("det(L) close to 0!\n Can't divide by 0...\n");
//        				System.exit(0);
//        			}
//        			upper[j][i] = (coefficients[j][i] - sum) / lower[j][j];
//        		}
//            }
//            
//
//            // Solve Ly = b using forward substitution
//            double[] y = new double[n];
//            for (int i = 0; i < n; i++) {
//                double sum = 0.0;
//                for (int j = 0; j < i; j++) {
//                    sum += lower[i][j] * y[j];
//                }
//                y[i] = (constants[i] - sum) / lower[i][i];
//            }
//
//            // Solve Ux = y using backward substitution
//            double[] x = new double[n];
//            for (int i = n - 1; i >= 0; i--) {
//                double sum = 0.0;
//                for (int j = i + 1; j < n; j++) {
//                    sum += upper[i][j] * x[j];
//                }
//                x[i] = (y[i] - sum) / upper[i][i];
//            }
//            for (int i = 0; i < n; i++) {
//                x[i] = roundToSignificantFigures(x[i], significantFigures);
//            }
//    		return x;
//        }
//    public static double[] jacobi(double[][] a,double[] b,double[] intialGuess, double relativeError, int iterLimit, int SF){
//        // intialize the solution with zeros
//        double[] nextIterValues = new double[b.length];
//        nextIterValues = intialGuess.clone();
//        // array to hold the solution
//        double[] x = new double[b.length];
//       // the loop with k is
//        for (int k = 0; k < iterLimit; k++) {
//            for (int i = 0; i < x.length; i++) {
//                double sum = 0;
//                for (int j = 0; j < x.length; j++) {
//                    if (i != j) {
//                        sum =sum + a[i][j]*nextIterValues[j];
//                    }
//                }
//                x[i] = (b[i]-sum)/a[i][i]; 
//                x[i] = roundToSignificantFigures(x[i], SF);
//            }
//            //System.out.println("iteration number is"+k+" sol "+Arrays.toString(x));
//            if (isConverged(x, nextIterValues, relativeError)) {
//                break;
//            }
//             nextIterValues = x.clone();
//            
//        }
//        
//        return x;
//    }
//    public static double[] gSediel(double[][] a,double[] b,double[] intialGuess, double relativeError, int iterLimit, int SF){
//        // intialize the solution with zeros
//        double[] lastValues = new double[b.length];
//        lastValues = intialGuess.clone();
//        // array to hold the solution
//        double[] x = intialGuess.clone();
//       // the loop with k is
//        for (int k = 0; k < iterLimit; k++) {
//            for (int i = 0; i < x.length; i++) {
//                double sum = 0;
//                for (int j = 0; j < x.length; j++) {
//                    if (i != j) {
//                        sum =sum + a[i][j]*x[j];
//                    }
//                }
//                x[i] = (b[i]-sum)/a[i][i]; 
//                x[i] = roundToSignificantFigures(x[i], SF);
//            }
//            System.out.println("iteration number is"+k+" sol "+Arrays.toString(x));
//            if (isConverged(x, lastValues, relativeError)) {
//                break;
//            }
//            lastValues = x.clone(); 
//        }  
//        return x;
//    }
    
//    protected static boolean isConverged(double[] xCurr, double[] xPrev, double relativeError){
//        for (int i = 0; i < xPrev.length; i++) {
//            double error = (xCurr[i]-xPrev[i])*100/xCurr[i];
//            if (error > relativeError) {
//                return false;
//            }
//        }
//        return true;
//    }
//    public static double[] LUch(double[][] matrix ,double[] constants ,  int sf) {
//       
//        double[][] lowerT = getLowerTriangular(matrix, sf);
//        double[][] upper = transpose(lowerT);
//        double[] y = special_forward_elimination(lowerT, constants, sf);
//        double[] x = special_backward_elimination(upper, y, sf);
//        return x;
//    }
//    private static double[][] getLowerTriangular(double[][] matrix, int sf) {
//        double[][] lowerTriangle = new double[matrix.length][matrix.length];
//        for (int i = 0; i < matrix.length; i++) {
//            for (int j = 0; j <=i ; j++) {
//                double sum = 0;
//                if (i != j) {
//                    for (int k = 0; k < j; k++) {
//                        sum += (lowerTriangle[j][k] * lowerTriangle[i][k]);
//                        sum = roundToSignificantFigures(sum, sf);
//                    }
//                    lowerTriangle[i][j] = roundToSignificantFigures((roundToSignificantFigures(matrix[i][j], sf) - sum) / lowerTriangle[j][j], sf);
//                } else {
//                    for (int k = 0; k < i; k++) {
//                        sum += Math.pow(lowerTriangle[i][k], 2);
//                        sum = roundToSignificantFigures(sum, sf);
//                    }
//                    lowerTriangle[i][i] = roundToSignificantFigures(Math.sqrt((roundToSignificantFigures(matrix[i][i], sf) - sum)), sf);
//                }
//            }
//        }
//        return lowerTriangle;
//    }
//    private static double[][] transpose(double[][] matrix) {
//        double[][] transpose = new double[matrix.length][matrix.length];
//        for (int i = 0; i < matrix.length; i++) {
//            for (int j = 0; j < matrix.length; j++) {
//                transpose[i][j] = matrix[j][i];
//            }
//        }
//        return transpose;
//    }
//    private static double[] special_forward_elimination(double[][] matrix, double[] constants, int sf) {
//        for (int i = 0; i < matrix.length; i++) {
//            double diagonalD = matrix[i][i];
//            for (int k = i + 1; k < matrix.length; k++) {
//                double multiplier = -roundToSignificantFigures(matrix[k][i] / diagonalD, sf);
//                constants[k] = roundToSignificantFigures(constants[i] * multiplier, sf) + roundToSignificantFigures(constants[k], sf);
//            }
//            constants[i] = roundToSignificantFigures(roundToSignificantFigures(constants[i],sf) / diagonalD,sf);
//        }
//        return constants;
//    }
//    private static double[] special_backward_elimination(double[][] matrix, double[] constants, int sf) {
//        for (int i = matrix.length-1; i >=0 ; i--) {
//            double diagonalD = matrix[i][i];
//            for (int k = i - 1; k >= 0 ; k--) {
//                double multiplier = -roundToSignificantFigures(matrix[k][i] / diagonalD, sf);
//                constants[k] = roundToSignificantFigures(constants[i] * multiplier, sf) + roundToSignificantFigures(constants[k], sf);
//            }
//            constants[i] = roundToSignificantFigures(roundToSignificantFigures(constants[i],sf) / diagonalD,sf);
//        }
//        return constants;
//    }
//    public static double[] gauss(int n, double[][] a, double[] b, int sf) {
//        // create the array that holds the scaling factors
//        double[] s = new double[n];
//        for (int i=0; i<n; i++) {
//            s[i] = Math.abs(a[i][1]);
//            for (int j=1; j<n; j++) {
//                double abs = Math.abs(a[i][j]);
//                if (abs>s[i]) {
//                    s[i] = abs;
//                }
//            }
//        }
//
//        forwardEleminate(n, a, b, s, sf);
//        double[] x = new double[n];
//        backwardSubstitute(n, a, b, x, sf);
//        return x;
//    } 
//    private static void partialPivot(int n, double[][] a, double[] b, double[] s, int r, int sf) {
//        int pivotRow = r;
//
//        double pivot = Math.abs(a[r][r]/s[r]);
//        pivot = roundToSignificantFigures(pivot, sf);
//
//        for (int i=r+1; i<n; i++) {
//            double num = Math.abs(a[i][r]/s[r]);
//            num = roundToSignificantFigures(num, sf);
//            if (num>pivot) {
//                pivot = num;
//                pivotRow = i;
//            }
//        }
//
//        if (pivotRow!=r) {
//            for (int j=r; j<n; j++) {
//                double temp = a[pivotRow][j];
//                a[pivotRow][j] = a[r][j];
//                a[r][j] = temp;
//            }
//
//            double temp = b[pivotRow];
//            b[pivotRow] = b[r];
//            b[r] = temp;
//
//            temp = s[pivotRow];
//            s[pivotRow] = s[r];
//            s[r] = temp;
//        }
//    }
//
//    private static void forwardEleminate(int n, double[][] a, double[] b, double[] s, int sf) {
//        for (int r=0; r<n-1; r++) {
//            partialPivot(n, a, b, s, r, sf);
//
//            for (int i=r+1; i<n; i++) {
//                double factor = a[i][r]/a[r][r];
//                factor = roundToSignificantFigures(factor, sf);
//                for (int j=r+1; j<n; j++) {
//                    a[i][j] = a[i][j] - roundToSignificantFigures((factor*a[r][j]), sf);
//                    a[i][j] = roundToSignificantFigures(a[i][j], sf);
//                }
//                b[i] = b[i] - roundToSignificantFigures((factor*b[r]), sf);
//                b[i] = roundToSignificantFigures(b[i], sf);
//            }
//        }
//    }
//
//    private static void backwardSubstitute(int n, double[][] a, double[] b, double[] x, int sf) {
//        x[n-1] = b[n-1]/a[n-1][n-1];
//        x[n-1] = roundToSignificantFigures(x[n-1], sf);
//        for (int i=n-2; i>=0; i--) {
//            double sum = 0;
//            for (int j=i+1; j<n; j++) {
//                sum += roundToSignificantFigures(a[i][j]*x[j], sf);
//                sum = roundToSignificantFigures(sum, sf);
//            }
//            x[i] = roundToSignificantFigures((b[i]-sum), sf)/a[i][i];
//            x[i] = roundToSignificantFigures(x[i], sf);
//        }
//    }
//
//    // guass jordon elimination 
//    public static double[] gaussJordan(int n, double[][] a, double[] b, int precision) {
//        fullEliminate(n, a, b, precision);
//        return b;
//    }
//    private static void partialPivot2(int n, double[][] a, double[] b, int r) {
//        int pivotRow = r;
//
//        double pivot = a[r][r];
//        for (int i=r+1; i<n; i++) {
//            double num = Math.abs(a[i][r]);
//            if (num>pivot) {
//                pivot = num;
//                pivotRow = i;
//            }
//        }
//
//        if (pivotRow!=r) {
//            for (int j=r; j<n; j++) {
//                double temp = a[pivotRow][j];
//                a[pivotRow][j] = a[r][j];
//                a[r][j] = temp;
//            }
//
//            double temp = b[pivotRow];
//            b[pivotRow] = b[r];
//            b[r] = temp;
//        }
//    }
//    private static void fullEliminate(int n, double[][] a, double[] b, int sf) {
//        for (int r=0; r<n; r++) {
//            partialPivot2(n, a, b, r);
//
//            //scaling the pivot row
//            for (int j=r+1; j<n; j++) {
//                a[r][j] = a[r][j] / a[r][r];
//                a[r][j] = roundToSignificantFigures(a[r][j], sf);
//            }
//            b[r] = b[r]/a[r][r];
//            b[r] = roundToSignificantFigures(b[r], sf);
//
//            for (int i=0; i<n; i++) {
//                if (i==r) continue;
//                double factor = a[i][r];
//                for (int j=r+1; j<n; j++) {
//                    a[i][j] = a[i][j] - roundToSignificantFigures((factor*a[r][j]), sf);
//                    a[i][j] = roundToSignificantFigures(a[i][j], sf);
//                }
//                b[i] = b[i] - roundToSignificantFigures((factor*b[r]), sf);
//                b[i] = roundToSignificantFigures(b[i], sf);
//            }
//        }
//    }
//    public static boolean hasUniqueSolu(double[][] coefficients, double[] constants) {
//        int numRows = coefficients.length;
//        int numCols = coefficients[0].length;
//
//        // Create an augmented matrix [A | B]
//        double[][] augmentedMatrix = new double[numRows][numCols + 1];
//        for (int i = 0; i < numRows; i++) {
//            System.arraycopy(coefficients[i], 0, augmentedMatrix[i], 0, numCols);
//            augmentedMatrix[i][numCols] = constants[i];
//        }
//
//        // Apply Gaussian elimination
//        for (int i = 0; i < numRows; i++) {
//            // Make the diagonal element 1
//            double diagonalElement = augmentedMatrix[i][i];
//            if (diagonalElement == 0) {
//                return false;  // The system has no unique solution
//            }
//            for (int j = 0; j <= numCols; j++) {
//                augmentedMatrix[i][j] /= diagonalElement;
//            }
//
//            // Eliminate other elements in the current column
//            for (int k = 0; k < numRows; k++) {
//                if (k != i) {
//                    double factor = augmentedMatrix[k][i];
//                    for (int j = 0; j <= numCols; j++) {
//                        augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
//                    }
//                }
//            }
//        }
//        // Check for inconsistent or dependent equations
//        for (int i = 0; i < numRows; i++) {
//            boolean allZero = true;
//            for (int j = 0; j < numCols; j++) {
//                if (augmentedMatrix[i][j] != 0) {
//                    allZero = false;
//                    break;
//                }
//            }
//            if (allZero && augmentedMatrix[i][numCols] != 0) {
//                return false;  // The system is inconsistent
//            }
//        }
//
//        return true;  // The system has a unique solution
//    }
//    private static double roundToSignificantFigures(double value, int significantFigures) {
//        if (value == 0) return 0; // Avoid issues with log(0)
//
//        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
//        return Math.round(value * magnitude) / magnitude;
//    }

