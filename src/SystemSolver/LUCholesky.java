package SystemSolver;

public class LUCholesky {
    public static double[] LUch(double[][] matrix ,double[] constants ,  int sf) {
        
        double[][] lowerT = getLowerTriangular(matrix, sf);
        double[][] upper = transpose(lowerT);
        double[] y = special_forward_elimination(lowerT, constants, sf);
        double[] x = special_backward_elimination(upper, y, sf);
        return x;
    }
    private static double[][] getLowerTriangular(double[][] matrix, int sf) {
        double[][] lowerTriangle = new double[matrix.length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j <=i ; j++) {
                double sum = 0;
                if (i != j) {
                    for (int k = 0; k < j; k++) {
                        sum += (lowerTriangle[j][k] * lowerTriangle[i][k]);
                        sum = roundToSignificantFigures(sum, sf);
                    }
                    lowerTriangle[i][j] = roundToSignificantFigures((roundToSignificantFigures(matrix[i][j], sf) - sum) / lowerTriangle[j][j], sf);
                } else {
                    for (int k = 0; k < i; k++) {
                        sum += Math.pow(lowerTriangle[i][k], 2);
                        sum = roundToSignificantFigures(sum, sf);
                    }
                    lowerTriangle[i][i] = roundToSignificantFigures(Math.sqrt((roundToSignificantFigures(matrix[i][i], sf) - sum)), sf);
                }
            }
        }
        return lowerTriangle;
    }
    private static double[][] transpose(double[][] matrix) {
        double[][] transpose = new double[matrix.length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                transpose[i][j] = matrix[j][i];
            }
        }
        return transpose;
    }
    private static double[] special_forward_elimination(double[][] matrix, double[] constants, int sf) {
        for (int i = 0; i < matrix.length; i++) {
            double diagonalD = matrix[i][i];
            for (int k = i + 1; k < matrix.length; k++) {
                double multiplier = -roundToSignificantFigures(matrix[k][i] / diagonalD, sf);
                constants[k] = roundToSignificantFigures(constants[i] * multiplier, sf) + roundToSignificantFigures(constants[k], sf);
            }
            constants[i] = roundToSignificantFigures(roundToSignificantFigures(constants[i],sf) / diagonalD,sf);
        }
        return constants;
    }
    private static double[] special_backward_elimination(double[][] matrix, double[] constants, int sf) {
        for (int i = matrix.length-1; i >=0 ; i--) {
            double diagonalD = matrix[i][i];
            for (int k = i - 1; k >= 0 ; k--) {
                double multiplier = -roundToSignificantFigures(matrix[k][i] / diagonalD, sf);
                constants[k] = roundToSignificantFigures(constants[i] * multiplier, sf) + roundToSignificantFigures(constants[k], sf);
            }
            constants[i] = roundToSignificantFigures(roundToSignificantFigures(constants[i],sf) / diagonalD,sf);
        }
        return constants;
    }
    private static double roundToSignificantFigures(double value, int significantFigures) {
        if (value == 0) return 0; // Avoid issues with log(0)

        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
        return Math.round(value * magnitude) / magnitude;
    }
}
