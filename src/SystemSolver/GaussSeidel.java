package SystemSolver;

public class GaussSeidel  {
	public static double[] gSediel(double[][] a,double[] b,double[] intialGuess, double relativeError, int iterLimit, int SF){
        // intialize the solution with zeros
        double[] lastValues = new double[b.length];
        lastValues = intialGuess.clone();
        // array to hold the solution
        double[] x = intialGuess.clone();
       // the loop with k is
        for (int k = 0; k < iterLimit; k++) {
            for (int i = 0; i < x.length; i++) {
                double sum = 0;
                for (int j = 0; j < x.length; j++) {
                    if (i != j) {
                        sum =sum + a[i][j]*x[j];
                    }
                }
                x[i] = (b[i]-sum)/a[i][i]; 
                x[i] = roundToSignificantFigures(x[i], SF);
            }
            //System.out.println("iteration number is"+k+" sol "+Arrays.toString(x));
            if (isConverged(x, lastValues, relativeError)) {
                break;
            }
            lastValues = x.clone(); 
        }  
        return x;
    }
	private static boolean isConverged(double[] xCurr, double[] xPrev, double relativeError){
        for (int i = 0; i < xPrev.length; i++) {
            double error = (xCurr[i]-xPrev[i])*100/xCurr[i];
            if (error > relativeError) {
                return false;
            }
        }
        return true;
    }
	private static double roundToSignificantFigures(double value, int significantFigures) {
        if (value == 0) return 0; // Avoid issues with log(0)

        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
        return Math.round(value * magnitude) / magnitude;
    }
}
