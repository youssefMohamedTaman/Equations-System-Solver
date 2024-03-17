package SystemSolver;
import java.util.function.Function;

public class secant {
	public static double[] secant(Function<Double, Double> g, double firstGuess, double secondGuess, int maxIterations , int sf , double tolerance) {
        double x0 = precision(firstGuess, sf);
        double x1 = precision(secondGuess, sf);
        double x2 = x1 - g.apply(x1)*(x1-x0)/(g.apply(x1)-g.apply(x0));
        double approximateRelativeError = 1;
        double iters = 0;
        int flag = 0;
        while(approximateRelativeError > tolerance && iters <= maxIterations) {
        	if(g.apply(x1)==g.apply(x0)) {
        		double[] result = {Double.POSITIVE_INFINITY,-1};
        		return result;
        	}
        	if(flag == 100) {
        		double[] result = {Double.NaN,-1};
        		return result;
        	}
            x2 = x1 - g.apply(x1)*(x1-x0)/(g.apply(x1)-g.apply(x0));
            x2 = precision(x2, sf);
            if(Math.abs(g.apply(x2))<= 1e-15) {
            	break;
            }
            approximateRelativeError = calcApproximateError(x2, x1);
            if(approximateRelativeError>1)
            	flag++;
            else 
            	flag = 0;
            x0 = x1;
            x1 = x2;
            iters++;
           
        }
        double[] result = {x2,iters};
        if(Math.round(g.apply(x2))!=0)
        	result[0] = Double.NaN;
        return result;
    }
	
	public static double calcApproximateError(double current, double prev) {
        return Math.abs((current - prev) / current) * 100;
    }

    public static boolean isTolerant(double approximateError, double tolerance) {
        return approximateError < tolerance;
    }
    public static double precision(double value, int digits) {
        int counter = 0;
        double x = value;
        while (x >= 1) {
            counter++;
            x /= 10;
        }
        int remain = digits - counter;
        x = value;
        x = value * Math.pow(10, remain);
        x = Math.round(x);
        x = x / Math.pow(10, remain);
        return (x);
    }

	public static void main(String[] args) {
		Function<Double, Double> f = x -> Math.pow(x, 2) - 2 * Math.pow(x, 1) +1;
		double[] sol = secant(f, 5, 11.5, 100, 10, .0000000000001);
			System.out.print(sol[0]+"\nReached num of iterations "+sol[1]);
	}

}
