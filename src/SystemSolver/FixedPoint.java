package SystemSolver;
import java.util.function.Function;

public class FixedPoint {

    public static double[] fixedPoint(Function<Double, Double> g, double initialInput, int maxIterations , int sf , double tolerance) {
        double xr = precision(initialInput, sf);
        double prevXr = xr;
        double approximateRelativeError=1;
        double iters = 0;
        int flag = 0;
        while(iters < maxIterations && approximateRelativeError > tolerance) {
        	if(flag==5 || Math.abs(xr-prevXr)>1e5) {
        		double[] result = {Double.NaN,-1};
        		return result;
        	}
        	if(Math.abs(g.apply(xr))<=1e-15) {
        		break;
        	}
            prevXr = xr;
            xr = precision(g.apply(xr), sf);
            approximateRelativeError = calcApproximateError(xr, prevXr);
            if(Math.abs(xr-prevXr)>10)
            	flag++;
            else 
            	flag = 0;
            iters++; 
        }
        double[] result = {xr,iters};
        return result;
    }

    public static double calcApproximateError(double current, double prev) {
        return Math.abs((current - prev) / current) * 100;
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
    public static void main(String arg[]) {
    	 //Function<Double, Double> f = x -> Math.exp(-x);
    	 //Function<Double, Double> f = x -> Math.pow(x, 2) - 2 * Math.pow(x, 1)+1 ;
    	 //Function<Double, Double> f = x -> Math.pow(x, 4)+ 3 * Math.pow(x, 1) - 4 ;
    	 //Function<Double, Double> f = x -> Math.sqrt(2*x+3);
    	 Function<Double, Double> f = x -> 3/(x-2);
         double[] result = fixedPoint(f,4,200,7,0.00001); 
         System.out.println("Root: " + result[0]+" \nReached iterations "+result[1]);
    }
}
