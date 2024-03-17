package SystemSolver;
import java.util.function.Function;

public class ModifiedNewtonRaphson {
    private static Function<Double,Double> derive(Function<Double,Double> f) {
        final double dx = 0.000001;
        return (x) -> (f.apply(x + dx) - f.apply(x)) / dx;
    }

    private static double roundToSignificantFigures(double value, int significantFigures) {
        if (value == 0) return 0; 
        
        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
        return Math.round(value * magnitude) / magnitude;
    }

    public static double[] modifiedNewtonRaphson(double a, double error,int iterations,int sf,Function<Double, Double> f) {
        double x = a;
        double err = 1;
        int iters = 0;

        Function<Double,Double> firstDerivFunction = derive(f);
        Function<Double,Double> secondDerivFunction = derive(firstDerivFunction);

        while (err>error && iters<iterations) {
            double firstDeriv = firstDerivFunction.apply(x);
            double secondDeriv = secondDerivFunction.apply(x);
           
            double oldX = x;
            double numerator = roundToSignificantFigures((f.apply(oldX)*firstDeriv), sf);
            double term1 = roundToSignificantFigures((firstDeriv*firstDeriv), sf);
            double term2 = roundToSignificantFigures((f.apply(oldX)*secondDeriv), sf);
            double denominator  = roundToSignificantFigures((term1-term2), sf);
            if(Math.abs(denominator)<=1e-15) {
            	double[] result = {Double.POSITIVE_INFINITY,iters};
            	return result;
            }
            double fraction = roundToSignificantFigures((numerator/denominator), sf);
            x = roundToSignificantFigures((oldX - fraction), sf);
            if(Math.abs(f.apply(x)) <= 1e-15){
            	break;
            }
            err = Math.abs((x-oldX)/x);
            iters++;
        }
        
        double[] result = {x,iters};
        if(Math.round(f.apply(x))!=0)
        	result[0] = Double.NaN;
        return result;
    }
    public static double[] modifiedNewtonRaphson2(double a,int m, double error,int iterations,int sf,Function<Double, Double> f) {
        double x = a;
        double err = 1;
        int iters = 0;

        Function<Double,Double> firstDerivFunction = derive(f);
        while (err>error && iters<iterations) {
            double firstDeriv = firstDerivFunction.apply(x);
            
            if(Math.abs(firstDeriv)<=1e-15) {
            	double[] result = {Double.POSITIVE_INFINITY,iters};
            	return result;
            }
            double oldX = x;
            double numerator = roundToSignificantFigures((f.apply(oldX)), sf);
            double denominator  = roundToSignificantFigures(firstDeriv, sf);
            double fraction = roundToSignificantFigures((numerator/denominator), sf);
            x = roundToSignificantFigures((oldX - m*fraction), sf);
            if(Math.abs(f.apply(x)) <=1e-15) {
            	break;
            }
            err = Math.abs((x-oldX)/x);
            iters++;
            
        }
        double[] result = {x,iters};
        if(Math.round(f.apply(x))!=0)
        	result[0] = Double.NaN;
        return result;
    }
    public static void main(String[] arg) {
    	
    	Function<Double,Double> f = x -> Math.pow(x, 5) -11*Math.pow(x, 4) +46*Math.pow(x, 3) -90*Math.pow(x, 2) +81*Math.pow(x, 1) -27;
    	//Function<Double, Double> f = x -> Math.pow(x, 2) - 2 * Math.pow(x, 1) +5;
    	double[] result = modifiedNewtonRaphson(10,.00001,100000,7,f);
    	//double[] result = modifiedNewtonRaphson2(5,3,.000001,500,7,f);
    	System.out.println("Root = "+result[0]+" \nReached itrations = "+result[1]);
    }
}

