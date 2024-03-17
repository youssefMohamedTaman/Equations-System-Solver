package SystemSolver;

import java.util.function.Function;

public class NewtonRaphason {
    private static Function<Double,Double> derive(Function<Double,Double> f) {
        final double dx = 0.000001;
        Function<Double ,Double> ff = x -> (f.apply(x + dx) - f.apply(x)) / dx;
        return ff;
    }

    private static double roundToSignificantFigures(double value, int significantFigures) {
        if (value == 0) return 0; 
        
        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
        return Math.round(value * magnitude) / magnitude;
    }

    public static double[] NewtonRaphson(double a, double error,int iterations,int sf,Function<Double, Double> f) {
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
            x = roundToSignificantFigures((oldX - fraction), sf);
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
    	
    	//Function<Double,Double> f = x -> Math.exp(-x)-x;
    	Function<Double, Double> f = x -> Math.pow(x, 2) - 2 * Math.pow(x, 1) +1;
    	double[] result = NewtonRaphson(0,.000001,500,7,f);
    	System.out.println("Root = "+result[0]+" \nReached itrations = "+result[1]);
    }


}
