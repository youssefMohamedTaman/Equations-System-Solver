package SystemSolver;

import java.util.function.Function;

public class Bisection {
	
    public static double[] bisection(double a, double b, double error,int iterations,int sf,Function<Double, Double> f){
    	if(a>b) {
    		double temp = a;
    		a = b;
    		b = temp;
    	}
        int cnt = 0;
        double fa = roundToSignificantFigures(f.apply(a),sf);
        while (Math.abs(b - a) / 2 > error && cnt < iterations) {
            double r = roundToSignificantFigures((a + b) / 2 ,sf);
            double fr = roundToSignificantFigures(f.apply(r),sf);
            if (fr == 0) {
                break;
            } else if (fr * fa < 0) {
                b = r;
            } else {
                a = r;
            }
            //System.out.println("a= "+a+" ,b= "+b);
            cnt++;
            
        }
        System.out.println("terminated after "+cnt+" iterations");
        double root = roundToSignificantFigures((a + b) / 2 ,sf);
        double[] result = {root,cnt};
        if(Math.round(f.apply(root))!=0)
        	result[0] = Double.NaN;
        return result;
        
    }
    private static double roundToSignificantFigures(double value, int significantFigures) {
        if (value == 0) return 0; 
        
        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
        return Math.round(value * magnitude) / magnitude;
    }
	public static void main(String[] args) {
		double a = -1;
        double b = 1;
        double error = 1e-5;
        Function<Double, Double> f = x -> Math.pow(x, 2) - 2 * Math.pow(x, 1) -2;
        double[] sol = bisection(a, b, error,500,10,f);
        System.out.println("Using bisection method \nthe root = " + sol[0]+" \nReached iterations = "+sol[1]);
        
	}

}
