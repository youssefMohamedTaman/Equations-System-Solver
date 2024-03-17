package SystemSolver;

import java.util.function.Function;

public class FalsePosition {
	    public static double[] falsePosition(double a, double b, double error,int iterations,int sf,Function<Double, Double> f){
	    	if(a>b) {
        		double temp = a;
        		a = b;
        		b = temp;
        	}
	        double fa = roundToSignificantFigures(f.apply(a),sf);
	        double fb = roundToSignificantFigures(f.apply(b),sf);
	        double cnt = 0;
	        int flag = 0;
	        double prev=0;
	        double approximateError =1;
	        while (approximateError > error && cnt < iterations) {
		         fa = roundToSignificantFigures(f.apply(a),sf);
		         fb = roundToSignificantFigures(f.apply(b),sf);
	        	if(fb==fa) {
	        		double result[] = {Double.POSITIVE_INFINITY,-1};
	        		return result;
	        	}
	            double r = roundToSignificantFigures((a * fb - b * fa),sf) / roundToSignificantFigures((fb - fa),sf);
	            r = roundToSignificantFigures(r,sf);
	            double fr = roundToSignificantFigures(f.apply(r),sf);
	            approximateError = roundToSignificantFigures((r-prev)/r,sf)*100;
	            if(Math.round(r*10000) == Math.round(prev*10000))
	            	flag++;
	            else 
	            	flag=0;
	            prev = r;
	            if(flag==5) {
	            	if (fr * fa < 0) {
	                    b = r+.01;
	                    
	                } else {
	                    a = r+.01;
	                }
	            	System.out.println("shifted");
	            	flag=0;
	            	continue;	
	            }
	            
	            if (fr == 0) {
	                break;
	            } else if (fr * fa < 0) {
	                b = r;
	            } else {
	                a = r;
	            }
	            cnt++;
	        }
	        //System.out.println("terminated after "+cnt+" iterations");
	        double xr = roundToSignificantFigures((a * fb - b * fa),sf)/roundToSignificantFigures((fb - fa),sf);
	        xr = roundToSignificantFigures(xr,sf);
	        //System.out.println("prev "+prev+" root = "+xr+" f(root) = "+f.apply(xr));
	        double []result = {xr,cnt};
	        if(Math.round(f.apply(xr))!=0)
	        	result[0] = Double.NaN;
	        return result;
	        
	    }
	    private static double roundToSignificantFigures(double value, int significantFigures) {
	        if (value == 0) return 0; // Avoid issues with log(0)

	        double magnitude = Math.pow(10, significantFigures - 1 - (int) Math.floor(Math.log10(Math.abs(value))));
	        return Math.round(value * magnitude) / magnitude;
	    }
	public static void main(String[] args) {
        double a = 0;
        double b = 3;
        double error = .00000000001;
		Function<Double, Double> f = z -> Math.pow(z, 4)+3*z-4;
        //Function<Double, Double> f = z -> Math.pow(z, 3)-.165*Math.pow(z, 2)+0.0003993;
        //Function<Double, Double> f = x -> Math.pow(x, 2) - 2 *x +5;
        //Function<Double, Double> f = x -> Math.pow(x, 3) - 1;
		double[] sol = falsePosition(a, b, error,200,10,f);
        System.out.println("Using false position method \nthe root = " + sol[0] +" \nReached iterations = " +sol[1]);

	}
}
