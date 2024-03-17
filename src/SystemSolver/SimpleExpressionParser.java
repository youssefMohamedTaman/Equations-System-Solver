package SystemSolver;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

public class SimpleExpressionParser {

    public static double evaluateExpression(double[] coefficients , double x) {
        // Parse the coefficients of the polynomial
        //double[] coefficients = parseCoefficients(expression);

        // Create a polynomial function
        PolynomialFunction polynomialFunction = new PolynomialFunction(coefficients);

        // Evaluate the polynomial function at the specified value of 'x'
        return polynomialFunction.value(x);
    }

//    private static double[] parseCoefficients(String expression) {
//        // Implement your custom logic to parse coefficients from the expression
//        // For simplicity, let's assume coefficients are space-separated numbers
//        String[] coefficientStrings = expression.split(" ");
//        double[] coefficients = new double[coefficientStrings.length];
//
//        for (int i = coefficientStrings.length-1; i>=0; i--) {
//            coefficients[i] = Double.parseDouble(coefficientStrings[i]);
//        }
//
//        return coefficients;
//    }

    public static void main(String[] args) {
        // Example usage
        //String expression = "3.0 6.5 4.0";  // Coefficients of the polynomial 4x^2 + 6.5x + 3
        double [] expression ={3.0, 6.5, 4};
        double x = 5;

        //double result = evaluateExpression(expression, 5);
        //System.out.println("The result of the expression at x = " + x + " is: " + result);
        System.out.println(Math.sin(.5235987756));
    }
}








