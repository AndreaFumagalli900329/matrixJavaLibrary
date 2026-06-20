package matrix.utils;
import org.ejml.simple.SimpleMatrix;

public class MatrixResult {
    private String methodName;
    private double relativeError;
    private int iterations;
    private boolean converged;
    private double validationError;
    private SimpleMatrix solution;

    public MatrixResult(String methodName, double relativeError, int iterations, boolean converged, double validationError, SimpleMatrix solution) {
        this.methodName = methodName;
        this.relativeError = relativeError;
        this.iterations = iterations;
        this.converged = converged;
        this.validationError = validationError;
        this.solution = solution;
    }

    public String getMethodName() {
        return methodName;
    }

    public double getRelativeError() {
        return relativeError;
    }

    public int getIterations() {
        return iterations;
    }

    public boolean isConverged() {
        return converged;
    }

    public double getValidationError() {
        return validationError;
    }

    public SimpleMatrix getSolution() {
        return solution;
    }

    public String toString() {
        String str = String.format("Method: %s\nRelative Error: %.6e\nIterations: %d\nConverged: %b\nValidation Error: %.6e",
                methodName, relativeError, iterations, converged, validationError);
        return str;
    }
}
