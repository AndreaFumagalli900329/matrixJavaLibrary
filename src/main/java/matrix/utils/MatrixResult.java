package matrix.utils;

import org.ejml.simple.SimpleMatrix;

public class MatrixResult {
    private String methodName;
    private double relativeError;
    private int iterations;
    private boolean converged;
    private double validationError;
    private SimpleMatrix solution;

    public MatrixResult(String methodName, double relativeError, int iterations, boolean converged,
            double validationError, SimpleMatrix solution) {
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
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("Method: %s\n", methodName));
        sb.append(String.format("Relative Error: %.6e\n", relativeError));
        sb.append(String.format("Iterations: %d\n", iterations));
        sb.append(String.format("Converged: %b\n", converged));
        sb.append(String.format("Validation Error: %.6e\n", validationError));

        sb.append("Solution Vector:\n");
        sb.append(ProjectMatrixUtils.vectorToString(solution));

        return sb.toString();
    }
}
