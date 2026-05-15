package matrix.solver;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;
import matrix.utils.MatrixResult;
import matrix.utils.ProjectMatrixUtils;

public class Jacobi implements Solver {

    public Jacobi() {
    }

    @Override
    public MatrixResult solve(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol) {
        int n = matrix.numRows;
        double[] b = vector.getDDRM().data;
        double[] x = new double[n];
        double[] nextX = new double[n];
        double validationError = 0.0;
        double relativeError = 0.0;

        double bNorm = 0.0;
        for (double v : b)
            bNorm += v * v;
        bNorm = Math.sqrt(bNorm);

        if (bNorm == 0) return new MatrixResult("Jacobi", relativeError, 0, 0.0, true, validationError, new SimpleMatrix(n, 1, true, x));

        double[] invDiag = ProjectMatrixUtils.inverseDiagonal(matrix);
        long startTime = System.nanoTime();

        for (int k = 0; k < MAX_ITER; k++) {
            double residualNormSq = 0.0;

            for (int i = 0; i < n; i++) {
                double Ax_i = 0.0;
                
                for (int j = 0; j < n; j++) {
                    double val = matrix.get(i, j);
                    if (val != 0) {
                        Ax_i += val * x[j];
                    }
                }

                double r_i = b[i] - Ax_i;
                residualNormSq += r_i * r_i;

                nextX[i] = x[i] + (r_i * invDiag[i]);
            }

            double[] temp = x;
            x = nextX;
            nextX = temp;

            relativeError = Math.sqrt(residualNormSq) / bNorm;

            if (relativeError < tol) {
                double executionTime = (System.nanoTime() - startTime) / 1e6;
                validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
                return new MatrixResult("Jacobi", relativeError, k + 1, executionTime, true, validationError, new SimpleMatrix(n, 1, true, x));
            }
        }
        validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
        double executionTime = (System.nanoTime() - startTime) / 1e6;
        return new MatrixResult("Jacobi", relativeError, MAX_ITER, executionTime, false, validationError, new SimpleMatrix(n, 1, true, x));
    }
}
