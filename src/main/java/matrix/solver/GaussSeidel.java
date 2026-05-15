package matrix.solver;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;
import matrix.utils.MatrixResult;
import matrix.utils.ProjectMatrixUtils;

public class GaussSeidel implements Solver {

    public GaussSeidel() {}

    @Override
    public MatrixResult solve(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol) {
        int n = matrix.numRows;

        double[] b = vector.getDDRM().data;
        double[] x = new double[n];
        double validationError = 0.0;
        double relativeError = 0.0;

        double bNorm = 0.0;
        for (double val : b) bNorm += val * val;
        bNorm = Math.sqrt(bNorm);

        if (bNorm == 0) {
            return new MatrixResult("Gauss-Seidel", relativeError, 0, 0.0, true, validationError, new SimpleMatrix(n, 1, true, x));
        }

        double[] invDiag = ProjectMatrixUtils.inverseDiagonal(matrix);
        long startTime = System.nanoTime();

        for (int k = 0; k < MAX_ITER; k++) {
            // Aggiorno x[i] e lo uso subito
            for (int i = 0; i < n; i++) {
                double Ax_i = 0.0;
                
                // Prodotto riga i-esima
                // Usiamo x[j] che contiene sia valori "vecchi" (se j > i) sia valori "nuovi" (se j < i)
                for (int j = 0; j < n; j++) {
                    double val = matrix.get(i, j);
                    if (val != 0) {
                        Ax_i += val * x[j];
                    }
                }

                // Correzione in-place: x_new = x_old + (b - Ax_old) / A_ii
                x[i] = x[i] + (b[i] - Ax_i) * invDiag[i];
            }

            // Calcolo Residuo Manuale (dopo aver completato l'aggiornamento di tutto il vettore x)
            double residualNormSq = 0.0;
            for (int i = 0; i < n; i++) {
                double currentAx_i = 0.0;
                for (int j = 0; j < n; j++) {
                    double val = matrix.get(i, j);
                    if (val != 0) currentAx_i += val * x[j];
                }
                double r = b[i] - currentAx_i;
                residualNormSq += r * r;
            }

            relativeError = Math.sqrt(residualNormSq) / bNorm;

            if (relativeError < tol) {
                double executionTime = (System.nanoTime() - startTime) / 1e6;
                validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
                return new MatrixResult("Gauss-Seidel", relativeError, k + 1, executionTime, true, validationError, new SimpleMatrix(n, 1, true, x));
            }
        }
        
        validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
        double executionTime = (System.nanoTime() - startTime) / 1e6;
        return new MatrixResult("Gauss-Seidel", relativeError, MAX_ITER, executionTime, false, validationError, new SimpleMatrix(n, 1, true, x));
    }
}