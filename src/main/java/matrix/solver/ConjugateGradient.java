package matrix.solver;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;
import org.ejml.sparse.csc.CommonOps_DSCC;

import matrix.utils.MatrixResult;
import matrix.utils.ProjectMatrixUtils;

public class ConjugateGradient implements Solver {
    public ConjugateGradient() {
    }

    @Override
    public MatrixResult solve(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol) {
        int n = matrix.numRows;
        double[] b = vector.getDDRM().data;
        double[] x = new double[n];
        double num = 0.0;
        double den = 0.0;
        double validationError = 0.0;
        double relativeError = 0.0;

        double a = 0.0;
        double[] r = new double[n];
        double[] d = new double[n];
        double[] y = new double[n];

        double bNorm = 0.0;
        for (double val : b)
            bNorm += val * val;
        bNorm = Math.sqrt(bNorm);

        if (bNorm == 0) {
            return new MatrixResult("Conjugate Gradient", relativeError, 0, true, validationError,
                    new SimpleMatrix(n, 1, true, x));
        }

        r = b.clone();
        d = r.clone();

        // Wrappiamo gli array in strutture EJML senza copiare dati in memoria
        DMatrixRMaj d_mat = DMatrixRMaj.wrap(n, 1, d);
        DMatrixRMaj y_mat = DMatrixRMaj.wrap(n, 1, y);

        for (int k = 0; k < MAX_ITER; k++) {
            double residualNormSq = 0.0;
            num = 0.0;
            den = 0.0;
            // libreria per calcolare efficientemente y = A * d
            CommonOps_DSCC.mult(matrix, d_mat, y_mat);

            for (int i = 0; i < n; i++) {
                num += d[i] * r[i];
                den += y[i] * d[i];
            }
            a = num / den;

            for (int i = 0; i < n; i++) {
                x[i] = x[i] + a * d[i];
                r[i] = r[i] - a * y[i];
                residualNormSq += r[i] * r[i];
            }

            // usando la regola della simmetria posso usare y invece di A*rk+1
            double beta = 0.0;
            num = 0.0;
            for (int l = 0; l < n; l++) {
                num += y[l] * r[l];
            }
            beta = num / den;

            for (int i = 0; i < n; i++) {
                d[i] = r[i] - beta * d[i];
            }

            relativeError = Math.sqrt(residualNormSq) / bNorm;
            if (relativeError < tol) {
                validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data)
                        : 0.0;
                return new MatrixResult("Conjugate Gradient", relativeError, k + 1, true,
                        validationError, new SimpleMatrix(n, 1, true, x));
            }
        }

        validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
        return new MatrixResult("Conjugate Gradient", relativeError, MAX_ITER, false,
                validationError,
                new SimpleMatrix(n, 1, true, x));
    }
}
