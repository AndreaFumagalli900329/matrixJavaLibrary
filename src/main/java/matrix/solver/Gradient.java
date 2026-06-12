package matrix.solver;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;
import org.ejml.sparse.csc.CommonOps_DSCC;

import matrix.utils.MatrixResult;
import matrix.utils.ProjectMatrixUtils;

public class Gradient implements Solver {

    public Gradient() {
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
        double[] p = new double[n];

        double bNorm = 0.0;
        for (double val : b)
            bNorm += val * val;
        bNorm = Math.sqrt(bNorm);

        if (bNorm == 0) {
            return new MatrixResult("Gradient", relativeError, 0, true, validationError,
                    new SimpleMatrix(n, 1, true, x));
        }

        r = b.clone();
        DMatrixRMaj p_mat = DMatrixRMaj.wrap(n, 1, p);
        DMatrixRMaj r_mat = DMatrixRMaj.wrap(n, 1, r);

        for (int k = 0; k < MAX_ITER; k++) {
            double residualNormSq = 0.0;
            num = 0.0;
            den = 0.0;

            CommonOps_DSCC.mult(matrix, r_mat, p_mat);
            for (int i = 0; i < n; i++) {
                num += r[i] * r[i];
                den += p[i] * r[i];
            }
            a = num / den;

            for (int i = 0; i < n; i++) {
                x[i] = x[i] + a * r[i];
                r[i] = r[i] - a * p[i];
                residualNormSq += r[i] * r[i];
            }

            relativeError = Math.sqrt(residualNormSq) / bNorm;
            if (relativeError < tol) {
                validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data)
                        : 0.0;
                return new MatrixResult("Gradient", relativeError, k + 1, true, validationError,
                        new SimpleMatrix(n, 1, true, x));
            }
        }
        
        validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
        return new MatrixResult("Gradient", relativeError, MAX_ITER, false, validationError,
                new SimpleMatrix(n, 1, true, x));
    }
}
