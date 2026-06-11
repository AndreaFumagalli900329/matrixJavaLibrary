package matrix.solver;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;
import org.ejml.sparse.csc.CommonOps_DSCC;

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

        if (bNorm == 0)
            return new MatrixResult("Jacobi", relativeError, 0, true, validationError,
                    new SimpleMatrix(n, 1, true, x));

        double[] invDiag = ProjectMatrixUtils.inverseDiagonal(matrix);

        double[] Ax = new double[n];
        DMatrixRMaj Ax_mat = DMatrixRMaj.wrap(n, 1, Ax);

        for (int k = 0; k < MAX_ITER; k++) {
            double residualNormSq = 0.0;

            DMatrixRMaj x_mat = DMatrixRMaj.wrap(n, 1, x);

            CommonOps_DSCC.mult(matrix, x_mat, Ax_mat);
            
            for (int i = 0; i < n; i++) {
                // double Ax_i = 0.0;

                // for (int j = 0; j < n; j++) {
                //     double val = matrix.get(i, j);
                //     if (val != 0) {
                //         Ax_i += val * x[j];
                //     }
                // }

                double r_i = b[i] - Ax[i];
                residualNormSq += r_i * r_i;

                nextX[i] = x[i] + (r_i * invDiag[i]);
            }

            double[] temp = x;
            x = nextX;
            nextX = temp;

            relativeError = Math.sqrt(residualNormSq) / bNorm;

            if (relativeError < tol) {
                validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data)
                        : 0.0;
                return new MatrixResult("Jacobi", relativeError, k + 1, true, validationError,
                        new SimpleMatrix(n, 1, true, x));
            }
        }
        validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
        return new MatrixResult("Jacobi", relativeError, MAX_ITER, false, validationError,
                new SimpleMatrix(n, 1, true, x));
    }
}
