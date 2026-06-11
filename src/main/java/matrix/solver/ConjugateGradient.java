package matrix.solver;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;
import org.ejml.sparse.csc.CommonOps_DSCC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;

import matrix.utils.MatrixResult;
import matrix.utils.ProjectMatrixUtils;

public class ConjugateGradient implements Solver {

    // array, hybrid o ejml
    private String mode; 
    public ConjugateGradient() {
        this.mode = "array";
    }

    public ConjugateGradient(String mode) {
        this.mode = mode;
    }

    @Override
    public MatrixResult solve(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol) {
        if ("array".equals(mode)) {
            return solveArray(matrix, vector, tol, exactSol);
        } else if ("hybrid".equals(mode)) {
            return solveHybrid(matrix, vector, tol, exactSol);
        } else if ("ejml".equals(mode)) {
            return solveEJML(matrix, vector, tol, exactSol);
        } else {
            throw new IllegalArgumentException("Modalità non supportata: " + mode);
        }
    }

    private MatrixResult solveArray(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol) {
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
            return new MatrixResult("Conjugate Gradient Array", relativeError, 0, true, validationError,
                    new SimpleMatrix(n, 1, true, x));
        }

        r = b.clone();
        d = r.clone();

        for (int k = 0; k < MAX_ITER; k++) {
            double residualNormSq = 0.0;
            num = 0.0;
            den = 0.0;
            // calcolo y = A * d
            for (int i = 0; i < n; i++) {
                double Ad_i = 0.0;
                for (int j = 0; j < n; j++) {
                    double val = matrix.get(i, j);
                    if (val != 0) {
                        Ad_i += val * d[j];
                    }
                }
                y[i] = Ad_i;
            }

            // calcolo a
            for (int i = 0; i < n; i++) {
                num += d[i] * r[i];
                den += y[i] * d[i];
            }
            a = num / den;

            // aggiorno x e r
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

            // aggiorno d
            for (int i = 0; i < n; i++) {
                d[i] = r[i] - beta * d[i];
            }

            relativeError = Math.sqrt(residualNormSq) / bNorm;
            if (relativeError < tol) {
                validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data)
                        : 0.0;
                return new MatrixResult("Conjugate Gradient Array", relativeError, k + 1, true,
                        validationError, new SimpleMatrix(n, 1, true, x));
            }
        }

        validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
        return new MatrixResult("Conjugate Gradient Array", relativeError, MAX_ITER, false, validationError,
                new SimpleMatrix(n, 1, true, x));
    }

    private MatrixResult solveHybrid(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol) {
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
            return new MatrixResult("Conjugate Gradient Hybrid", relativeError, 0, true, validationError,
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
                return new MatrixResult("Conjugate Gradient Hybrid", relativeError, k + 1, true,
                        validationError, new SimpleMatrix(n, 1, true, x));
            }
        }

        validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
        return new MatrixResult("Conjugate Gradient Hybrid", relativeError, MAX_ITER, false,
                validationError,
                new SimpleMatrix(n, 1, true, x));
    }

    private MatrixResult solveEJML(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol) {
        int n = matrix.numRows;

        // inizializzo sfruttando la libreria
        DMatrixRMaj b = vector.getDDRM();
        DMatrixRMaj x = new DMatrixRMaj(n, 1);
        DMatrixRMaj r = b.copy(); // r_0 = b
        DMatrixRMaj d = r.copy(); // d_0 = r_0
        DMatrixRMaj y = new DMatrixRMaj(n, 1);

        double bNorm = NormOps_DDRM.normF(b);
        if (bNorm == 0)
            return new MatrixResult("Conjugate Gradient EJML", 0.0, 0, true, 0.0, new SimpleMatrix(x));


        for (int k = 0; k < MAX_ITER; k++) {
            // y = A * d
            CommonOps_DSCC.mult(matrix, d, y);

            // calcolo a
            double numAlpha = CommonOps_DDRM.dot(d, r);
            double den = CommonOps_DDRM.dot(y, d);
            double a = numAlpha / den;

            // aggiorno x
            CommonOps_DDRM.addEquals(x, a, d);

            // aggiorno r
            CommonOps_DDRM.addEquals(r, -a, y);

            // caloclo beta
            double numBeta = CommonOps_DDRM.dot(y, r);
            double beta = numBeta / den;

            // aggiorno d
            // add(r, -beta, d, d) --> d = r - beta * d
            CommonOps_DDRM.add(r, -beta, d, d);

            // Controllo Convergenza
            double residualNorm = NormOps_DDRM.normF(r);
            double relativeError = residualNorm / bNorm;

            if (relativeError < tol) {
                double valError = (exactSol != null)
                        ? ProjectMatrixUtils.validationError(x.data, exactSol.getDDRM().data)
                        : 0.0;
                return new MatrixResult("Conjugate Gradient EJML", relativeError, k + 1, true, valError,
                        new SimpleMatrix(x));
            }

        }

        double valError = (exactSol != null) ? ProjectMatrixUtils.validationError(x.data, exactSol.getDDRM().data)
                : 0.0;
        return new MatrixResult("Conjugate Gradient EJML", NormOps_DDRM.normF(r) / bNorm, MAX_ITER,
                false,
                valError, new SimpleMatrix(x));
    }
}
