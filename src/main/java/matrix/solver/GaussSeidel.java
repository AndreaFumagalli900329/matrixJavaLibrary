package matrix.solver;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;
import org.ejml.sparse.csc.CommonOps_DSCC;

import matrix.utils.MatrixResult;
import matrix.utils.ProjectMatrixUtils;

public class GaussSeidel implements Solver {

    public GaussSeidel() {
    }

    @Override
    public MatrixResult solve(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol) {
        int n = matrix.numRows;

        double[] b = vector.getDDRM().data;
        double[] x = new double[n];
        double validationError = 0.0;
        double relativeError = 0.0;

        double bNorm = 0.0;
        for (double val : b)
            bNorm += val * val;
        bNorm = Math.sqrt(bNorm);

        if (bNorm == 0) {
            return new MatrixResult("Gauss-Seidel", relativeError, 0, true, validationError,
                    new SimpleMatrix(n, 1, true, x));
        }

        double[] invDiag = ProjectMatrixUtils.inverseDiagonal(matrix);

        
        // Salvo quanti elementi non zero ci sono per ogni riga (per costruire la matrice in formato CSR)
        int[] rowPtr = new int[n + 1];

        for(int p = 0; p < matrix.nz_length; p++) {
            rowPtr[matrix.nz_rows[p]]++;
        }

        // Calcolo dei puntatori cumulativi per ogni riga
        int sum = 0;
        for (int i = 0; i < n; i++) {
            int temp = rowPtr[i];
            rowPtr[i] = sum;
            sum += temp;
        }
        rowPtr[n] = sum; // Ultimo elemento punta alla fine dell'array nz_rows

        // Array di colonne corrispondenti agli elementi non zero, ordinati per riga
        int[] csrColIdx = new int[matrix.nz_length];
        double[] csrValues = new double[matrix.nz_length];
        int[] currentOffsets = rowPtr.clone();

        for (int c = 0; c < n; c++) {
            int start = matrix.col_idx[c];
            int end = matrix.col_idx[c + 1];
            for (int p = start; p < end; p++) {
                int r = matrix.nz_rows[p];
                int destPos = currentOffsets[r];
                csrColIdx[destPos] = c;
                csrValues[destPos] = matrix.nz_values[p];
                currentOffsets[r]++;
            }
        }

        double[] Ax = new double[n];
        DMatrixRMaj Ax_final_mat = DMatrixRMaj.wrap(n, 1, Ax);
        DMatrixRMaj x_mat = DMatrixRMaj.wrap(n, 1, x);

        for (int k = 0; k < MAX_ITER; k++) {
            // Aggiorno x[i] e lo uso subito
            for (int i = 0; i < n; i++) {
                double Ax_i = 0.0;

                /* Prodotto riga i-esima
                // Usiamo x[j] che contiene sia valori "vecchi" (se j > i) sia valori "nuovi"
                // (se j < i)
                for (int j = 0; j < n; j++) {
                    double val = matrix.get(i, j);
                    if (val != 0) {
                        Ax_i += val * x[j];
                    }
                }*/

                int start = rowPtr[i];
                int end = rowPtr[i + 1];
                for (int p = start; p < end; p++) {
                    int j = csrColIdx[p];
                    Ax_i += csrValues[p] * x[j];
                }

                // Correzione in-place: x_new = x_old + (b - Ax_old) / A_ii
                x[i] = x[i] + (b[i] - Ax_i) * invDiag[i];
            }

            //Uso la libreria per calcolare la Ax con le x aggiornate
            CommonOps_DSCC.mult(matrix, x_mat, Ax_final_mat);
            // Calcolo Residuo Manuale (dopo aver completato l'aggiornamento di tutto il
            // vettore x)
            double residualNormSq = 0.0;
            for (int i = 0; i < n; i++) {
                double r = b[i] - Ax[i];
                residualNormSq += r * r;
            }

            relativeError = Math.sqrt(residualNormSq) / bNorm;

            if (relativeError < tol) {
                validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data)
                        : 0.0;
                return new MatrixResult("Gauss-Seidel", relativeError, k + 1, true, validationError,
                        new SimpleMatrix(n, 1, true, x));
            }
        }

        validationError = (exactSol != null) ? ProjectMatrixUtils.validationError(x, exactSol.getDDRM().data) : 0.0;
        return new MatrixResult("Gauss-Seidel", relativeError, MAX_ITER, false, validationError,
                new SimpleMatrix(n, 1, true, x));
    }
}