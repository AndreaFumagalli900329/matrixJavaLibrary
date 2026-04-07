package matrix.solver;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;

import matrix.utils.MatrixResult;

public interface Solver {
    int MAX_ITER  = 100000;

    MatrixResult solve(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol);

    default MatrixResult solve(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol) {
        return solve(matrix, vector, tol, null);
    }
}
