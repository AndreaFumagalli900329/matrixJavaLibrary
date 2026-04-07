package matrix.solver;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;

import matrix.utils.MatrixResult;

public class Gradient implements Solver {

    public Gradient() {
    }

    @Override
    public MatrixResult solve(DMatrixSparseCSC matrix, SimpleMatrix vector, double tol, SimpleMatrix exactSol) {
        SimpleMatrix x = new SimpleMatrix(matrix.getNumRows(), 1);
        return new MatrixResult("Gradient", 0.0, 0, 0.0, false, 0.0, x);
    }
}
