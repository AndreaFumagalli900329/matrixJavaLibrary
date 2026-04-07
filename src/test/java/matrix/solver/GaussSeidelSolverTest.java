package matrix.solver;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.simple.SimpleMatrix;
import org.junit.Before;
import org.junit.Test;
import matrix.utils.MatrixResult;
import static org.junit.Assert.*;

public class GaussSeidelSolverTest {

    private SimpleMatrix matrix2x2;
    private SimpleMatrix diagonallyDominantMatrix;
    private static final double TOL = 1e-8;

    @Before
    public void setUp() {
        // Matrice simmetrica e definita positiva [cite: 5]
        matrix2x2 = new SimpleMatrix(new double[][]{
            {4, 1},
            {1, 3}
        });

        diagonallyDominantMatrix = new SimpleMatrix(new double[][]{
            {10, 1, -1},
            {1, 10, 2},
            {-1, 2, 10}
        });
    }

    private DMatrixSparseCSC toSparse(SimpleMatrix sm) {
        DMatrixSparseCSC sparse = new DMatrixSparseCSC(sm.getNumRows(), sm.getNumCols(), 0);
        DConvertMatrixStruct.convert(sm.getDDRM(), sparse);
        return sparse;
    }

    private SimpleMatrix createExactSol(int n) {
        SimpleMatrix x = new SimpleMatrix(n, 1);
        x.fill(1.0); // Step 1: soluzione con tutte entrate pari a 1 [cite: 33, 34]
        return x;
    }

    @Test
    public void testGaussSeidelConvergenceStepByStep() {
        GaussSeidel solver = new GaussSeidel();
        int n = matrix2x2.getNumRows();

        // Step 1 & 2: b = Ax [cite: 35, 36]
        SimpleMatrix exactX = createExactSol(n);
        SimpleMatrix b = matrix2x2.mult(exactX);

        // Step 3: Calcolo soluzione [cite: 37]
        MatrixResult result = solver.solve(toSparse(matrix2x2), b, TOL, exactX);

        // Step 4: Verifica errore relativo soluzione [cite: 38]
        assertTrue("Dovrebbe convergere", result.isConverged());
        assertTrue("L'errore di validazione deve essere basso", result.getValidationError() < 1e-5);
        assertEquals("Gauss-Seidel", result.getMethodName());
    }

    @Test(expected = ArithmeticException.class)
    public void testGaussSeidelZeroDiagonal() {
        GaussSeidel solver = new GaussSeidel();
        SimpleMatrix zeroDiag = new SimpleMatrix(new double[][]{{0, 1}, {1, 2}});
        // exactSol può essere null qui perché l'eccezione deve scattare prima
        solver.solve(toSparse(zeroDiag), new SimpleMatrix(2, 1, true, new double[]{1, 1}), TOL, null);
    }

    @Test
    public void testGaussSeidelWithZeroVector() {
        GaussSeidel solver = new GaussSeidel();
        int n = 3;
        SimpleMatrix zeroB = new SimpleMatrix(n, 1);
        SimpleMatrix exactX = new SimpleMatrix(n, 1); // Soluzione attesa: tutto zero
        
        MatrixResult result = solver.solve(toSparse(diagonallyDominantMatrix), zeroB, TOL, exactX);

        assertTrue(result.isConverged());
        assertEquals(0, result.getIterations()); // [cite: 19]
        assertEquals(0.0, result.getSolution().normF(), 1e-15);
    }

    @Test
    public void testFasterThanJacobi() {
        DMatrixSparseCSC sparse = toSparse(diagonallyDominantMatrix);
        int n = diagonallyDominantMatrix.getNumRows();
        SimpleMatrix exactX = createExactSol(n);
        SimpleMatrix b = diagonallyDominantMatrix.mult(exactX);

        // Entrambi devono partire da vettore nullo [cite: 19]
        MatrixResult resJacobi = new Jacobi().solve(sparse, b, TOL, exactX);
        MatrixResult resGS = new GaussSeidel().solve(sparse, b, TOL, exactX);

        if (resJacobi.isConverged() && resGS.isConverged()) {
            // Teoricamente GS converge circa il doppio più velocemente di Jacobi
            assertTrue("Gauss-Seidel dovrebbe richiedere meno iterazioni di Jacobi", 
                       resGS.getIterations() <= resJacobi.getIterations());
        }
    }
}