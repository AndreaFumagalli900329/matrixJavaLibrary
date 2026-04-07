package matrix.solver;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.simple.SimpleMatrix;
import org.junit.Before;
import org.junit.Test;
import matrix.utils.MatrixResult;
import static org.junit.Assert.*;

public class JacobiSolverTest {

    private SimpleMatrix matrix2x2;
    private SimpleMatrix matrix3x3;
    private SimpleMatrix matrix1x1;
    private SimpleMatrix zeroDiagonalMatrix;
    private SimpleMatrix nonConvergentMatrix;

    private static final double TOL = 1e-8; // Tolleranza richiesta dal progetto
    private static final double ACCEPTABLE_ERROR = 1e-5;

    @Before
    public void setUp() {
        // Matrice 2x2 simmetrica e definita positiva
        matrix2x2 = new SimpleMatrix(new double[][]{
            {4, 1},
            {1, 3}
        });

        // Matrice 3x3 simmetrica e definita positiva (diagonale dominante)
        matrix3x3 = new SimpleMatrix(new double[][]{
            {10, 1, -1},
            {1, 10, 2},
            {-1, 2, 10}
        });

        // Matrice 1x1
        matrix1x1 = new SimpleMatrix(new double[][]{{5.0}});

        // Matrice con zero sulla diagonale (per testare l'eccezione)
        zeroDiagonalMatrix = new SimpleMatrix(new double[][]{
            {0, 1},
            {1, 2}
        });

        // Matrice che non converge con Jacobi (non diagonale dominante)
        nonConvergentMatrix = new SimpleMatrix(new double[][]{
            {1, 2},
            {2, 1}
        });
    }

    /**
     * Helper per convertire SimpleMatrix in DMatrixSparseCSC come richiesto dal Solver
     */
    private DMatrixSparseCSC toSparse(SimpleMatrix sm) {
        DMatrixSparseCSC sparse = new DMatrixSparseCSC(sm.getNumRows(), sm.getNumCols(), 0);
        DConvertMatrixStruct.convert(sm.getDDRM(), sparse);
        return sparse;
    }

    /**
     * Helper per creare la soluzione esatta di soli 1 (Step 1 della consegna)
     */
    private SimpleMatrix createExactSol(int n) {
        SimpleMatrix x = new SimpleMatrix(n, 1);
        x.fill(1.0);
        return x;
    }

    @Test
    public void testJacobiConvergenceStepByStep() {
        Jacobi solver = new Jacobi();
        int n = matrix3x3.getNumRows();

        // Step 1: Creazione soluzione esatta [1, 1, 1]
        SimpleMatrix exactX = createExactSol(n);

        // Step 2: Creazione termine noto b = Ax
        SimpleMatrix b = matrix3x3.mult(exactX);

        // Step 3: Calcolo soluzione approssimata
        MatrixResult result = solver.solve(toSparse(matrix3x3), b, TOL, exactX);

        // Step 4: Validazione risultati
        assertTrue("Il metodo dovrebbe convergere", result.isConverged());
        assertTrue("L'errore di validazione (soluzione) deve essere basso", result.getValidationError() < ACCEPTABLE_ERROR);
        assertTrue("L'errore relativo (residuo) deve essere < tolleranza", result.getRelativeError() < TOL);
        System.out.println("Test Jacobi 3x3 - Iterazioni: " + result.getIterations());
    }

    @Test
    public void testJacobi1x1() {
        Jacobi solver = new Jacobi();
        SimpleMatrix exactX = createExactSol(1);
        SimpleMatrix b = matrix1x1.mult(exactX);

        MatrixResult result = solver.solve(toSparse(matrix1x1), b, TOL, exactX);

        assertTrue(result.isConverged());
        assertEquals(1.0, result.getSolution().get(0, 0), 1e-10);
        assertEquals(0.0, result.getValidationError(), 1e-15);
    }

    @Test(expected = ArithmeticException.class)
    public void testJacobiZeroDiagonal() {
        Jacobi solver = new Jacobi();
        SimpleMatrix b = new SimpleMatrix(2, 1, true, new double[]{1, 1});
        // Non serve la soluzione esatta qui, deve lanciare eccezione prima
        solver.solve(toSparse(zeroDiagonalMatrix), b, TOL, null);
    }

    @Test
    public void testJacobiWithZeroVector() {
        Jacobi solver = new Jacobi();
        SimpleMatrix zeroB = new SimpleMatrix(2, 1); // b = [0, 0]
        SimpleMatrix exactX = new SimpleMatrix(2, 1); // Soluzione attesa è [0, 0]
        
        MatrixResult result = solver.solve(toSparse(matrix2x2), zeroB, TOL, exactX);

        assertTrue("Convergenza istantanea con vettore nullo", result.isConverged());
        assertEquals(0, result.getIterations());
        assertEquals(0.0, result.getSolution().normF(), 1e-15);
    }

    @Test
    public void testJacobiNonConvergent() {
        Jacobi solver = new Jacobi();
        int n = nonConvergentMatrix.getNumRows();
        SimpleMatrix exactX = createExactSol(n);
        SimpleMatrix b = nonConvergentMatrix.mult(exactX);

        MatrixResult result = solver.solve(toSparse(nonConvergentMatrix), b, 1e-10, exactX);

        // Jacobi non converge su questa matrice perché non è diagonale dominante
        assertFalse(result.isConverged());
        assertEquals(Solver.MAX_ITER, result.getIterations()); // Deve fermarsi a 20.000 iterazioni 
    }

    @Test
    public void testJacobiDeterministic() {
        Jacobi solver = new Jacobi();
        DMatrixSparseCSC sparse = toSparse(matrix2x2);
        SimpleMatrix exactX = createExactSol(2);
        SimpleMatrix b = matrix2x2.mult(exactX);

        MatrixResult res1 = solver.solve(sparse, b, TOL, exactX);
        MatrixResult res2 = solver.solve(sparse, b, TOL, exactX);

        assertEquals("Il numero di iterazioni deve essere identico", res1.getIterations(), res2.getIterations());
        assertEquals("L'errore finale deve essere identico", res1.getRelativeError(), res2.getRelativeError(), 1e-18);
    }
}