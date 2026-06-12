package matrix.solver;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.simple.SimpleMatrix;
import org.junit.Before;
import org.junit.Test;

import matrix.utils.MatrixResult;

public class ConjugateGradientSolverTest {

    private SimpleMatrix matrix2x2;
    private SimpleMatrix matrix3x3;
    private SimpleMatrix matrix1x1;
    private SimpleMatrix spdWithNegativeOffDiag;
    private SimpleMatrix nonSymmetricMatrix; // Aggiunta per test di robustezza
    
    private static final double TOL = 1e-8;
    private static final double ACCEPTABLE_VALIDATION_ERROR = 1e-6;

    @Before
    public void setUp() {
        matrix2x2 = new SimpleMatrix(new double[][]{
            {4, 1},
            {1, 3}
        });

        matrix3x3 = new SimpleMatrix(new double[][]{
            {10, 1, -1},
            {1, 10, 2},
            {-1, 2, 10}
        });

        // Rinominata per chiarezza: è SPD ma ha elementi negativi fuori diagonale
        spdWithNegativeOffDiag = new SimpleMatrix(new double[][]{
            {4, -1, 0},
            {-1, 4, -1},
            {0, -1, 4}
        });

        matrix1x1 = new SimpleMatrix(new double[][]{{5.0}});

        // Matrice non simmetrica: il CG puro DEVE fallire o non garantire la convergenza
        nonSymmetricMatrix = new SimpleMatrix(new double[][]{
            {4, 5},
            {1, 3}
        });
    }

    private DMatrixSparseCSC toSparse(SimpleMatrix sm) {
        DMatrixSparseCSC sparse = new DMatrixSparseCSC(sm.getNumRows(), sm.getNumCols(), 0);
        DConvertMatrixStruct.convert(sm.getDDRM(), sparse);
        return sparse;
    }

    private SimpleMatrix createExactSol(int n) {
        SimpleMatrix x = new SimpleMatrix(n, 1);
        x.fill(1.0);
        return x;
    }

    @Test
    public void testConjugateGradientConvergenceStepByStep() {
        ConjugateGradient solver = new ConjugateGradient();
        int n = matrix3x3.getNumRows();

        SimpleMatrix exactX = createExactSol(n);
        SimpleMatrix b = matrix3x3.mult(exactX);

        MatrixResult result = solver.solve(toSparse(matrix3x3), b, TOL, exactX);

        assertTrue("Il metodo Conjugate Gradient dovrebbe convergere su matrice SPD", result.isConverged());
        assertEquals("Il nome del metodo deve essere Conjugate Gradient", "Conjugate Gradient", result.getMethodName());
        assertTrue("L'errore relativo deve essere minore della tolleranza", result.getRelativeError() < TOL);
        assertTrue("L'errore di validazione deve essere basso", result.getValidationError() < ACCEPTABLE_VALIDATION_ERROR);
        
        // FISSA TEORICA: Il CG deve convergere in <= n iterazioni (tolleriamo +1 per arrotondamenti)
        assertTrue("Il numero di iterazioni deve essere positivo", result.getIterations() > 0);
        assertTrue("Il CG su matrici piccole deve convergere in <= n+1 iterazioni", result.getIterations() <= n + 1);
        
        assertNotNull("La soluzione non deve essere null", result.getSolution());
    }

    @Test
    public void testConjugateGradient1x1() {
        ConjugateGradient solver = new ConjugateGradient();
        SimpleMatrix exactX = createExactSol(1);
        SimpleMatrix b = matrix1x1.mult(exactX);

        MatrixResult result = solver.solve(toSparse(matrix1x1), b, TOL, exactX);

        assertTrue(result.isConverged());
        assertEquals(1.0, result.getSolution().get(0, 0), 1e-10);
        assertEquals(0.0, result.getValidationError(), 1e-15);
        assertEquals("Nel caso 1x1 deve metterci esattamente 1 iterazione", 1, result.getIterations());
    }

    @Test
    public void testConjugateGradientWithZeroVector() {
        ConjugateGradient solver = new ConjugateGradient();
        SimpleMatrix zeroB = new SimpleMatrix(2, 1);
        SimpleMatrix exactX = new SimpleMatrix(2, 1);

        MatrixResult result = solver.solve(toSparse(matrix2x2), zeroB, TOL, exactX);

        assertTrue("Con vettore nullo il metodo deve considerarsi convergente", result.isConverged());
        assertEquals("Il numero di iterazioni deve essere zero", 0, result.getIterations());
        assertEquals("La soluzione attesa è il vettore nullo", 0.0, result.getSolution().normF(), 1e-15);
    }

    @Test
    public void testConjugateGradientWithoutExactSolution() {
        ConjugateGradient solver = new ConjugateGradient();
        int n = matrix2x2.getNumRows();
        SimpleMatrix exactX = createExactSol(n);
        SimpleMatrix b = matrix2x2.mult(exactX);

        MatrixResult result = solver.solve(toSparse(matrix2x2), b, TOL, null);

        assertTrue(result.isConverged());
        assertEquals("Senza exactSol validationError deve essere 0.0", 0.0, result.getValidationError(), 1e-15);
    }

    @Test
    public void testConjugateGradientDeterministic() {
        ConjugateGradient solver = new ConjugateGradient();
        DMatrixSparseCSC sparse = toSparse(matrix2x2);
        SimpleMatrix exactX = createExactSol(2);
        SimpleMatrix b = matrix2x2.mult(exactX);

        MatrixResult result1 = solver.solve(sparse, b, TOL, exactX);
        MatrixResult result2 = solver.solve(sparse, b, TOL, exactX);

        assertEquals("Il numero di iterazioni deve essere identico", result1.getIterations(), result2.getIterations());
        assertEquals("L'errore relativo deve essere identico", result1.getRelativeError(), result2.getRelativeError(), 1e-18);
    }

    @Test
    public void testConjugateGradientWithNegativeValuesSPD() {
        ConjugateGradient solver = new ConjugateGradient();
        int n = spdWithNegativeOffDiag.getNumRows();
        SimpleMatrix exactX = createExactSol(n);
        SimpleMatrix b = spdWithNegativeOffDiag.mult(exactX);

        MatrixResult result = solver.solve(toSparse(spdWithNegativeOffDiag), b, TOL, exactX);

        assertTrue("Il metodo deve convergere anche su una matrice SPD con valori negativi fuori diagonale", result.isConverged());
        assertTrue(result.getIterations() <= n + 1);
    }

    // --- NUOVO TEST CRITICO ---
    @Test
    public void testConjugateGradientFailsOnNonSymmetricMatrix() {
        ConjugateGradient solver = new ConjugateGradient();
        int n = nonSymmetricMatrix.getNumRows();
        SimpleMatrix exactX = createExactSol(n);
        SimpleMatrix b = nonSymmetricMatrix.mult(exactX);

        MatrixResult result = solver.solve(toSparse(nonSymmetricMatrix), b, TOL, exactX);

        // A seconda di come gestisci il fallimento nel solver, questo test si assicura 
        // che l'algoritmo non restituisca un flag di finta convergenza corretta.
        assertFalse("Il CG non deve dichiararsi felicemente convergente su matrici non asimmetriche", 
                    result.isConverged() && result.getValidationError() < ACCEPTABLE_VALIDATION_ERROR);
    }
}