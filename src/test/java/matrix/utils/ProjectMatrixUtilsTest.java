package matrix.utils;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.simple.SimpleMatrix;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

public class ProjectMatrixUtilsTest {
    
    private DMatrixSparseCSC symmetricMatrix;
    private DMatrixSparseCSC asymmetricMatrix;
    private DMatrixSparseCSC positiveDefiniteMatrix;
    private DMatrixSparseCSC nonSymmetricMatrix;
    private DMatrixSparseCSC negativeMatrix;
    private DMatrixSparseCSC sparseSymmetric;
    private SimpleMatrix vector;

    @Before
    public void setUp() {
        // Matrice simmetrica 3x3
        symmetricMatrix = createSparseMatrix(new double[][] {
            {4, 1, 2},
            {1, 3, 0},
            {2, 0, 5}
        });

        // Matrice non simmetrica
        asymmetricMatrix = createSparseMatrix(new double[][] {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}
        });

        // Matrice simmetrica definita positiva
        positiveDefiniteMatrix = createSparseMatrix(new double[][] {
            {2, 0, 0},
            {0, 2, 0},
            {0, 0, 2}
        });

        // Matrice non simmetrica (quindi non può essere definita positiva)
        nonSymmetricMatrix = createSparseMatrix(new double[][] {
            {1, 1, 0},
            {0, 1, 0},
            {0, 0, 1}
        });

        // Matrice definita negativa
        negativeMatrix = createSparseMatrix(new double[][] {
            {-2, 0},
            {0, -2}
        });

        // Matrice sparsa 100x100 simmetrica
        DMatrixSparseTriplet sparseTriplet = new DMatrixSparseTriplet(100, 100, 4);
        sparseTriplet.addItem(0, 0, 10.0);
        sparseTriplet.addItem(99, 99, 10.0);
        sparseTriplet.addItem(0, 99, 1.5);
        sparseTriplet.addItem(99, 0, 1.5);
        sparseSymmetric = new DMatrixSparseCSC(100, 100, 4);
        DConvertMatrixStruct.convert(sparseTriplet, sparseSymmetric);

        // Vettore di test
        vector = new SimpleMatrix(3, 1, true, new double[]{1.5, 2.5, 3.5});
    }

    // Test per isSymmetric
    @Test
    public void testIsSymmetricWithSymmetricMatrix() {
        assertTrue("La matrice dovrebbe essere simmetrica", ProjectMatrixUtils.isSymmetric(symmetricMatrix));
    }

    @Test
    public void testIsSymmetricWithAsymmetricMatrix() {
        assertFalse("La matrice non dovrebbe essere simmetrica", ProjectMatrixUtils.isSymmetric(asymmetricMatrix));
    }

    @Test
    public void testIsSymmetricWithIdentityMatrix() {
        DMatrixSparseCSC identity = createSparseMatrix(new double[][] {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        });
        assertTrue("La matrice identità dovrebbe essere simmetrica", ProjectMatrixUtils.isSymmetric(identity));
    }

    @Test
    public void testIsSymmetricWithSingleElementMatrix() {
        DMatrixSparseCSC singleElement = createSparseMatrix(new double[][]{{5}});
        assertTrue("Una matrice 1x1 dovrebbe sempre essere simmetrica", 
                   ProjectMatrixUtils.isSymmetric(singleElement));
    }

    @Test
    public void testIsSymmetricWithNonSquareMatrix() {
        DMatrixSparseCSC nonSquare = createSparseMatrix(new double[][] {
            {1, 2, 3},
            {4, 5, 6}
        });
        assertFalse("Una matrice non quadrata non può essere simmetrica", 
                    ProjectMatrixUtils.isSymmetric(nonSquare));
    }

    // Test per isPositiveDefinite
    @Test
    public void testIsPositiveDefiniteWithPositiveDefiniteMatrix() {
        assertTrue("La matrice dovrebbe essere definita positiva", 
                   ProjectMatrixUtils.isPositiveDefinite(positiveDefiniteMatrix));
    }

    @Test
    public void testIsPositiveDefiniteWithNonsymmetricMatrix() {
        assertFalse("Una matrice non simmetrica non può essere definita positiva", 
                    ProjectMatrixUtils.isPositiveDefinite(nonSymmetricMatrix));
    }

    @Test
    public void testIsPositiveDefiniteWithIdentityMatrix() {
        DMatrixSparseCSC identity = createSparseMatrix(new double[][] {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        });
        assertTrue("La matrice identità dovrebbe essere definita positiva", 
                   ProjectMatrixUtils.isPositiveDefinite(identity));
    }

    @Test
    public void testIsPositiveDefiniteWithDiagonalMatrix() {
        DMatrixSparseCSC diagonalPositive = createSparseMatrix(new double[][] {
            {5, 0, 0},
            {0, 3, 0},
            {0, 0, 2}
        });
        assertTrue("Una matrice diagonale simmetrica con elementi positivi dovrebbe essere definita positiva", 
                   ProjectMatrixUtils.isPositiveDefinite(diagonalPositive));
    }

    @Test
    public void testIsPositiveDefiniteWithNegativeMatrix() {
        assertFalse("Una matrice definita negativa non è definita positiva", 
                    ProjectMatrixUtils.isPositiveDefinite(negativeMatrix));
    }

    // Test per vectorToString
    @Test
    public void testVectorToStringNotNull() {
        String result = ProjectMatrixUtils.vectorToString(vector);
        assertNotNull("vectorToString non dovrebbe ritornare null", result);
    }

    @Test
    public void testVectorToStringWithNullVector() {
        String result = ProjectMatrixUtils.vectorToString(null);
        assertEquals("vectorToString con null dovrebbe ritornare stringa vuota", "", result);
    }

    @Test
    public void testVectorToStringNotEmpty() {
        String result = ProjectMatrixUtils.vectorToString(vector);
        assertFalse("vectorToString non dovrebbe ritornare una stringa vuota", result.isEmpty());
    }

    @Test
    public void testVectorToStringFormatting() {
        String result = ProjectMatrixUtils.vectorToString(vector);
        assertTrue("La stringa dovrebbe contenere il formato .6f", result.matches(".*\\d+\\.\\d{6}.*"));
    }

    @Test
    public void testVectorToStringWithZeroVector() {
        SimpleMatrix zeroVector = new SimpleMatrix(3, 1);
        String result = ProjectMatrixUtils.vectorToString(zeroVector);
        assertTrue("La stringa dovrebbe contenere zeri", result.contains("0.000000"));
    }

    @Test
    public void testVectorToStringContainsExtraSpaces() {
        // Verifica che la nuova implementazione aggiunge spazi tra gli elementi
        String result = ProjectMatrixUtils.vectorToString(vector);
        // Ogni elemento è formattato con %12.6f seguito da uno spazio
        // Poi c'è uno spazio extra se non raggiunge 10 elementi
        int spaceCount = result.split(" ").length - 1;
        assertTrue("Dovrebbe contenere molteplici spazi tra gli elementi", spaceCount >= vector.getNumRows());
    }

    // Test per isSymmetric con matrice densa rappresentata come CSC
    @Test
    public void testIsSymmetricWithDenseSparseMatrix() {
        DMatrixSparseCSC dense = createSparseMatrix(new double[][] {
            {1, 2, 3},
            {2, 4, 5},
            {3, 5, 6}
        });
        assertTrue("Matrice densa simmetrica dovrebbe essere riconosciuta", 
                   ProjectMatrixUtils.isSymmetric(dense));
    }

    // Test per isSymmetric con matrice sparsa simulata
    @Test
    public void testIsSymmetricWithSparseMatrix() {
        assertTrue("La matrice sparsa simmetrica dovrebbe essere riconosciuta", 
                ProjectMatrixUtils.isSymmetric(sparseSymmetric));
    }

    @Test
    public void testIsSymmetricWithTolerance() {
        DMatrixSparseCSC m = createSparseMatrix(new double[][]{
            {1, 2},
            {2 + 1e-13, 1}
        });
        assertTrue(ProjectMatrixUtils.isSymmetric(m));
    }

    // Test per isPositiveDefiniteWithSparseMatrix
    @Test
    public void testIsPositiveDefiniteWithSparseMatrix() {
        // Una matrice identità sparsa simulata è definita positiva
        DMatrixSparseCSC sparseIdentity = createSparseMatrix(new double[50][50]);
        for (int i = 0; i < 50; i++) {
            sparseIdentity.set(i, i, 1.0);
        }
        assertTrue("L'identità sparsa deve essere definita positiva", 
                ProjectMatrixUtils.isPositiveDefinite(sparseIdentity));
    }
    
    @Test(expected = ArithmeticException.class)
    public void testInverseDiagonalWithZeroOnDiagonal() {
        DMatrixSparseCSC zeroDiag = createSparseMatrix(new double[][]{
            {1, 0},
            {0, 0} // Lo zero qui deve far scattare l'errore
        });
        ProjectMatrixUtils.inverseDiagonal(zeroDiag);
    }

    private DMatrixSparseCSC createSparseMatrix(double[][] values) {
        int rows = values.length;
        int cols = rows == 0 ? 0 : values[0].length;
        int nonZero = 0;
        for (double[] row : values) {
            for (double value : row) {
                if (Math.abs(value) > 0.0) {
                    nonZero++;
                }
            }
        }

        DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(rows, cols, nonZero);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double value = values[i][j];
                if (Math.abs(value) > 0.0) {
                    triplet.addItem(i, j, value);
                }
            }
        }

        DMatrixSparseCSC matrix = new DMatrixSparseCSC(rows, cols, nonZero);
        DConvertMatrixStruct.convert(triplet, matrix);
        return matrix;
    }
}

