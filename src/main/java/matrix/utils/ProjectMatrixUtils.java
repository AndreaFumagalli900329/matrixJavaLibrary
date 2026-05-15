package matrix.utils;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.sparse.csc.MatrixFeatures_DSCC;
import org.ejml.sparse.csc.factory.DecompositionFactory_DSCC;
import org.ejml.interfaces.decomposition.CholeskySparseDecomposition_F64;
import org.ejml.simple.SimpleMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Locale;
import java.util.Scanner;

public final class ProjectMatrixUtils {

    private ProjectMatrixUtils() {}
    
    public static DMatrixSparseCSC importMatrix(String fileName){
        File file = new File(fileName);
        
        try (Scanner sc = new Scanner(file)){
            sc.useLocale(Locale.US);
            int rows = 0, cols = 0, nonZero = 0;

            while(sc.hasNextLine()){
                String line = sc.nextLine().trim();
                if(line.startsWith("%") || line.isEmpty()) continue;

                Scanner lineScanner = new Scanner(line);
                lineScanner.useLocale(Locale.US);
                rows = lineScanner.nextInt();
                cols = lineScanner.nextInt();
                nonZero = lineScanner.nextInt();
                lineScanner.close();
                break;
            }

            // OTTIMIZZAZIONE: Usiamo Triplet per caricamento veloce
            DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(rows, cols, nonZero);

            int count = 0;
            while(sc.hasNext() && count < nonZero) {
                int r = sc.nextInt() - 1;
                int c = sc.nextInt() - 1;
                double value = sc.nextDouble();
                triplet.addItem(r, c, value); // Molto più veloce di matrix.set
                count++;
            }

            // Convertiamo da Triplet a CSC (formato efficiente per calcoli)
            DMatrixSparseCSC matrix = new DMatrixSparseCSC(rows, cols, nonZero);
            DConvertMatrixStruct.convert(triplet, matrix);
            
            return matrix;
        } catch (FileNotFoundException e) {
            throw new IllegalArgumentException("File non trovato: " + fileName, e);
        }
    }

    public static boolean isSymmetric(DMatrixSparseCSC matrix) {
        return MatrixFeatures_DSCC.isSymmetric(matrix, 1e-12);
    }

    @SuppressWarnings("unchecked")
    public static boolean isPositiveDefinite(DMatrixSparseCSC matrix) {
        if (!isSymmetric(matrix)) return false;

        CholeskySparseDecomposition_F64<DMatrixSparseCSC> chol = DecompositionFactory_DSCC.cholesky();
        return chol.decompose(matrix);
    }

    public static double[] inverseDiagonal(DMatrixSparseCSC matrix) {
        int n = matrix.numRows;
        double[] invDiag = new double[n];
        for (int i = 0; i < n; i++) {
            double diag = matrix.get(i, i);
            if (Math.abs(diag) < 1e-15) {
                throw new ArithmeticException("Zero o valore troppo piccolo sulla diagonale alla riga " + i);
            }
            invDiag[i] = 1.0 / diag;
        }
        return invDiag;
    }

    public static double validationError(double[] xComputed, double[] xExact) {
        double diffNormSq = 0.0;
        double exactNormSq = 0.0;
        for (int i = 0; i < xComputed.length; i++) {
            double diff = xExact[i] - xComputed[i];
            diffNormSq += diff * diff;
            exactNormSq += xExact[i] * xExact[i];
        }
        return Math.sqrt(diffNormSq) / Math.sqrt(exactNormSq);
    }

    public static String vectorToString(SimpleMatrix v) {
        if (v == null) return "";
        
        double[] data = v.getDDRM().data;
        int numElements = v.getNumElements(); 
        
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < numElements; i++) {
            sb.append(String.format(Locale.US, "%12.6f", data[i]));
            
            if ((i + 1) % 5 == 0) {
                sb.append("\n");
            } else if (i < numElements - 1) {
                sb.append(" ");
            }
        }
        return sb.toString();
    }
}