package matrix;

import java.util.Scanner;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;
import matrix.solver.Jacobi;
import matrix.solver.GaussSeidel;
import matrix.utils.MatrixResult;
import matrix.utils.ProjectMatrixUtils;

public class Main {
    private static final String FILE_NAME = "spa1.mtx";

    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);
        
        System.out.print("Inserire l'esponente della tolleranza (es. -4 per 10^-4): ");
        if (!input.hasNextInt()) {
            System.err.println("Errore: Inserire un numero intero.");
            input.close();
            return;
        }
        int exp = input.nextInt();
        double tol = Math.pow(10, exp);
        System.out.println("Configurata tolleranza: " + tol);

        try {
            DMatrixSparseCSC matrix = ProjectMatrixUtils.importMatrix(FILE_NAME);
            int n = matrix.numRows;
            System.out.println("Matrice " + FILE_NAME + " caricata con successo (" + n + "x" + n + ")");

            SimpleMatrix exactSol = new SimpleMatrix(n, 1);
            exactSol.fill(1.0);
            SimpleMatrix b = SimpleMatrix.wrap(matrix).mult(exactSol);

            if (ProjectMatrixUtils.isPositiveDefinite(matrix)) {
                Jacobi jacobiSolver = new Jacobi();
                MatrixResult jacobiResult = jacobiSolver.solve(matrix, b, tol, exactSol);
                System.out.println(jacobiResult);

                System.out.println("------------------------------------------------------------");

                GaussSeidel gsSolver = new GaussSeidel();
                MatrixResult gsResult = gsSolver.solve(matrix, b, tol, exactSol);
                System.out.println(gsResult);

            } else {
                System.out.println("La matrice non rispetta le condizioni (Simmetria/Positività).");
            }

        } catch (Exception e) {
            System.err.println("ERRORE DURANTE L'ESECUZIONE: " + e.getMessage());
            e.printStackTrace();
        } finally {
            input.close();
        }
    }
}