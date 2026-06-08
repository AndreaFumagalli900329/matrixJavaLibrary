// package matrix;

// import java.util.Scanner;
// import org.ejml.data.DMatrixSparseCSC;
// import org.ejml.simple.SimpleMatrix;
// import matrix.solver.*;
// import matrix.utils.MatrixResult;
// import matrix.utils.ProjectMatrixUtils;

// public class Main {
//     private static final String FILE_NAME = "matrixJavaLibrary/spa2.mtx";

//     public static void main(String[] args) {
//         Scanner input = new Scanner(System.in);

//         System.out.print("Inserire l'esponente della tolleranza (es. -4 per 10^-4): ");
//         if (!input.hasNextInt()) {
//             System.err.println("Errore: Inserire un numero intero.");
//             input.close();
//             return;
//         }
//         int exp = input.nextInt();
//         double tol = Math.pow(10, exp);
//         System.out.println("Configurata tollera-4nza: " + tol);

//         try {
//             DMatrixSparseCSC matrix = ProjectMatrixUtils.importMatrix(FILE_NAME);
//             int n = matrix.numRows;
//             System.out.println("Matrice " + FILE_NAME + " caricata con successo (" + n + "x" + n + ")");

//             SimpleMatrix exactSol = new SimpleMatrix(n, 1);
//             exactSol.fill(1.0);
//             SimpleMatrix b = SimpleMatrix.wrap(matrix).mult(exactSol);

//             if (ProjectMatrixUtils.isPositiveDefinite(matrix)) {
//                 Jacobi jacobiSolver = new Jacobi();
//                 MatrixResult jacobiResult = jacobiSolver.solve(matrix, b, tol, exactSol);
//                 System.out.println(jacobiResult);

//                 System.out.println("------------------------------------------------------------");

//                 GaussSeidel gsSolver = new GaussSeidel();
//                 MatrixResult gsResult = gsSolver.solve(matrix, b, tol, exactSol);
//                 System.out.println(gsResult);

//                 System.out.println("------------------------------------------------------------");

//                 Gradient gradientSolver = new Gradient();
//                 MatrixResult gradientResult = gradientSolver.solve(matrix, b, tol, exactSol);
//                 System.out.println(gradientResult);

//                 System.out.println("------------------------------------------------------------");

//                 ConjugateGradient conjugateGradientSolver = new ConjugateGradient("array");
//                 MatrixResult conjugateGradientResult = conjugateGradientSolver.solve(matrix, b, tol, exactSol);
//                 System.out.println(conjugateGradientResult);

//                 System.out.println("------------------------------------------------------------");

//                 conjugateGradientSolver = new ConjugateGradient("hybrid");
//                 conjugateGradientResult = conjugateGradientSolver.solve(matrix, b, tol, exactSol);
//                 System.out.println(conjugateGradientResult);

//                 System.out.println("------------------------------------------------------------");

//                 conjugateGradientSolver = new ConjugateGradient("ejml");
//                 conjugateGradientResult = conjugateGradientSolver.solve(matrix, b, tol, exactSol);
//                 System.out.println(conjugateGradientResult);

//             } else {
//                 System.out.println("La matrice non rispetta le condizioni (Simmetria/Positività).");
//             }

//         } catch (Exception e) {
//             System.err.println("ERRORE DURANTE L'ESECUZIONE: " + e.getMessage());
//             e.printStackTrace();
//         } finally {
//             input.close();
//         }
//     }
// }

package matrix;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;
import matrix.solver.*;
import matrix.utils.MatrixResult;
import matrix.utils.ProjectMatrixUtils;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.io.File;

public class Main extends JFrame {

    private File selectedFile = null;
    private JTextArea consoleOutput;
    private JTextField tolInputField;
    private JTextField fileField;
    private JButton runButton;

    public Main() {
        setTitle("Metodi del Calcolo Scientifico - Progetto 1 Cosnegna Alternativa");
        setSize(800, 600);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLocationRelativeTo(null);
        setLayout(new BorderLayout(10, 10));

        // Pannello Input (File, Tol e Run)
        JPanel topPanel = new JPanel();
        topPanel.setLayout(new GridLayout(3, 1, 5, 5));
        topPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        // Pannello File Riga 1
        JPanel filePanel = new JPanel(new BorderLayout(10, 0));
        JButton selectFileBtn = new JButton("Seleziona File .mtx");
        selectFileBtn.setFont(new Font("Arial", Font.ITALIC, 14));
        fileField = new JTextField("Nessun file selezionato");
        fileField.setEditable(false);

        selectFileBtn.addActionListener(e -> chooseFile());

        filePanel.add(selectFileBtn, BorderLayout.WEST);
        filePanel.add(fileField, BorderLayout.CENTER);

        // Pannello Tolleranza Riga 2
        JPanel tolPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        JLabel tolLabel = new JLabel("Esponente Tolleranza (es. -4 per 10^-4): ");
        tolLabel.setFont(new Font("Arial", Font.BOLD, 14));
        tolInputField = new JTextField(5);
        tolInputField.setText("-4");

        tolPanel.add(tolLabel);
        tolPanel.add(tolInputField);

        // Pannello Run Riga 3
        JPanel runPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        runButton = new JButton("Avvia Calcolo");
        runButton.setFont(new Font("Arial", Font.BOLD, 16));
        runButton.setBackground(new Color(60, 130, 200));

        runButton.addActionListener(e -> startCalculations());

        runPanel.add(runButton);

        topPanel.add(filePanel);
        topPanel.add(tolPanel);
        topPanel.add(runPanel);

        // Pannello Output
        consoleOutput = new JTextArea();
        consoleOutput.setEditable(false);
        consoleOutput.setFont(new Font("Monospaced", Font.PLAIN, 12)); // Font fisso per allineare bene i risultati
        JScrollPane scrollPane = new JScrollPane(consoleOutput);
        scrollPane.setBorder(BorderFactory.createTitledBorder("Risultati dei Metodi Iterativi"));

        // Aggiungo i pannelli alla finestra principale
        add(topPanel, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);
    }

    private void chooseFile() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

        FileNameExtensionFilter filter = new FileNameExtensionFilter("Matrix Market Files (*.mtx)", "mtx");
        fileChooser.setFileFilter(filter);

        int result = fileChooser.showOpenDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
            selectedFile = fileChooser.getSelectedFile();
            fileField.setText(selectedFile.getAbsolutePath());
            fileField.setFont(new Font("Arial", Font.PLAIN, 12));
            log("File selezionato: " + selectedFile.getName());
        }
    }

    private void log(String message) {
        SwingUtilities.invokeLater(() -> {
            consoleOutput.append(message + "\n");
            consoleOutput.setCaretPosition(consoleOutput.getDocument().getLength());
        });
    }

    private void startCalculations() {
        if (selectedFile == null) {
            JOptionPane.showMessageDialog(this, "Si prega di selezionare un file .mtx", "Errore File",
                    JOptionPane.WARNING_MESSAGE);
            return;
        }

        int exp;
        try {
            exp = Integer.parseInt(tolInputField.getText().trim());
        } catch (NumberFormatException ex) {
            JOptionPane.showMessageDialog(this, "L'esponente della tolleranza deve essere un numero intero.",
                    "Errore Input", JOptionPane.ERROR_MESSAGE);
            return;
        }

        double tol = Math.pow(10, exp);

        runButton.setEnabled(false);
        consoleOutput.setText(""); 
        log("Configurata tolleranza: " + tol);
        log("Avvio elaborazione in corso, attendere...\n");

        new Thread(() -> {
            try {
                DMatrixSparseCSC matrix = ProjectMatrixUtils.importMatrix(selectedFile.getAbsolutePath());
                int n = matrix.numRows;
                log("Matrice caricata con successo (" + n + "x" + n + ")");

                SimpleMatrix exactSol = new SimpleMatrix(n, 1);
                exactSol.fill(1.0);
                SimpleMatrix b = SimpleMatrix.wrap(matrix).mult(exactSol);

                if (ProjectMatrixUtils.isPositiveDefinite(matrix)) {
                    log("La matrice è simmetrica e definita positiva. Avvio solutori...\n");

                    Jacobi jacobiSolver = new Jacobi();
                    MatrixResult jacobiResult = jacobiSolver.solve(matrix, b, tol, exactSol);
                    log(jacobiResult.toString());
                    log("------------------------------------------------------------");

                    GaussSeidel gsSolver = new GaussSeidel();
                    MatrixResult gsResult = gsSolver.solve(matrix, b, tol, exactSol);
                    log(gsResult.toString());
                    log("------------------------------------------------------------");

                    Gradient gradientSolver = new Gradient();
                    MatrixResult gradientResult = gradientSolver.solve(matrix, b, tol, exactSol);
                    log(gradientResult.toString());
                    log("------------------------------------------------------------");

                    ConjugateGradient conjugateGradientSolver = new ConjugateGradient("array");
                    MatrixResult conjugateGradientResult = conjugateGradientSolver.solve(matrix, b, tol, exactSol);
                    log(conjugateGradientResult.toString());
                    log("------------------------------------------------------------");

                    conjugateGradientSolver = new ConjugateGradient("hybrid");
                    conjugateGradientResult = conjugateGradientSolver.solve(matrix, b, tol, exactSol);
                    log(conjugateGradientResult.toString());
                    log("------------------------------------------------------------");

                    conjugateGradientSolver = new ConjugateGradient("ejml");
                    conjugateGradientResult = conjugateGradientSolver.solve(matrix, b, tol, exactSol);
                    log(conjugateGradientResult.toString());
                    log("------------------------------------------------------------");

                    log("\n--- ELABORAZIONE COMPLETATA ---");

                } else {
                    log("ATTENZIONE: La matrice non rispetta le condizioni (Simmetria/Positività).");
                }

            } catch (Exception e) {
                log("ERRORE DURANTE L'ESECUZIONE: " + e.getMessage());
                e.printStackTrace();
            } finally {
                SwingUtilities.invokeLater(() -> runButton.setEnabled(true));
            }
        }).start();
    }

    public static void main(String[] args) {
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            e.printStackTrace();
        }

        SwingUtilities.invokeLater(() -> {
            Main app = new Main();
            app.setVisible(true);
        });
    }
}