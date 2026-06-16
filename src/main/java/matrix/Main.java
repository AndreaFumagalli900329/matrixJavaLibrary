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
                long matrixImportStartTime = System.nanoTime();
                DMatrixSparseCSC matrix = ProjectMatrixUtils.importMatrix(selectedFile.getAbsolutePath());
                int n = matrix.numRows;
                long matrixImportTime = System.nanoTime() - matrixImportStartTime;
                log("Matrice caricata con successo (" + n + "x" + n + ") in " + String.format("%.2f ms", matrixImportTime / 1e6));

                SimpleMatrix exactSol = new SimpleMatrix(n, 1);
                exactSol.fill(1.0);
                SimpleMatrix b = SimpleMatrix.wrap(matrix).mult(exactSol);

                if (ProjectMatrixUtils.isPositiveDefinite(matrix)) {
                    log("La matrice è simmetrica e definita positiva. Avvio solutori...\n");

                    long jacobiStartTime = System.nanoTime();
                    Jacobi jacobiSolver = new Jacobi();
                    MatrixResult jacobiResult = jacobiSolver.solve(matrix, b, tol, exactSol);
                    double jacobiTime = (System.nanoTime() - jacobiStartTime) / 1e6;
                    log(jacobiResult.toString() + "\n" + String.format(" (Tempo totale: %.2f ms)", jacobiTime));
                    log("------------------------------------------------------------");

                    long gsStartTime = System.nanoTime();
                    GaussSeidel gsSolver = new GaussSeidel();
                    MatrixResult gsResult = gsSolver.solve(matrix, b, tol, exactSol);
                    double gsTime = (System.nanoTime() - gsStartTime) / 1e6;
                    log(gsResult.toString() + "\n" + String.format(" (Tempo totale: %.2f ms)", gsTime));
                    log("------------------------------------------------------------");

                    long gradientStartTime = System.nanoTime();
                    Gradient gradientSolver = new Gradient();
                    MatrixResult gradientResult = gradientSolver.solve(matrix, b, tol, exactSol);
                    double gradientTime = (System.nanoTime() - gradientStartTime) / 1e6;
                    log(gradientResult.toString() + "\n" + String.format(" (Tempo totale: %.2f ms)", gradientTime));
                    log("------------------------------------------------------------");

                    long cgHybridStartTime = System.nanoTime();
                    ConjugateGradient conjugateGradientSolver = new ConjugateGradient();
                    MatrixResult conjugateGradientResult = conjugateGradientSolver.solve(matrix, b, tol, exactSol);
                    double cgHybridTime = (System.nanoTime() - cgHybridStartTime) / 1e6;
                    log(conjugateGradientResult.toString() + "\n" + String.format(" (Tempo totale: %.2f ms)", cgHybridTime));
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