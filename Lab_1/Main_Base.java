import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

public class Main {
    public static void main(String[] args) {
//        generateValues();
        double[] B = readVectorFromFile("B.txt");
        double[][] MC = readMatrixFromFile("MC.txt");
        double[] D = readVectorFromFile("D.txt");
        double[][] MM = readMatrixFromFile("MM.txt");

        double[][] MZ = readMatrixFromFile("MZ.txt");
        double a = readValueFromFile("a.txt");

        int amount = 3;
        System.out.println("Vector B:");
        output(B);

        System.out.println("\nMatrix MC:");
        output(MC, amount);

        System.out.println("\nVector D:");
        output(D);

        System.out.println("\nMatrix MM:");
        output(MM, amount);

        System.out.println("\nMatrix MZ:");
        output(MZ, amount);

        System.out.print("\nSkalar a = ");
        output(a);


        /*double[] arr = new double[10];
//        Arrays.fill(arr, 0.1);
//        double sum = 0;
//        for (double v : arr) {
//            sum += v;
//        }
//        System.out.println("Sum for array of 10 elements with value 0.1 (should be 0.1*10=1):");
//        System.out.println("Sum = " + sum);
//        System.out.println("Kahan sum = "+kahanSum(arr));
//        System.out.println("Kahan-Babushka sum = "+kahanBabushkaSum(arr));
//        System.out.println("Sum for [0.1; 10^100; 1; -10^100] (Should be 1.1):");
//        System.out.println("Sum = "+(0.1+Math.pow(10,100)+1-Math.pow(10,100)));
//        System.out.println("Kahan sum = "+kahanSum(0.1,Math.pow(10,100),1,-Math.pow(10,100)));
//        System.out.println("Kahan-Babushka sum = "+kahanBabushkaSum(0.1,Math.pow(10,100),1,-Math.pow(10,100)));*/

        double timeVectorSimple, timeVectorKahan, timeVectorKahanBabushka;
        double timeMatrixSimple, timeMatrixKahan, timeMatrixKahanBabushka;

        //simple method
        System.out.println("\n\nVector C:");
        timeVectorSimple = System.currentTimeMillis();
        double[] C = subtract(multiply(MC, B), multiply(MM, D));
        timeVectorSimple = System.currentTimeMillis()- timeVectorSimple;
        output(C);
        writeToFile("C_simple.txt", C);
        System.out.println("\nMatrix MF:");
        timeMatrixSimple = System.currentTimeMillis();
        double[][] MF = add(multiply(multiply(MC, MZ), min(add(B, D))), multiply(multiply(MM, add(MC, MM)), a));
        timeMatrixSimple = System.currentTimeMillis()-timeMatrixSimple;
        output(MF, amount);
        writeToFile("MF_simple.txt", MF);


        //kahan method
        System.out.println("\n\nVector C_Kahan:");
        timeVectorKahan = System.currentTimeMillis();
        double[] C_Kahan = subtractKahan(multiplyKahan(MC, B), multiplyKahan(MM, D));
        timeVectorKahan= System.currentTimeMillis()- timeVectorKahan;
        output(C_Kahan);
        writeToFile("C_kahan.txt", C_Kahan);
        System.out.println("\nMatrix MF_Kahan:");
        timeMatrixKahan = System.currentTimeMillis();
        double[][] MF_Kahan = addKahan(multiply(multiplyKahan(MC, MZ), min(addKahan(B, D))), multiply(multiplyKahan(MM, addKahan(MC, MM)), a));
        timeMatrixKahan = System.currentTimeMillis()- timeMatrixKahan;
        output(MF_Kahan, amount);
        writeToFile("MF_kahan.txt", MF_Kahan);


        //kahan-babushka
        System.out.println("\n\nVector C_Kahan-Babushka:");
        timeVectorKahanBabushka= System.currentTimeMillis();
        double[] C_KahanBabushka = subtractKahanBabushka(multiplyKahanBabushka(MC, B), multiplyKahanBabushka(MM, D));
        timeVectorKahanBabushka= System.currentTimeMillis()- timeVectorKahanBabushka;
        output(C_KahanBabushka);
        writeToFile("C_kahan_babushka.txt", C_KahanBabushka);
        System.out.println("\nMatrix MF_Kahan-Babushka:");
        timeMatrixKahanBabushka = System.currentTimeMillis();
        double[][] MF_KahanBabushka = addKahanBabushka(multiply(multiplyKahanBabushka(MC, MZ), min(addKahanBabushka(B, D))), multiply(multiplyKahanBabushka(MM, addKahanBabushka(MC, MM)), a));
        timeMatrixKahanBabushka = System.currentTimeMillis()- timeMatrixKahanBabushka;
        output(MF_KahanBabushka, amount);
        writeToFile("MF_kahan_babushka.txt", MF_KahanBabushka);

        System.out.println("Time vector simple = "+timeVectorSimple);
        System.out.println("Time matrix simple = "+timeMatrixSimple);
        System.out.println("Time vector Kahan = "+timeVectorKahan);
        System.out.println("Time matrix Kahan = "+timeMatrixKahan);
        System.out.println("Time vector Kahan-Babushka = "+timeVectorKahanBabushka);
        System.out.println("Time matrix Kahan-Babushka = "+timeMatrixKahanBabushka);
        //BigDecimal type for matrix 10x10
//        BigDecimal[] B_BigDecimal = readBigDecimalVectorFromFile("B.txt");
//        BigDecimal[][] MC_BigDecimal = readBigDecimalMatrixFromFile("MC.txt");
//        BigDecimal[] D_BigDecimal = readBigDecimalVectorFromFile("D.txt");
//        BigDecimal[][] MM_BigDecimal = readBigDecimalMatrixFromFile("MM.txt");
//
//        BigDecimal[][] MZ_BigDecimal = readBigDecimalMatrixFromFile("MZ.txt");
//        BigDecimal a_BigDecimal = readBigDecimalValueFromFile("a.txt");
//
//        System.out.println("\n\nVector C_BigDecimal:");
//        BigDecimal[] C_BigDecimal = subtractBigDecimal(multiplyBigDecimal(MC_BigDecimal, B_BigDecimal), multiplyBigDecimal(MM_BigDecimal, D_BigDecimal));
//        outputBigDecimal(C_BigDecimal);
//        writeToFile("C_BigDecimal.txt", C_BigDecimal);
//        System.out.println("\nMatrix MF_BigDecimal:");
//        BigDecimal[][] MF_BigDecimal = addBigDecimal(multiplyBigDecimal(multiplyBigDecimal(MC_BigDecimal, MZ_BigDecimal), minBigDecimal(addBigDecimal(B_BigDecimal, D_BigDecimal))), multiplyBigDecimal(multiplyBigDecimal(MM_BigDecimal, addBigDecimal(MC_BigDecimal, MM_BigDecimal)), a_BigDecimal));
//        outputBigDecimal(MF_BigDecimal);
//        writeToFile("MF_BigDecimal.txt", MF_BigDecimal);
    }

    public static void generateValues() {
        int dimension = 1000;
        int magnitudeMatrix = 150;
        int magnitudeVector = -150;
        double[] B = generate(dimension, magnitudeVector);
        writeToFile("B.txt", B);

        double[][] MC = generate(dimension, dimension, magnitudeMatrix);
        writeToFile("MC.txt", MC);

        double[] D = generate(dimension, magnitudeVector);
        writeToFile("D.txt", D);

        double[][] MM = generate(dimension, dimension, magnitudeMatrix);
        writeToFile("MM.txt", MM);

        double[][] MZ = generate(dimension, dimension, magnitudeMatrix);
        writeToFile("MZ.txt", MZ);

        double a = generate();
        writeToFile("a.txt", a);
    }

    /**----------------------------Simple Method----------------------------------**/
    public static double[] add(double[] a, double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    public static double[][] add(double[][] a, double[][] b) {
        double[][] result = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                result[i][j] = a[i][j] + b[i][j];
            }
        }
        return result;
    }

    public static double[] subtract(double[] a, double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    public static double[] multiply(double[][] matrix, double[] vector) {
        double[] result = new double[matrix.length];
        double sum;
        for (int i = 0; i < matrix.length; i++) {
            sum = 0;
            for (int j = 0; j < matrix[i].length; j++) {
                sum += matrix[i][j] * vector[j];
            }
            result[i] = sum;
        }
        return result;
    }

    public static double[][] multiply(double[][] a, double[][] b) {
        double[][] result = new double[a.length][b[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b[0].length; j++) {
                double sum = 0;
                for (int k = 0; k < a[i].length; k++) {
                    sum += a[i][k] * b[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    public static double[][] multiply(double[][] a, double b) {
        double[][] result = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                result[i][j] = a[i][j] * b;
            }
        }
        return result;
    }

    public static double min(double[] vector) {
        double min = vector[0];
        for (double v : vector) {
            if (v < min) {
                min = v;
            }
        }
        return min;
    }


    /**-------------------------------Kahan Method--------------------------------**/
    private static double kahanSum(double... fa) {
        double sum = 0.0;
        double c = 0.0;

        for (double f : fa) {
            double y = f - c;
            double t = sum + y;
            double z = t-sum;
            c = z - y;

            sum = t;
        }

        return sum;
    }

    public static double[] addKahan(double[] a, double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = kahanSum(a[i], b[i]);
        }
        return result;
    }

    public static double[][] addKahan(double[][] a, double[][] b) {
        double[][] result = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                result[i][j] = kahanSum(a[i][j], b[i][j]);
            }
        }
        return result;
    }

    public static double[] subtractKahan(double[] a, double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = kahanSum(a[i], -1*b[i]);
        }
        return result;
    }

    public static double[] multiplyKahan(double[][] matrix, double[] vector) {
        double[] result = new double[matrix.length];
        double sum, c, y, t;
        for(int i = 0; i <matrix.length; i++) {
            sum = 0.0;
            c = 0.0;

            for (int j = 0; j < vector.length; j++) {
                y = vector[j] * matrix[i][j] - c;
                t = sum + y;
                c = (t-sum)-y;

                sum = t;
            }
            result[i] = sum+c;
        }
        return result;
    }

    public static double[][] multiplyKahan(double[][] a, double[][] b) {
        double[][] result = new double[a.length][b[0].length];
        double sum, c, y, t;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b[0].length; j++) {
                sum = 0;
                c = 0;
                for (int k = 0; k < a[i].length; k++) {
                    y = a[i][k] * b[k][j] - c;
                    t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    /**----------------------------Kahan-Babushka Method-------------------------**/
    private static double kahanBabushkaSum(double... fa) {
        double sum = 0.0;
        double ccs = 0.0;
        double cs = 0.0;
        double cc, c, t;

        for (double f : fa) {
            t = sum + f;
            if (Math.abs(sum) >= Math.abs(f)) {
                c = (sum - t) + f;
            } else {
                c = (f - t) + sum;
            }
            sum = t;
            t = cs + c;
            if (Math.abs(cs) >= Math.abs(c)) {
                cc = (cs - t) + c;
            } else {
                cc = (c - t) + cs;
            }
            cs = t;
            ccs = ccs + cc;
        }
        return sum + cs + ccs;
    }

    public static double[] addKahanBabushka(double[] a, double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = kahanBabushkaSum(a[i], b[i]);
        }
        return result;
    }

    public static double[][] addKahanBabushka(double[][] a, double[][] b) {
        double[][] result = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                result[i][j] = kahanBabushkaSum(a[i][j], b[i][j]);
            }
        }
        return result;
    }

    public static double[] subtractKahanBabushka(double[] a, double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = kahanBabushkaSum(a[i], -1*b[i]);
        }
        return result;
    }

    public static double[] multiplyKahanBabushka(double[][] matrix, double[] vector) {
        double[] result = new double[matrix.length];
        double sum, ccs, cs, cc, c, t;
        for (int i = 0; i < matrix.length; i++) {
            sum = 0.0;
            ccs = 0.0;
            cs = 0.0;

            for (int j = 0; j < vector.length; j++) {
                t = sum + matrix[i][j] * vector[j];
                if (Math.abs(sum) >= Math.abs(matrix[i][j] * vector[j])) {
                    c = (sum - t) + matrix[i][j] * vector[j];
                } else {
                    c = (matrix[i][j] * vector[j] - t) + sum;
                }
                sum = t;
                t = cs + c;
                if (Math.abs(cs) >= Math.abs(c)) {
                    cc = (cs - t) + c;
                } else {
                    cc = (c - t) + cs;
                }
                cs = t;
                ccs = ccs + cc;
            }
            result[i] = sum + cs + ccs;
        }
        return result;
    }

    public static double[][] multiplyKahanBabushka(double[][] a, double[][] b) {
        double[][] result = new double[a.length][b[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b[0].length; j++) {
                double sum = 0;
                double ccs = 0.0;
                double cs = 0.0;
                double cc, c, t = 0, product = 0;
                for (int k = 0; k<a[i].length;k++) {
                    product = a[i][k] * b[k][j];
                    t = sum + product;
                    if (Math.abs(sum) >= Math.abs(product)) {
                        c = (sum - t) + product;
                    } else {
                        c = (product - t) + sum;
                    }
                    sum = t;
                    t = cs + c;
                    if (Math.abs(cs) >= Math.abs(c)) {
                        cc = (cs - t) + c;
                    } else {
                        cc = (c - t) + cs;
                    }
                    cs = t;
                    ccs = ccs + cc;
                }
                result[i][j] = sum + cs + ccs;
            }
        }
        return result;
    }


    /**--------------------------Big Decimal type----------------------------**/
    public static BigDecimal[] readBigDecimalVectorFromFile(String filename) {
        List<BigDecimal> list = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = reader.readLine()) != null) {
                list.add(new BigDecimal(line));
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
        return list.toArray(new BigDecimal[0]);
    }

    public static BigDecimal[][] readBigDecimalMatrixFromFile(String filename) {
        List<BigDecimal[]> list = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] values = line.split(" ");
                BigDecimal[] row = new BigDecimal[values.length];
                for (int i = 0; i < values.length; i++) {
                    row[i] = new BigDecimal(values[i]);
                }
                list.add(row);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
        return list.toArray(new BigDecimal[0][]);
    }

    public static BigDecimal readBigDecimalValueFromFile(String filename) {
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            return new BigDecimal(reader.readLine());
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
        return BigDecimal.ZERO;
    }

    public static void outputBigDecimal(BigDecimal[] vector) {
        for (BigDecimal v : vector) {
            System.out.println(v);
        }
    }

    public static void outputBigDecimal(BigDecimal[][] matrix) {
        int numCols = matrix[0].length;

        int[] maxColumnWidths = maxWidthInCol(matrix);
        for (BigDecimal[] decimals : matrix) {
            for (int col = 0; col < numCols; col++) {
                System.out.printf("%-" + (maxColumnWidths[col] + 2) + "s", decimals[col]);
            }
            System.out.println();
        }
    }

    public static int[] maxWidthInCol(BigDecimal[][] matrix) {
        int numCols = matrix[0].length;

        int[] maxColumnWidths = new int[numCols];
        int maxColumnWidth;
        int width;
        for (int col = 0; col < numCols; col++) {
            maxColumnWidth = 0;
            for (BigDecimal[] decimals : matrix) {
                width = String.valueOf(decimals[col]).length();
                if (width > maxColumnWidth) {
                    maxColumnWidth = width;
                }
            }
            maxColumnWidths[col] = maxColumnWidth;
        }
        return maxColumnWidths;
    }

    public static BigDecimal[] addBigDecimal(BigDecimal[] a, BigDecimal[] b) {
        BigDecimal[] result = new BigDecimal[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i].add(b[i]);
        }
        return result;
    }

    public static BigDecimal[][] addBigDecimal(BigDecimal[][] a, BigDecimal[][] b) {
        BigDecimal[][] result = new BigDecimal[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                result[i][j] = a[i][j].add(b[i][j]);
            }
        }
        return result;
    }

    public static BigDecimal[] subtractBigDecimal(BigDecimal[] a, BigDecimal[] b) {
        BigDecimal[] result = new BigDecimal[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i].subtract(b[i]);
        }
        return result;
    }

    public static BigDecimal[] multiplyBigDecimal(BigDecimal[][] matrix, BigDecimal[] vector) {
        BigDecimal[] result = new BigDecimal[matrix.length];
        BigDecimal sum;
        for (int i = 0; i < matrix.length; i++) {
            sum = BigDecimal.ZERO;
            for (int j = 0; j < matrix[i].length; j++) {
                sum = sum.add(matrix[i][j].multiply(vector[j]));
            }
            result[i] = sum;
        }
        return result;
    }

    public static BigDecimal[][] multiplyBigDecimal(BigDecimal[][] a, BigDecimal[][] b) {
        BigDecimal[][] result = new BigDecimal[a.length][b[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b[0].length; j++) {
                BigDecimal sum = BigDecimal.ZERO;
                for (int k = 0; k < a[i].length; k++) {
                    sum = sum.add(a[i][k].multiply(b[k][j]));
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    public static BigDecimal[][] multiplyBigDecimal(BigDecimal[][] a, BigDecimal b) {
        BigDecimal[][] result = new BigDecimal[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                result[i][j] = a[i][j].multiply(b);
            }
        }
        return result;
    }

    public static BigDecimal minBigDecimal(BigDecimal[] vector) {
        BigDecimal min = vector[0];
        for (BigDecimal v : vector) {
            if (v.compareTo(min) < 0) {
                min = v;
            }
        }
        return min;
    }

    public static void writeToFile(String filename, BigDecimal[] vector) {
        try (FileWriter writer = new FileWriter(filename)) {
            StringBuilder content = new StringBuilder();
            for (BigDecimal v : vector) {
                content.append(v).append("\n");
            }
            writer.write(content.toString());
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void writeToFile(String filename, BigDecimal[][] matrix) {
        try (FileWriter writer = new FileWriter(filename)) {
            StringBuilder content = new StringBuilder();
            for (BigDecimal[] row : matrix) {
                for (BigDecimal v : row) {
                    content.append(v).append(" ");
                }
                content.append("\n");
            }
            writer.write(content.toString());
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    /**--------------------------Standard Methods----------------------------**/

    public static void output(double a) {
        System.out.println(a);
    }

    public static void output(double[] vector) {
        for (double v : vector) {
            System.out.println(v);
        }
    }

    public static void output(double[][] matrix) {
        int numCols = matrix[0].length;

        int[] maxColumnWidths = maxWidthInCol(matrix);
        for (double[] doubles : matrix) {
            for (int col = 0; col < numCols; col++) {
                System.out.printf("%-" + (maxColumnWidths[col] + 2) + "s", doubles[col]);
            }
            System.out.println();
        }
    }

    public static void output(double[][] matrix, int amount) {
        if(matrix.length<=amount+1) {
            output(matrix);
            return;
        }

        double[][] matrixDemo = cutMatrix(matrix, amount);
        int[] maxColumnWidths = maxWidthInCol(matrixDemo);

        for (int row = 0; row < amount+1; row++) {
            for (int col = 0; col < amount+1; col++) {
                if (row + col == amount) {
                    System.out.printf("%-" + (maxColumnWidths[col] + 2) + "s", "         ∙∙∙");
                }
                else {
                    System.out.printf("%-" + (maxColumnWidths[col] + 2) + "s", matrixDemo[row][col]);
                }
            }
            System.out.println();
        }
    }

    public static int[] maxWidthInCol(double[][] matrix) {
        int numCols = matrix[0].length;

        int[] maxColumnWidths = new int[numCols];
        int maxColumnWidth;
        int width;
        for (int col = 0; col < numCols; col++) {
            maxColumnWidth = 0;
            for (double[] doubles : matrix) {
                width = String.valueOf(doubles[col]).length();
                if (width > maxColumnWidth) {
                    maxColumnWidth = width;
                }
            }
            maxColumnWidths[col] = maxColumnWidth;
        }
        return maxColumnWidths;
    }

    public static double[][] cutMatrix(double[][] matrix, int amount) {
        int numRows = matrix.length;
        int numCols = matrix[0].length;
        int n = amount + 1;

        double[][] matrixDemo = new double[n][n];
        for (int row = 0; row < amount; row++) {
            for (int col = 0; col < n; col++) {
                if (row + col < amount) {
                    matrixDemo[row][col] = matrix[row][col];
                } else {
                    matrixDemo[row+1][col] = matrix[numRows - amount + row][numCols - n + col];
                }
            }
        }
        return matrixDemo;
    }

    public static double[] readVectorFromFile(String filename) {
        List<Double> list = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = reader.readLine()) != null) {
                list.add(Double.parseDouble(line));
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
        return list.stream().mapToDouble(d -> d).toArray();
    }

    public static double[][] readMatrixFromFile(String filename) {
        List<double[]> list = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] values = line.split(" ");
                double[] row = new double[values.length];
                for (int i = 0; i < values.length; i++) {
                    row[i] = Double.parseDouble(values[i]);
                }
                list.add(row);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
        return list.toArray(new double[0][]);
    }

    public static double readValueFromFile(String filename) {
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            return Double.parseDouble(reader.readLine());
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
        return 0;
    }

    public static void writeToFile(String filename, double[] vector) {
        try (FileWriter writer = new FileWriter(filename)) {
            StringBuilder content = new StringBuilder();
            for (double v : vector) {
                content.append(v).append("\n");
            }
            writer.write(content.toString());
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void writeToFile(String filename, double[][] matrix) {
        try (FileWriter writer = new FileWriter(filename)) {
            StringBuilder content = new StringBuilder();
            for (double[] row : matrix) {
                for (double v : row) {
                    content.append(v).append(" ");
                }
                content.append("\n");
            }
            writer.write(content.toString());
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void writeToFile(String filename, double value) {
        try (FileWriter writer = new FileWriter(filename)) {
            writer.write(String.valueOf(value));
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static double generate() {
        return 9*Math.random()+1;
    }

    public static double[] generate(int size, int magnitude) {
        double[] vector = new double[size];
        for (int i = 0; i < size; i++) {
            // Math.random()*(max-min)+min
            // max-min = Math.pow(10, magnitude+1)-Math.pow(10, magnitude) = 9*Math.pow(10, magnitude)
            vector[i] = Math.random()*(9*Math.pow(10, magnitude))+Math.pow(10, magnitude);
        }
        return vector;
    }

    public static double[][] generate(int rows, int cols, int magnitude) {
        double[][] matrix = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            matrix[i] = generate(cols, magnitude);
        }
        return matrix;
    }
}
