import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReentrantLock;

public class Main {
    private static final ReentrantLock lock = new ReentrantLock();
    private static final ExecutorService executor = Executors.newFixedThreadPool(12);
    private static final AtomicLong timeVectorSimple = new AtomicLong();
    private static final AtomicLong timeVectorKahan = new AtomicLong();
    private static final AtomicLong timeVectorKahanBabushka = new AtomicLong();
    private static final AtomicLong timeMatrixSimple = new AtomicLong();
    private static final AtomicLong timeMatrixKahan = new AtomicLong();
    private static final AtomicLong timeMatrixKahanBabushka = new AtomicLong();

    public static void main(String[] args) {
        String path = "C:\\Users\\Lenovo\\Documents\\KPI\\III курс\\ПРО\\Data\\";
        double[] B = readVectorFromFile(path+"B.txt");
        double[][] MC = readMatrixFromFile(path+"MC.txt");
        double[] D = readVectorFromFile(path+"D.txt");
        double[][] MM = readMatrixFromFile(path+"MM.txt");

        double[][] MZ = readMatrixFromFile(path+"MZ.txt");
        double a = readValueFromFile(path+"a.txt");

        int amount = 3;
        System.out.println("Vector B:");
        output(B, amount);

        System.out.println("\nMatrix MC:");
        output(MC, amount);

        System.out.println("\nVector D:");
        output(D, amount);

        System.out.println("\nMatrix MM:");
        output(MM, amount);

        System.out.println("\nMatrix MZ:");
        output(MZ, amount);

        System.out.print("\nSkalar a = ");
        output(a);

        long timeS = System.currentTimeMillis();
        Callable<double[]> taskVectorSimple = () -> {
            long startTime = System.currentTimeMillis();
            double[] C = subtract(multiply(MC, B), multiply(MM, D)); // Assume computeVector() is a method that computes something
            timeVectorSimple.set(System.currentTimeMillis() - startTime);
            lock.lock();
            try {
                System.out.println("\n\nVector C with modified threads:");
                output(C, amount);
                writeToFile("C_simple.txt", C);
            } finally {
                lock.unlock();
            }
            return C;
        };

        Callable<double[][]> taskMatrixSimple = () -> {
            long startTime = System.currentTimeMillis();
            Future<double[][]> future1 = executor.submit(() -> multiply(multiply(MC, MZ), min(add(B, D))));
            Future<double[][]> future2 = executor.submit(() -> multiply(multiply(MM, add(MC, MM)), a));

            double[][] MF = add(future1.get(), future2.get());
            timeMatrixSimple.set(System.currentTimeMillis() - startTime);
            lock.lock();
            try {
                System.out.println("\nMatrix MF with modified threads:");
                output(MF, amount);
                writeToFile("MF_simple.txt", MF);
            } finally {
                lock.unlock();
            }
            return MF;
        };

        Callable<double[]> taskVectorKahan = () -> {
            long startTime = System.currentTimeMillis();
            double[] C_Kahan = subtract(multiplyKahan(MC, B), multiplyKahan(MM, D));
            timeVectorKahan.set(System.currentTimeMillis() - startTime);
            lock.lock();
            try {
                System.out.println("\n\nVector C_Kahan with modified threads:");
                output(C_Kahan, amount);
                writeToFile("C_kahan.txt", C_Kahan);
            } finally {
                lock.unlock();
            }
            return C_Kahan;
        };

        Callable<double[][]> taskMatrixKahan = () -> {
            long startTime = System.currentTimeMillis();
            Future<double[][]> future1 = executor.submit(() -> multiply(multiplyKahan(MC, MZ), min(add(B, D))));
            Future<double[][]> future2 = executor.submit(() -> multiply(multiplyKahan(MM, add(MC, MM)), a));

            double[][] MF_Kahan = add(future1.get(), future2.get());
            timeMatrixKahan.set(System.currentTimeMillis() - startTime);
            lock.lock();
            try {
                System.out.println("\nMatrix MF_Kahan with modified threads:");
                output(MF_Kahan, amount);
                writeToFile("MF_kahan.txt", MF_Kahan);
            } finally {
                lock.unlock();
            }
            return MF_Kahan;
        };

        Callable<double[]> taskVectorKahanBabushka = () -> {
            long startTime = System.currentTimeMillis();
            double[] C_KahanBabushka = subtract(multiplyKahanBabushka(MC, B), multiplyKahanBabushka(MM, D));
            timeVectorKahanBabushka.set(System.currentTimeMillis() - startTime);
            lock.lock();
            try {
                System.out.println("\n\nVector C_Kahan-Babushka with threads:");
                output(C_KahanBabushka, amount);
                writeToFile("C_kahan_babushka.txt", C_KahanBabushka);
            } finally {
                lock.unlock();
            }
            return C_KahanBabushka;
        };

        Callable<double[][]> taskMatrixKahanBabushka = () -> {
            long startTime = System.currentTimeMillis();
            Future<double[][]> future1 = executor.submit(() -> multiply(multiplyKahanBabushka(MC, MZ), min(add(B, D))));
            Future<double[][]> future2 = executor.submit(() -> multiply(multiplyKahanBabushka(MM, add(MC, MM)), a));

            double[][] MF_KahanBabushka = add(future1.get(), future2.get());
            timeMatrixKahanBabushka.set(System.currentTimeMillis() - startTime);
            lock.lock();
            try {
                System.out.println("\nMatrix MF_Kahan-Babushka with threads:");
                output(MF_KahanBabushka, amount);
                writeToFile("MF_kahan_babushka.txt", MF_KahanBabushka);
            } finally {
                lock.unlock();
            }
            return MF_KahanBabushka;
        };

        Future<double[]> futureVectorSimple = executor.submit(taskVectorSimple);
        Future<double[][]> futureMatrixSimple = executor.submit(taskMatrixSimple);
        Future<double[]> futureVectorKahan = executor.submit(taskVectorKahan);
        Future<double[][]> futureMatrixKahan = executor.submit(taskMatrixKahan);
        Future<double[]> futureVectorKahanBabushka = executor.submit(taskVectorKahanBabushka);
        Future<double[][]> futureMatrixKahanBabushka = executor.submit(taskMatrixKahanBabushka);
        try {
            futureVectorSimple.get();
            futureMatrixSimple.get();
            futureVectorKahan.get();
            futureMatrixKahan.get();
            futureVectorKahanBabushka.get();
            futureMatrixKahanBabushka.get();
            executor.shutdown();
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }

        /*Runnable taskVectorSimple = createSimpleVectorTask(B, MC, D, MM, amount);
        Runnable taskMatrixSimple = createSimpleMatrixTask(MC, MZ, B, D, MM, a, amount);
        Runnable taskVectorKahan = createKahanVectorTask(B, MC, D, MM, amount);
        Runnable taskMatrixKahan = createKahanMatrixTask(MC, MZ, B, D, MM, a, amount);
        Runnable taskVectorKahanBabushka = createKahanBabushkaVectorTask(B, MC, D, MM, amount);
        Runnable taskMatrixKahanBabushka = createKahanBabushkaMatrixTask(MC, MZ, B, D, MM, a, amount);

        scheduler.schedule(taskVectorSimple, 0, TimeUnit.MILLISECONDS);
        scheduler.schedule(taskMatrixSimple, 0, TimeUnit.MILLISECONDS);
        scheduler.schedule(taskVectorKahan, 0, TimeUnit.MILLISECONDS);
        scheduler.schedule(taskMatrixKahan, 0, TimeUnit.MILLISECONDS);
        scheduler.schedule(taskVectorKahanBabushka, 0, TimeUnit.MILLISECONDS);
        scheduler.schedule(taskMatrixKahanBabushka, 0, TimeUnit.MILLISECONDS);
        scheduler.shutdown();
        try {
            if (!scheduler.awaitTermination(60, TimeUnit.SECONDS)) {
                scheduler.shutdownNow();
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }*/
        /*executor.execute(taskVectorSimple);
        executor.execute(taskMatrixSimple);
        executor.execute(taskVectorKahan);
        executor.execute(taskMatrixKahan);
        executor.execute(taskVectorKahanBabushka);
        executor.execute(taskMatrixKahanBabushka);

        executor.shutdown();
        try {
            if (!executor.awaitTermination(60, TimeUnit.SECONDS)) {
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }*/

        //simple method with modified treads
        /*Thread threadVectorSimpleModified = new Thread(() -> {
            double[] C;
            try {
                long startTime = System.currentTimeMillis();
                C = subtract(multiply(MC, B), multiply(MM, D));
                timeVectorSimple.set(System.currentTimeMillis() - startTime);
            } finally {
                latch.countDown();
            }

            try {
                outputSemaphore.acquire();
                System.out.println("\n\nVector C with modified threads:");
                output(C, amount);
                writeToFile("C_simple.txt", C);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            } finally {
                outputSemaphore.release();
            }
        });*/

        /*Thread threadMatrixSimpleModified = new Thread(() -> {
            final CyclicBarrier barrier = new CyclicBarrier(3);
//                double[][] MF = add(multiply(multiply(MC, MZ), min(add(B, D))), multiply(multiply(MM, add(MC, MM)), a));
            double[][][] matrix = new double[2][][];
            double[][] MF;
            try {
                long startTime = System.currentTimeMillis();
                Thread tread1 = new Thread(() -> {
                    try {
                        matrix[0] = multiply(multiply(MC, MZ), min(add(B, D)));
                        barrier.await();
                    } catch (Exception e) {
                        Thread.currentThread().interrupt();
                        e.printStackTrace();
                    }
                });// matrix[0] = multiply(multiply(MC, MZ), min(add(B, D))
                Thread tread2 = new Thread(() -> {
                    try {
                        matrix[1] = multiply(multiply(MM, add(MC, MM)), a);
                        barrier.await();
                    } catch (Exception e) {
                        Thread.currentThread().interrupt();
                        e.printStackTrace();
                    }
                });// matrix[1] = multiply(multiply(MM, add(MC, MM)), a)
                tread1.start();
                tread2.start();

                try {
                    barrier.await();
                    MF = add(matrix[0], matrix[1]);
                    timeMatrixSimple.set(System.currentTimeMillis() - startTime);
                    try {
                        outputSemaphore.acquire();
                        System.out.println("\nMatrix MF with modified threads:");
                        output(MF, amount);
                        writeToFile("MF_simple.txt", MF);
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                    } finally {
                        outputSemaphore.release();
                    }
                } catch (Exception e) {
                    Thread.currentThread().interrupt();
                    e.printStackTrace();
                }
            } finally {
                latch.countDown();
            }
        });*/


        //kahan method with modified threads
        /*Thread threadVectorKahanModified = new Thread(() -> {
//                double[] C_Kahan = subtractKahan(multiplyKahan(MC, B), multiplyKahan(MM, D));
            double[] C_Kahan;
            try {
                long startTime = System.currentTimeMillis();
                C_Kahan = subtract(multiplyKahan(MC, B), multiplyKahan(MM, D));
                timeVectorKahan.set(System.currentTimeMillis() - startTime);
            } finally {
                latch.countDown();
            }

            try {
                outputSemaphore.acquire();
                System.out.println("\n\nVector C_Kahan with modified threads:");
                output(C_Kahan, amount);
                writeToFile("C_kahan.txt", C_Kahan);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            } finally {
                outputSemaphore.release();
            }
        });*/

        /*Thread threadMatrixKahanModified = new Thread(() -> {
            final CyclicBarrier barrier = new CyclicBarrier(3);
            double[][][] matrix = new double[2][][];
            double[][] MF_Kahan;
            try {
                long startTime = System.currentTimeMillis();
//                double[][] MF_Kahan = addKahan(multiply(multiplyKahan(MC, MZ), min(addKahan(B, D))), multiply(multiplyKahan(MM, addKahan(MC, MM)), a));
                Thread tread1 = new Thread(() -> {
                    try {
                        matrix[0] = multiply(multiplyKahan(MC, MZ), min(add(B, D)));
                        barrier.await();
                    } catch (Exception e) {
                        Thread.currentThread().interrupt();
                        e.printStackTrace();
                    }
                });// matrix[0] = multiply(multiplyKahan(MC, MZ), min(addKahan(B, D))
                Thread tread2 = new Thread(() -> {
                    try {
                        matrix[1] = multiply(multiplyKahan(MM, add(MC, MM)), a);
                        barrier.await();
                    } catch (Exception e) {
                        Thread.currentThread().interrupt();
                        e.printStackTrace();
                    }
                });// matrix[1] = multiply(multiplyKahan(MM, addKahan(MC, MM)), a)
                tread1.start();
                tread2.start();
                try {
                    barrier.await();
                    MF_Kahan = add(matrix[0], matrix[1]);
                    timeMatrixKahan.set(System.currentTimeMillis() - startTime);
                    try {
                        outputSemaphore.acquire();
                        System.out.println("\nMatrix MF_Kahan with modified threads:");
                        output(MF_Kahan, amount);
                        writeToFile("MF_kahan.txt", MF_Kahan);
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                    } finally {
                        outputSemaphore.release();
                    }
                } catch (Exception e) {
                    Thread.currentThread().interrupt();
                    e.printStackTrace();
                }
            } finally {
                latch.countDown();
            }
        });*/


        //kahan-babushka method with modified threads
        /*Thread threadVectorKahanBabushkaModified = new Thread(() -> {
            double[] C_KahanBabushka;
            try {
                long startTime = System.currentTimeMillis();
                C_KahanBabushka = subtract(multiplyKahanBabushka(MC, B), multiplyKahanBabushka(MM, D));
                timeVectorKahanBabushka.set(System.currentTimeMillis() - startTime);
            } finally {
                latch.countDown();
            }

            try {
                outputSemaphore.acquire();
                System.out.println("\n\nVector C_Kahan-Babushka with threads:");
                output(C_KahanBabushka, amount);
                writeToFile("C_kahan_babushka.txt", C_KahanBabushka);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            } finally {
                outputSemaphore.release();
            }
        });*/

        /*Thread threadMatrixKahanBabushkaModified = new Thread(() -> {
            final CyclicBarrier barrier = new CyclicBarrier(3);
            double[][][] matrix = new double[2][][];
            double[][] MF_KahanBabushka;
            try {
                long startTime = System.currentTimeMillis();
                Thread tread1 = new Thread(() -> {
                    try {
                        matrix[0] = multiply(multiplyKahanBabushka(MC, MZ), min(add(B, D)));
                        barrier.await();
                    } catch (Exception e) {
                        Thread.currentThread().interrupt();
                        e.printStackTrace();
                    }
                });// matrix[0] = multiply(multiplyKahanBabushka(MC, MZ), min(addKahanBabushka(B, D))
                Thread tread2 = new Thread(() -> {
                    try {
                        matrix[1] = multiply(multiplyKahanBabushka(MM, add(MC, MM)), a);
                        barrier.await();
                    } catch (Exception e) {
                        Thread.currentThread().interrupt();
                        e.printStackTrace();
                    }
                });// matrix[1] = multiply(multiplyKahanBabushka(MM, add(MC, MM)), a)
                tread1.start();
                tread2.start();
                try {
                    barrier.await();
                    MF_KahanBabushka = add(matrix[0], matrix[1]);
                    timeMatrixKahanBabushka.set(System.currentTimeMillis() - startTime);
                    try {
                        outputSemaphore.acquire();
                        System.out.println("\nMatrix MF_Kahan-Babushka with threads:");
                        output(MF_KahanBabushka, amount);
                        writeToFile("MF_kahan_babushka.txt", MF_KahanBabushka);
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                    } finally {
                        outputSemaphore.release();
                    }
                } catch (Exception e) {
                    Thread.currentThread().interrupt();
                    e.printStackTrace();
                }
            } finally {
                latch.countDown();
            }
        });*/


        /*threadVectorSimpleModified.start();
        threadMatrixSimpleModified.start();
        threadVectorKahanModified.start();
        threadMatrixKahanModified.start();
        threadVectorKahanBabushkaModified.start();
        threadMatrixKahanBabushkaModified.start();*/

        /*try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }*/

        System.out.println("All time "+ (System.currentTimeMillis()-timeS));
        System.out.println("Time vector simple = "+timeVectorSimple);
        System.out.println("Time matrix simple = "+timeMatrixSimple);
        System.out.println("Time vector Kahan = "+timeVectorKahan);
        System.out.println("Time matrix Kahan = "+timeMatrixKahan);
        System.out.println("Time vector Kahan-Babushka = "+timeVectorKahanBabushka);
        System.out.println("Time matrix Kahan-Babushka = "+timeMatrixKahanBabushka);
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

    /**--------------------------Standard Methods----------------------------**/

    public static void output(double a) {
        System.out.println(a);
    }

    public static void output(double[] vector) {
        for (double v : vector) {
            System.out.println(v);
        }
    }

    public static void output(double[] vector, int amount) {
        if(vector.length<=2*amount) {
            output(vector);
            return;
        }

        for(int i = 0;i<amount; i++) {
            System.out.println(vector[i]);
        }
        System.out.println("∙∙∙");
        for(int i = vector.length-amount; i<vector.length; i++) {
            System.out.println(vector[i]);
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
