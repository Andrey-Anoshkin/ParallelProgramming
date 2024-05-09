import java.util.ArrayList;

class myThread extends Thread {
    myThread(double a, double b, double dx, int threadId, int threadCount) {
        this.a = a;
        this.b = b;
        this.dx = dx;
        this.threadId = threadId;
        this.threadCount = threadCount;
    }

    public double getSum() {
        return this.sum;
    }

    private double f(double x) {
        return 1.0 / (1 + x * x);
    }

    private double integral(double a, double b, double dx, int threadId, int threadCount) {
        double sum = 0;
        int n = (int)((b - a) / dx);

        for (int i = threadId; i < n; i += threadCount) {
            double x = a + i * dx + dx / 2;
            sum += this.f(x) * dx;
        }

        return sum;
    }

    @Override
    public void run() {
        this.sum = this.integral(this.a, this.b, this.dx, this.threadId, this.threadCount);
    }

    private int threadId = 0;
    private int threadCount = 0;
    private double sum = 0;
    private double a = 0;
    private double b = 0;
    private double dx = 0;
}

public class Main {
    public static double integral(int threadCount, ArrayList<myThread> threads) {
        double sum = 0;
        // Запуск потоков
        for (int i = 0; i < threadCount; ++i)
            threads.get(i).start();

        boolean isAlive;

        // Ожидание завершения выполнения всех потоков
        do {
            isAlive = false;
            for (int i = 0; i < threadCount; ++i)
                isAlive = isAlive || threads.get(i).isAlive();
        }
        while (isAlive);

        // Сбор информации из всех потоков
        for (int i = 0; i < threadCount; ++i)
            sum += threads.get(i).getSum();

        return sum;
    }

    public static void main(String[] args) {
        double sum = 0;
        int threadCount = 1;
        long startTime, finishTime;

        // Объявление потоков
        ArrayList<myThread> threads = new ArrayList<myThread>();

        for (int i = 0; i < threadCount; ++i)
            threads.add(new myThread(0, 1e+6, 0.001, i, threadCount));

        // Начало выполнения
        startTime = System.nanoTime();

        // Вычисление
        sum = integral(threadCount, threads);

        // Конец выполнения
        finishTime = System.nanoTime();

        // Завершение потоков
        for (int i = 0; i < threadCount; ++i)
            threads.get(i).interrupt();

        // Вывод
        System.out.println("Integral value = " + sum);
        System.out.println("Number of threads = " + threadCount);
        System.out.println("Execution time = " + (finishTime - startTime) / 1e+9);
    }
}