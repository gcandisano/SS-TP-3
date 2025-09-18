package ar.edu.itba.ss;

public class Main {
    public static void main(String[] args) {
        Double radius = 0.0015;
        Double v = 0.01;
        double[] L_values = {0.03, 0.05, 0.07, 0.09};
        int N = 200; // Number of particles
        double totalTime = 60.0; // Total simulation time in seconds

        for (double L : L_values) {
            System.out.println("Simulating for L = " + L);
            GasDiffusion simulation = new GasDiffusion(L, v, radius, N);
            simulation.runSimulation(totalTime);
            simulation.outputToFile("simulation_L_" + L + ".txt");
            simulation.outputWallCollisions("wall_collisions_L_" + L + ".txt");
        }
    }
}