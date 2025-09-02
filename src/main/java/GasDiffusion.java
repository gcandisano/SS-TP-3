public class GasDiffusion {
    Particle[] particles;
    private final Double L;
    private final Double v;
    private final Double radious;
    private final int N;
    private final Double deltaT;

    public GasDiffusion(Double L, Double v, Double radious, int N, Double deltaT) {
        this.L = L;
        this.v = v;
        this.radious = radious;
        this.N = N;
        this.deltaT = deltaT;
        this.particles = new Particle[N];
        initializeParticles();
    }

    private void initializeParticles() {
        for (int i = 0; i < N; i++) {
            double x = Math.random() * L;
            double y = Math.random() * L;
            for (int j = 0; j < i; j++) {
                if (particles[j].getX() == x && particles[j].getY() == y) {
                    i--;
                    break;
                }
            }
            particles[i] = new Particle(x, y, 0.0, 0.0, radious);
        }
    }

}
