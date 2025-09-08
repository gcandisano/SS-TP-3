public class GasDiffusion {
    Particle[] particles;
    private final Double L;
    private final Double v;
    private final Double radious;
    private final int N;
    private final Double deltaT;
    private Double fixedLength;

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
                if (particles[j].distanceTo(new Particle(x, y, 0.0, 0.0, radious)) < 2 * radious) {
                    i--;
                    break;
                }
            }
            Double vx = Math.random() * v;
            Double vy = Math.sqrt(v * v - vx * vx);
            particles[i] = new Particle(x, y, vx, vy, radious);
        }
    }

    /**
     * Calcula el tiempo al próximo evento de colisión para una partícula dada
     * @param particle La partícula para la cual calcular el próximo evento
     * @return El tiempo al próximo evento de colisión
     */
    public double getNextCollisionTime(Particle particle) {
        double minTime = Double.MAX_VALUE;
        
        // 1. Calcular tiempo de colisión con otras partículas
        for (Particle other : particles) {
            if (other != particle) {
                double collisionTime = getParticleCollisionTime(particle, other);
                if (collisionTime > 0 && collisionTime < minTime) {
                    minTime = collisionTime;
                }
            }
        }
        
        // 2. Calcular tiempo de colisión con las paredes
        double wallCollisionTime = getWallCollisionTime(particle);
        if (wallCollisionTime > 0 && wallCollisionTime < minTime) {
            minTime = wallCollisionTime;
        }
        
        return minTime == Double.MAX_VALUE ? -1 : minTime;
    }
    
    /**
     * Calcula el tiempo de colisión entre dos partículas
     */
    private double getParticleCollisionTime(Particle p1, Particle p2) {
        // Vector de posición relativa
        double dx = p2.getX() - p1.getX();
        double dy = p2.getY() - p1.getY();
        
        // Vector de velocidad relativa
        double dvx = p2.getVx() - p1.getVx();
        double dvy = p2.getVy() - p1.getVy();
        
        // Coeficientes de la ecuación cuadrática: a*t² + b*t + c = 0
        double a = dvx * dvx + dvy * dvy;
        double b = 2 * (dx * dvx + dy * dvy);
        double c = dx * dx + dy * dy - 4 * radious * radious;
        
        // Discriminante
        double discriminant = b * b - 4 * a * c;
        
        if (discriminant < 0) {
            return -1; // No hay colisión
        }
        
        double t1 = (-b - Math.sqrt(discriminant)) / (2 * a);
        double t2 = (-b + Math.sqrt(discriminant)) / (2 * a);
        
        // Retornar el tiempo positivo más pequeño
        if (t1 > 0 && t2 > 0) {
            return Math.min(t1, t2);
        } else if (t1 > 0) {
            return t1;
        } else if (t2 > 0) {
            return t2;
        } else {
            return -1; // Colisión en el pasado
        }
    }
    
    /**
     * Calcula el tiempo de colisión con las paredes
     */
    private double getWallCollisionTime(Particle particle) {
        double minTime = Double.MAX_VALUE;
        
        // Colisión con pared izquierda (x = 0)
        if (particle.getVx() < 0) {
            double t = -particle.getX() / particle.getVx();
            if (t > 0 && t < minTime) {
                minTime = t;
            }
        }
        
        // Colisión con pared derecha (x = L)
        if (particle.getVx() > 0) {
            double t = (L - particle.getX()) / particle.getVx();
            if (t > 0 && t < minTime) {
                minTime = t;
            }
        }
        
        // Colisión con pared inferior (y = 0)
        if (particle.getVy() < 0) {
            double t = -particle.getY() / particle.getVy();
            if (t > 0 && t < minTime) {
                minTime = t;
            }
        }
        
        // Colisión con pared superior (y = L)
        if (particle.getVy() > 0) {
            double t = (L - particle.getY()) / particle.getVy();
            if (t > 0 && t < minTime) {
                minTime = t;
            }
        }
        
        return minTime == Double.MAX_VALUE ? -1 : minTime;
    }
}
