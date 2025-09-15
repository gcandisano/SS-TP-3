import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class Particle {
    private int id;
    private Double x;
    private Double y;
    private Double vy;
    private Double vx;
    private Double radius;

    public Particle(int id, double x, double y, double vx, double vy, double radius) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.radius = radius;
    }

    public Double distanceTo(Particle particle) {
        return Math.sqrt(Math.pow(this.x - particle.getX(), 2) + Math.pow(this.y - particle.getY(), 2));
    }

    public void move(double dt) {
        x += vx * dt;
        y += vy * dt;
    }

    @Override
    public String toString() {
        return String.format("%.6f %.6f %.6f %.6f", x, y, vx, vy);
    }
}