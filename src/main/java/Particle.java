
import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class Particle {
    Double x;
    Double y;
    Double vy;
    Double vx;
    Double radious;

    public Particle(double x, double y, double vx, double vy, double radious) {
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.radious = radious;
    }

    public Double distanceTo(Particle particle) {
        return Math.sqrt(Math.pow(this.x - particle.getX(), 2) + Math.pow(this.y - particle.getY(), 2));
    }
}
