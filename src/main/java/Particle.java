
import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class Particle {
    Double x;
    Double y;

    double vy;
    Double vx;
    Double radious;


    public Particle(double x, double y, double vx, double vy, double radious) {
        this.x = x;
        this.y = y;

    }

}
