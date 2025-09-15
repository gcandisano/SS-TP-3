package ar.edu.itba.ss;

import lombok.Getter;
import lombok.Setter;

import java.util.Locale;

@Getter
@Setter
public class Particle {
    private int id;
    private Double x;
    private Double y;
    private Double vy;
    private Double vx;
    private Double radius;
    private int collided = 0;

    public Particle(int id, double x, double y, double vx, double vy, double radius) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.radius = radius;
    }

    public void move(double dt) {
        x += vx * dt;
        y += vy * dt;
    }

    @Override
    public String toString() {
        return String.format(Locale.ENGLISH, "%.6f %.6f %.6f %.6f %d", x, y, vx, vy, collided);
    }
}