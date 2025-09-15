import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class GasDiffusion {
    private final Particle[] particles;
    private final Double v;
    private final Double radius;
    private final int N;
    private double currentTime;
    private final List<String> outputData;

    // Domain dimensions
    private final double leftWidth = 0.09;
    private final double leftHeight = 0.09;
    private final double rightWidth = 0.09;
    private final double rightTopY;
    private final double rightBottomY;

    public GasDiffusion(Double L, Double v, Double radius, int N) {
        this.v = v;
        this.radius = radius;
        this.N = N + 2;
        this.particles = new Particle[N + 2];
        this.outputData = new ArrayList<>();
        this.currentTime = 0;
        this.rightTopY = (leftHeight + L) / 2;
        this.rightBottomY = (leftHeight - L) / 2;
        initializeParticles();
    }

    private void initializeParticles() {
        Random random = new Random();

        for (int i = 0; i < N - 2; i++) {
            boolean validPosition = false;
            double x = 0, y = 0;

            // Place particles only in the left chamber initially
            while (!validPosition) {
                x = radius + random.nextDouble() * (leftWidth - 2 * radius);
                y = radius + random.nextDouble() * (leftHeight - 2 * radius);

                validPosition = true;
                for (int j = 0; j < i; j++) {
                    double dx = x - particles[j].getX();
                    double dy = y - particles[j].getY();
                    double distance = Math.sqrt(dx * dx + dy * dy);
                    if (distance < 2 * radius) {
                        validPosition = false;
                        break;
                    }
                }
            }

            // Random velocity with fixed magnitude
            double angle = 2 * Math.PI * random.nextDouble();
            double vx = v * Math.cos(angle);
            double vy = v * Math.sin(angle);

            particles[i] = new Particle(i, x, y, vx, vy, radius);
        }
        particles[N-1] = new Particle(N-1, leftWidth, rightBottomY, 0, 0, 0.0001);
        particles[N-2] = new Particle(N-2, leftWidth, rightTopY, 0, 0, 0.0001);
    }

    public void runSimulation(double totalTime) {
        outputData.add("t0");
        for (Particle p : particles) {
            outputData.add(p.toString());
        }

        PriorityQueue<Event> eventQueue = new PriorityQueue<>();

        // Initialize event queue with all possible events
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double collisionTime = getParticleCollisionTime(particles[i], particles[j]);
                if (collisionTime > 0) {
                    eventQueue.add(new Event(collisionTime, i, j, EventType.PARTICLE_COLLISION));
                }
            }

            List<WallCollision> wallCollisions = getWallCollisionTimes(particles[i]);
            for (WallCollision wc : wallCollisions) {
                if (wc.time > 0) {
                    eventQueue.add(new Event(wc.time, i, -1, wc.type));
                }
            }
        }

        double step = 0.1; // Output every 0.01 seconds
        double nextOutputTime = step;
        int outputCount = 1;

        while (currentTime < totalTime && !eventQueue.isEmpty()) {
            Event nextEvent = eventQueue.poll();

            // Skip invalid events
            if (nextEvent.time < currentTime) continue;

            // Move all particles to the event time
            double dt = nextEvent.time - currentTime;
            for (Particle p : particles) {
                p.move(dt);
            }
            currentTime = nextEvent.time;

            // Handle the event
            if (nextEvent.type == EventType.PARTICLE_COLLISION) {
                handleParticleCollision(particles[nextEvent.particle1], particles[nextEvent.particle2]);
            } else {
                handleWallCollision(particles[nextEvent.particle1], nextEvent.type);
            }

            // Update events for the particles involved
            updateEvents(eventQueue, nextEvent);

            // Output if needed
            if (currentTime >= nextOutputTime) {
                outputData.add("t" + outputCount);
                for (Particle p : particles) {
                    outputData.add(p.toString());
                }
                nextOutputTime += step;
                outputCount++;
            }
        }
    }

    private void updateEvents(PriorityQueue<Event> eventQueue, Event handledEvent) {
        // Remove all events involving the particles from the handled event
        eventQueue.removeIf(event -> (event.particle1 != -1 && (event.particle1 == handledEvent.particle1 || event.particle1 == handledEvent.particle2)) ||
                (event.particle2 != -1 && (event.particle2 == handledEvent.particle1 || event.particle2 == handledEvent.particle2)));

        // Add new events for the particles involved
        int p1 = handledEvent.particle1;
        int p2 = handledEvent.particle2;

        // For particle-particle events
        if (p2 != -1) {
            for (int i = 0; i < N; i++) {
                if (i != p1 && i != p2) {
                    double collisionTime = getParticleCollisionTime(particles[p1], particles[i]);
                    if (collisionTime > 0) {
                        eventQueue.add(new Event(currentTime + collisionTime, p1, i, EventType.PARTICLE_COLLISION));
                    }

                    collisionTime = getParticleCollisionTime(particles[p2], particles[i]);
                    if (collisionTime > 0) {
                        eventQueue.add(new Event(currentTime + collisionTime, p2, i, EventType.PARTICLE_COLLISION));
                    }
                }
            }
        }

        // For wall events
        List<WallCollision> wallCollisions1 = getWallCollisionTimes(particles[p1]);
        for (WallCollision wc : wallCollisions1) {
            if (wc.time > 0) {
                eventQueue.add(new Event(currentTime + wc.time, p1, -1, wc.type));
            }
        }

        if (p2 != -1) {
            List<WallCollision> wallCollisions2 = getWallCollisionTimes(particles[p2]);
            for (WallCollision wc : wallCollisions2) {
                if (wc.time > 0) {
                    eventQueue.add(new Event(currentTime + wc.time, p2, -1, wc.type));
                }
            }
        }
    }

    private double getParticleCollisionTime(Particle p1, Particle p2) {
        // Vector of relative position
        double dx = p2.getX() - p1.getX();
        double dy = p2.getY() - p1.getY();

        // Vector of relative velocity
        double dvx = p2.getVx() - p1.getVx();
        double dvy = p2.getVy() - p1.getVy();

        // Quadratic equation coefficients: a*t² + b*t + c = 0
        double a = dvx * dvx + dvy * dvy;
        double b = 2 * (dx * dvx + dy * dvy);
        double c = dx * dx + dy * dy - 4 * radius * radius;

        // Discriminant
        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0 || a == 0) {
            return -1; // No collision
        }

        double t1 = (-b - Math.sqrt(discriminant)) / (2 * a);
        double t2 = (-b + Math.sqrt(discriminant)) / (2 * a);

        // Return the smallest positive time
        if (t1 > 0 && t2 > 0) {
            return Math.min(t1, t2);
        } else if (t1 > 0) {
            return t1;
        } else if (t2 > 0) {
            return t2;
        } else {
            return -1; // Collision in the past
        }
    }

    private List<WallCollision> getWallCollisionTimes(Particle particle) {
        List<WallCollision> collisions = new ArrayList<>();
        final double EPS = 1e-12;

        double x = particle.getX();
        double y = particle.getY();
        double vx = particle.getVx();
        double vy = particle.getVy();
        double r = particle.getRadius();

        // 1) Leftmost wall (x = 0) - always solid
        if (vx < -EPS) {
            double t = (r - x) / vx; // vx < 0 -> t positive if heading left
            if (t > EPS) collisions.add(new WallCollision(t, EventType.LEFT_WALL));
        }

        // 2) Rightmost wall (x = leftWidth + rightWidth) - always solid
        double rightmostX = leftWidth + rightWidth;
        if (vx > EPS) {
            double t = (rightmostX - r - x) / vx;
            if (t > EPS) collisions.add(new WallCollision(t, EventType.RIGHT_WALL));
        }

        // 3) Middle vertical wall at x = leftWidth (the separating wall with an opening)
        //    Opening vertical interval is [rightBottomY, rightTopY].
        //    If the particle would intersect the solid part of the middle wall, we produce a middle-wall event.
        double middleX = leftWidth;
        // From LEFT -> test crossing middle from left side
        if (vx > EPS && x < middleX - r) {
            double t = (middleX - r - x) / vx;
            if (t > EPS) {
                double yAt = y + vy * t;
                // if yAt lies outside opening (consider particle radius), it's a collision
                if (yAt < rightBottomY + r || yAt > rightTopY - r) {
                    collisions.add(new WallCollision(t, EventType.MIDDLE_WALL));
                }
            }
        }
        // From RIGHT -> test crossing middle from right side
        if (vx < -EPS && x > middleX + r) {
            double t = (middleX + r - x) / vx; // vx negative -> t positive if heading left
            if (t > EPS) {
                double yAt = y + vy * t;
                if (yAt < rightBottomY + r || yAt > rightTopY - r) {
                    collisions.add(new WallCollision(t, EventType.MIDDLE_WALL));
                }
            }
        }

        // 4) Top/bottom walls are always solid
        // For left chamber
        if (x < leftWidth - r + EPS) {
            if (vy < -EPS) {
                double t = (r - y) / vy;
                if (t > EPS) collisions.add(new WallCollision(t, EventType.BOTTOM_WALL));
            } else if (vy > EPS) {
                double t = (leftHeight - r - y) / vy;
                if (t > EPS) collisions.add(new WallCollision(t, EventType.TOP_WALL));
            }
        }
        // For right chamber
        else if (x > leftWidth + r - EPS) {
            if (vy < -EPS) {
                double t = (rightBottomY + r - y) / vy;
                if (t > EPS) collisions.add(new WallCollision(t, EventType.BOTTOM_WALL));
            } else if (vy > EPS) {
                double t = (rightTopY - r - y) / vy;
                if (t > EPS) collisions.add(new WallCollision(t, EventType.TOP_WALL));
            }
        }
        // For particles near the middle wall, check both sets of top/bottom
        /* else {
            if (vy < -EPS) {
                double t1 = (r - y) / vy;
                if (t1 > EPS) collisions.add(new WallCollision(t1, EventType.BOTTOM_WALL));
                double t2 = (rightBottomY + r - y) / vy;
                if (t2 > EPS) collisions.add(new WallCollision(t2, EventType.BOTTOM_WALL));
            } else if (vy > EPS) {
                double t1 = (leftHeight - r - y) / vy;
                if (t1 > EPS) collisions.add(new WallCollision(t1, EventType.TOP_WALL));
                double t2 = (rightTopY - r - y) / vy;
                if (t2 > EPS) collisions.add(new WallCollision(t2, EventType.TOP_WALL));
            }
        } */

        return collisions;
    }

    private void handleParticleCollision(Particle p1, Particle p2) {
        // Normal vector from p1 to p2
        double nx = p2.getX() - p1.getX();
        double ny = p2.getY() - p1.getY();
        double distance = Math.sqrt(nx * nx + ny * ny);
        nx /= distance;
        ny /= distance;

        // Tangent vector
        double tx = -ny;
        double ty = nx;

        // Dot product of velocity with normal and tangent
        double v1n = p1.getVx() * nx + p1.getVy() * ny;
        double v1t = p1.getVx() * tx + p1.getVy() * ty;
        double v2n = p2.getVx() * nx + p2.getVy() * ny;
        double v2t = p2.getVx() * tx + p2.getVy() * ty;

        // Convert back to x,y coordinates
        if (p1.getId() < N - 2) {
            p1.setVx(v2n * nx + v1t * tx);
            p1.setVy(v2n * ny + v1t * ty);
        }
        if (p2.getId() < N - 2) {
            p2.setVx(v1n * nx + v2t * tx);
            p2.setVy(v1n * ny + v2t * ty);
        }
    }

    private void handleWallCollision(Particle particle, EventType wallType) {
        double r = particle.getRadius();

        switch (wallType) {
            case LEFT_WALL:
                particle.setX(r);
                particle.setVx(Math.abs(particle.getVx()));
                break;
            case RIGHT_WALL:
                particle.setX(leftWidth + rightWidth - r);
                particle.setVx(-Math.abs(particle.getVx()));
                break;
            case TOP_WALL:
                // determine which top: if particle.x < leftWidth => left top at leftHeight else right top at rightTopY
                if (particle.getX() < leftWidth) {
                    particle.setY(leftHeight - r);
                } else {
                    particle.setY(rightTopY - r);
                }
                particle.setVy(-Math.abs(particle.getVy()));
                break;
            case BOTTOM_WALL:
                if (particle.getX() < leftWidth) {
                    particle.setY(r);
                } else {
                    particle.setY(rightBottomY + r);
                }
                particle.setVy(Math.abs(particle.getVy()));
                break;
            case MIDDLE_WALL:
                particle.setX(leftWidth - r);
                particle.setVx(-Math.abs(particle.getVx()));
                break;
            default:
                break;
        }

    }


    public void outputToFile(String filename) {
        try (FileWriter writer = new FileWriter(filename)) {
            for (String line : outputData) {
                writer.write(line + "\n");
            }
        } catch (IOException e) {
            System.err.println("Error writing to file: " + e.getMessage());
        }
    }

    // Helper classes
    private enum EventType {
        PARTICLE_COLLISION, LEFT_WALL, RIGHT_WALL, TOP_WALL, BOTTOM_WALL, MIDDLE_WALL
    }

    private static class Event implements Comparable<Event> {
        double time;
        int particle1;
        int particle2;
        EventType type;

        Event(double time, int particle1, int particle2, EventType type) {
            this.time = time;
            this.particle1 = particle1;
            this.particle2 = particle2;
            this.type = type;
        }

        @Override
        public int compareTo(Event other) {
            return Double.compare(this.time, other.time);
        }
    }

    private static class WallCollision {
        double time;
        EventType type;

        WallCollision(double time, EventType type) {
            this.time = time;
            this.type = type;
        }
    }

}
