package ar.edu.itba.ss;

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
        particles[N - 1] = new Particle(N - 1, leftWidth, rightBottomY, 0, 0, 0.0001);
        particles[N - 2] = new Particle(N - 2, leftWidth, rightTopY, 0, 0, 0.0001);
    }

    public void runSimulation(double totalTime) {
        outputData.add("t0");
        for (Particle p : particles) {
            outputData.add(p.toString());
        }

        PriorityQueue<Event> eventQueue = new PriorityQueue<>();

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

            if (nextEvent.type == EventType.PARTICLE_COLLISION) {
                handleParticleCollision(particles[nextEvent.particle1], particles[nextEvent.particle2]);
            } else {
                handleWallCollision(particles[nextEvent.particle1], nextEvent.type);
            }

            // Update events for the particles involved
            updateEvents(eventQueue, nextEvent);

            // Validate no overlapping particles (safety check)
            validateNoOverlaps();

            String timeStr = String.format(Locale.ENGLISH, "t %.6f", currentTime);
            outputData.add(timeStr);
            for (Particle p : particles) {
                if (p.getId() == nextEvent.particle1 || p.getId() == nextEvent.particle2) {
                    p.setCollided(1);
                }
                outputData.add(p.toString());
                p.setCollided(0);
            }

        }
    }

    private void updateEvents(PriorityQueue<Event> eventQueue, Event handledEvent) {
        // Remove all events involving the particles from the handled event
        eventQueue.removeIf(event -> (event.particle1 != -1 && (event.particle1 == handledEvent.particle1 || event.particle1 == handledEvent.particle2)) ||
                (event.particle2 != -1 && (event.particle2 == handledEvent.particle1 || event.particle2 == handledEvent.particle2)));

        // Add new events for ALL particles (not just the ones involved in the current event)
        // This ensures no collisions are missed
        for (int i = 0; i < N - 2; i++) { // Exclude virtual wall particles
            // Particle-particle collisions
            for (int j = i + 1; j < N - 2; j++) {
                double collisionTime = getParticleCollisionTime(particles[i], particles[j]);
                if (collisionTime > 0) {
                    eventQueue.add(new Event(currentTime + collisionTime, i, j, EventType.PARTICLE_COLLISION));
                }
            }
            
            // Wall collisions
            List<WallCollision> wallCollisions = getWallCollisionTimes(particles[i]);
            for (WallCollision wc : wallCollisions) {
                if (wc.time > 0) {
                    eventQueue.add(new Event(currentTime + wc.time, i, -1, wc.type));
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

        // Return the smallest positive time with safety margin
        double collisionTime = -1;
        if (t1 > 1e-10 && t2 > 1e-10) {
            collisionTime = Math.min(t1, t2);
        } else if (t1 > 1e-10) {
            collisionTime = t1;
        } else if (t2 > 1e-10) {
            collisionTime = t2;
        }
        
        // Additional safety check: ensure particles are not already overlapping
        double currentDistance = Math.sqrt(dx * dx + dy * dy);
        if (currentDistance < 2 * radius - 1e-10) {
            return 1e-10; // Force immediate collision if already overlapping
        }
        
        return collisionTime;
    }

    private List<WallCollision> getWallCollisionTimes(Particle particle) {
        List<WallCollision> collisions = new ArrayList<>();
        final double EPS = 1e-8;

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
        double dx = p2.getX() - p1.getX();
        double dy = p2.getY() - p1.getY();
        double distance = Math.sqrt(dx * dx + dy * dy);

        // Ensure particles are properly separated after collision
        if (distance < 2 * radius) {
            double overlap = 2 * radius - distance;
            double separationX = (dx / distance) * overlap / 2;
            double separationY = (dy / distance) * overlap / 2;
            
            // Move particles apart
            p1.setX(p1.getX() - separationX);
            p1.setY(p1.getY() - separationY);
            p2.setX(p2.getX() + separationX);
            p2.setY(p2.getY() + separationY);
            
            // Recalculate relative position after separation
            dx = p2.getX() - p1.getX();
            dy = p2.getY() - p1.getY();
        }

        double dvx = p2.getVx() - p1.getVx();
        double dvy = p2.getVy() - p1.getVy();

        double J = 2 * (dvx * dx + dvy * dy) / (4 * radius);

        // Convert back to x,y coordinates
        if (p1.getId() < N - 2) {
            p1.setVx(J * dx / (2 * radius) + p1.getVx());
            p1.setVy(J * dy / (2 * radius) + p1.getVy());
        }
        if (p2.getId() < N - 2) {
            p2.setVx(-J * dx / (2 * radius) + p2.getVx());
            p2.setVy(-J * dy / (2 * radius) + p2.getVy());
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

    private void validateNoOverlaps() {
        // Safety check to ensure no particles are overlapping
        for (int i = 0; i < N - 2; i++) {
            for (int j = i + 1; j < N - 2; j++) {
                double dx = particles[j].getX() - particles[i].getX();
                double dy = particles[j].getY() - particles[i].getY();
                double distance = Math.sqrt(dx * dx + dy * dy);
                
                if (distance < 2 * radius - 1e-10) {
                    // Force separation if particles are overlapping
                    double overlap = 2 * radius - distance;
                    double separationX = (dx / distance) * overlap / 2;
                    double separationY = (dy / distance) * overlap / 2;
                    
                    particles[i].setX(particles[i].getX() - separationX);
                    particles[i].setY(particles[i].getY() - separationY);
                    particles[j].setX(particles[j].getX() + separationX);
                    particles[j].setY(particles[j].getY() + separationY);
                }
            }
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
