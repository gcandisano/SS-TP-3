package ar.edu.itba.ss;

public class Event implements Comparable<Event> {
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

    @Override
    public String toString() {
        return "Event{" +
                "time=" + time +
                ", particle1=" + particle1 +
                ", particle2=" + particle2 +
                ", type=" + type +
                '}';
    }
}