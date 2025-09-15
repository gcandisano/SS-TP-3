import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle, Circle
import re
import os

class ParticleAnimation:
    def __init__(self, filename, right_height=0.09):
        self.filename = filename
        self.data = self.parse_file()
        self.times = list(self.data.keys())
        self.num_particles = len(self.data[self.times[0]])
        
        # Box geometry
        self.left_width = 0.09
        self.left_height = 0.09
        self.right_width = 0.09
        self.right_height = right_height
        self.right_y0 = (self.left_height - self.right_height) / 2  # bottom of right box
        self.right_y1 = self.right_y0 + self.right_height          # top of right box
        
        # Set up the figure and axis
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        self.ax.set_xlim(-0.01, 0.19)  # Slightly beyond walls
        self.ax.set_ylim(-0.01, 0.10)  # Slightly beyond walls
        self.ax.set_aspect('equal')
        self.ax.set_title('Gas Diffusion Simulation')
        self.ax.set_xlabel('X position (m)')
        self.ax.set_ylabel('Y position (m)')
        
        # Draw the chamber walls
        self.draw_chamber()
        
        # Create particle objects
        self.particles = []
        for i in range(self.num_particles):
            circle = Circle((0, 0), radius=0.0015, color='blue', alpha=0.7)
            self.particles.append(self.ax.add_patch(circle))
        
        # Add time text and event counter
        self.time_text = self.ax.text(0.02, 0.095, '', fontsize=12)
        self.event_text = self.ax.text(0.02, 0.085, '', fontsize=10, color='red')
        
        # Add legend for particle colors
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='blue', alpha=0.7, label='Partículas normales'),
            Patch(facecolor='red', alpha=1.0, label='Partículas colisionadas')
        ]
        self.ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
        
    def draw_chamber(self):
        # Left chamber (0.09 x 0.09)
        left_chamber = Rectangle((0, 0), self.left_width, self.left_height,
                                 fill=False, color='black', linewidth=2)
        self.ax.add_patch(left_chamber)
        
        # Right chamber, vertically centered
        right_chamber = Rectangle((self.left_width, self.right_y0),
                                  self.right_width, self.right_height,
                                  fill=False, color='black', linewidth=2)
        self.ax.add_patch(right_chamber)
        
        # Middle opening (the passage) – same height as right chamber
        middle_opening = Rectangle((self.left_width, self.right_y0),
                                   0.0001, self.right_height,
                                   fill=True, color='white', linewidth=0)
        self.ax.add_patch(middle_opening)

    
    def parse_file(self):
        data = {}
        current_time = None
        current_particles = []
        
        with open(self.filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('t'):
                    if current_time is not None:
                        data[current_time] = current_particles
                    current_time = line
                    current_particles = []
                else:
                    # Parse particle data: x y vx vy collided
                    parts = line.split()
                    if len(parts) == 5:
                        try:
                            x, y, vx, vy = map(float, parts[:4])
                            collided = parts[4].lower() == '1'
                            current_particles.append((x, y, vx, vy, collided))
                        except ValueError:
                            continue
        
        if current_time is not None and current_particles:
            data[current_time] = current_particles
        
        return data
    
    def update(self, frame):
        time_idx = frame % len(self.times)
        time_str = self.times[time_idx]
        particles_data = self.data[time_str]
        
        # Update particle positions
        for i, particle_data in enumerate(particles_data):
            if i < len(self.particles):
                x, y, vx, vy, collided = particle_data
                self.particles[i].center = (x, y)
                
                # Color particles based on collision status
                if collided:
                    # Red color for particles that just collided
                    self.particles[i].set_color('red')
                    self.particles[i].set_alpha(1.0)  # Full opacity for collided particles
                else:
                    # Blue color for normal particles
                    self.particles[i].set_color('blue')
                    self.particles[i].set_alpha(0.7)  # Slightly transparent for normal particles
        
        # Update time text with actual event time
        if time_str.startswith('t '):
            # Extract numeric time value for display
            time_value = time_str.split(' ')[1]
            self.time_text.set_text(f'Event Time: {time_value}s')
        else:
            self.time_text.set_text(f'Time: {time_str}')
        
        # Update event counter
        self.event_text.set_text(f'Event #{time_idx + 1} of {len(self.times)}')
        
        return self.particles + [self.time_text, self.event_text]
    
    def create_animation(self, output_gif='particle_animation.gif', fps=10, pause_duration=800):
        print(f"Creating fast event-based animation with {len(self.times)} events...")
        
        # Faster animation with shorter pauses between events
        anim = animation.FuncAnimation(
            self.fig, 
            self.update, 
            frames=len(self.times),
            interval=pause_duration,  # ms per frame - faster transitions
            blit=True,
            repeat=True
        )
        
        # Save as GIF with higher fps for faster viewing
        print("Saving fast event-based animation...")
        anim.save(output_gif, writer='pillow', fps=fps, dpi=100)
        print(f"Fast event-based animation saved as {output_gif}")
        
        plt.close()
        
        return anim

def create_fast_animation(filename, output_gif=None):
    """Create a very fast animation for quick overview"""
    return create_animation_from_file(filename, output_gif, fps=20, pause_duration=50)

def create_medium_animation(filename, output_gif=None):
    """Create a medium speed animation for balanced viewing"""
    return create_animation_from_file(filename, output_gif, fps=15, pause_duration=100)

def create_slow_animation(filename, output_gif=None):
    """Create a slow animation for detailed event analysis"""
    return create_animation_from_file(filename, output_gif, fps=8, pause_duration=300)

def create_animation_from_file(filename, output_gif=None, fps=15, pause_duration=100):
    if output_gif is None:
        output_gif = f"{os.path.splitext(filename)[0]}_animation.gif"
    
    try:
        # Extract L value from filename for right_height parameter
        import re
        match = re.search(r'L_(\d+\.\d+)', filename)
        if match:
            right_height = float(match.group(1))
        else:
            right_height = 0.03  # default
        
        animator = ParticleAnimation(filename, right_height=right_height)
        animator.create_animation(output_gif, fps=fps, pause_duration=pause_duration)
        print(f"Successfully created event-based animation: {output_gif}")
        return True
    except Exception as e:
        print(f"Error creating animation: {e}")
        return False

# Batch processing function for multiple files
def create_animations_for_all_simulations():
    simulation_files = [
        "simulation_L_0.03.txt",
        "simulation_L_0.05.txt", 
        "simulation_L_0.07.txt",
        "simulation_L_0.09.txt"
    ]
    
    for file in simulation_files:
        if os.path.exists(file):
            print(f"Processing {file}...")
            create_animation_from_file(file)
        else:
            print(f"File {file} not found, skipping...")

if __name__ == "__main__":
    # Example usage:
    # Create fast event-based animation for a specific file
    animator = ParticleAnimation("simulation_L_0.03.txt", right_height=0.03)
    animator.create_animation("gas_diffusion_L_0.03.gif", fps=10, pause_duration=800)
    
    # Or process all simulation files
    # create_animations_for_all_simulations()