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
        
        # Add time text
        self.time_text = self.ax.text(0.02, 0.095, '', fontsize=12)
        
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
                    # Parse particle data: x y vx vy
                    parts = line.split()
                    if len(parts) == 4:
                        try:
                            x, y, vx, vy = map(float, parts)
                            current_particles.append((x, y, vx, vy))
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
        for i, (x, y, vx, vy) in enumerate(particles_data):
            if i < len(self.particles):
                self.particles[i].center = (x, y)
                
                # Optional: color particles based on velocity magnitude
                speed = np.sqrt(vx**2 + vy**2)
                norm_speed = min(speed / 0.02, 1.0)  # Normalize to 0-1
                color = plt.cm.plasma(norm_speed)
                self.particles[i].set_color(color)
        
        # Update time text
        self.time_text.set_text(f'Time: {time_str}')
        
        return self.particles + [self.time_text]
    
    def create_animation(self, output_gif='particle_animation.gif', fps=10):
        print(f"Creating animation with {len(self.times)} frames...")
        
        # Create animation
        anim = animation.FuncAnimation(
            self.fig, 
            self.update, 
            frames=len(self.times),
            interval=1000/fps,  # ms per frame
            blit=True
        )
        
        # Save as GIF
        print("Saving animation...")
        anim.save(output_gif, writer='pillow', fps=fps, dpi=100)
        print(f"Animation saved as {output_gif}")
        
        plt.close()
        
        return anim

def create_animation_from_file(filename, output_gif=None):
    if output_gif is None:
        output_gif = f"{os.path.splitext(filename)[0]}_animation.gif"
    
    try:
        animator = ParticleAnimation(filename)
        animator.create_animation(output_gif)
        print(f"Successfully created animation: {output_gif}")
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
    # Create animation for a specific file
    animator = ParticleAnimation("simulation_L_0.03.txt", right_height=0.03)
    animator.create_animation("gas_diffusion_L0.03.gif")
    
    # Or process all simulation files
    # create_animations_for_all_simulations()