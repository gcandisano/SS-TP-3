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
    
    def create_animation(self, output_gif='particle_animation.gif', output_mp4=None, fps=10):
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
        print("Saving GIF animation...")
        anim.save(output_gif, writer='pillow', fps=fps, dpi=100)
        print(f"GIF animation saved as {output_gif}")
        
        # Save as MP4 video if requested
        if output_mp4:
            try:
                print("Saving MP4 video...")
                # Use ffmpeg writer for MP4
                anim.save(output_mp4, writer='ffmpeg', fps=fps, 
                         bitrate=1800, dpi=100,
                         extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])
                print(f"MP4 video saved as {output_mp4}")
            except Exception as e:
                print(f"Error saving MP4: {e}")
                print("Make sure ffmpeg is installed:")
                print("Windows: Download from https://ffmpeg.org/")
                print("Mac: brew install ffmpeg")
                print("Linux: sudo apt-get install ffmpeg")
        
        plt.close()
        
        return anim

def create_animation_from_file(filename, output_gif=None, output_mp4=None, right_height=0.09):
    if output_gif is None:
        output_gif = f"{os.path.splitext(filename)[0]}_animation.gif"
    
    if output_mp4 is None:
        output_mp4 = f"{os.path.splitext(filename)[0]}_animation.mp4"
    
    try:
        animator = ParticleAnimation(filename, right_height=right_height)
        animator.create_animation(output_gif, output_mp4)
        print(f"Successfully created animations: {output_gif} and {output_mp4}")
        return True
    except Exception as e:
        print(f"Error creating animation: {e}")
        return False

# Batch processing function for multiple files
def create_animations_for_all_simulations():
    simulation_files = [
        ("simulation_L_0.03.txt", 0.03),
        ("simulation_L_0.05.txt", 0.05), 
        ("simulation_L_0.07.txt", 0.07),
        ("simulation_L_0.09.txt", 0.09)
    ]
    
    for file, height in simulation_files:
        if os.path.exists(file):
            print(f"Processing {file} with height {height}...")
            create_animation_from_file(file, right_height=height)
        else:
            print(f"File {file} not found, skipping...")

# Function to check if ffmpeg is available
def check_ffmpeg():
    try:
        import subprocess
        result = subprocess.run(['ffmpeg', '-version'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ ffmpeg is available for MP4 creation")
            return True
        else:
            print("✗ ffmpeg not found. MP4 videos will not be created.")
            print("Install ffmpeg:")
            print("Windows: Download from https://ffmpeg.org/")
            print("Mac: brew install ffmpeg")
            print("Linux: sudo apt-get install ffmpeg")
            return False
    except:
        print("✗ ffmpeg not found. MP4 videos will not be created.")
        return False

if __name__ == "__main__":
    # Check if ffmpeg is available
    has_ffmpeg = check_ffmpeg()
    
    # Example usage for a specific file
    filename = "simulation_L_0.03.txt"
    right_height = 0.03
    
    if os.path.exists(filename):
        print(f"Creating animation for {filename}...")
        if has_ffmpeg:
            create_animation_from_file(filename, right_height=right_height)
        else:
            # Create only GIF if ffmpeg is not available
            create_animation_from_file(filename, output_mp4=None, right_height=right_height)
    else:
        print(f"File {filename} not found!")
    
    # Or process all simulation files
    # create_animations_for_all_simulations()