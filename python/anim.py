import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle, Circle
import re
import os
import subprocess

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
        # Removed the title as requested
        self.ax.set_xlabel('X position (m)', fontsize=16)  # Increased font size
        self.ax.set_ylabel('Y position (m)', fontsize=16)  # Increased font size
        
        # Draw the chamber walls
        self.draw_chamber()
        
        # Create particle objects
        self.particles = []
        for i in range(self.num_particles):
            circle = Circle((0, 0), radius=0.0015, color='blue', alpha=0.7)
            self.particles.append(self.ax.add_patch(circle))
        
        # Add time text with the requested format
        self.time_text = self.ax.text(0.02, 0.095, '', fontsize=14)
        
    def draw_chamber(self):
        # Draw a single continuous chamber in L-shape (inverted L)
        # Left chamber: full height
        # Right chamber: only at the top, creating an inverted L shape
        
        from matplotlib.patches import Polygon
        
        # Define the points of the L-shaped chamber
        # The right chamber is at the TOP, not in the middle
        points = [
            (0, 0),  # Bottom left
            (self.left_width, 0),  # Bottom right of left chamber
            (self.left_width, self.right_y0),  # Bottom of opening
            (self.left_width + self.right_width, self.right_y0),  # Bottom right of right chamber
            (self.left_width + self.right_width, self.right_y1),  # Top right of right chamber
            (self.left_width, self.right_y1),  # Top of opening
            (self.left_width, self.left_height),  # Top left of left chamber
            (0, self.left_height),  # Top left
            (0, 0)  # Close the polygon
        ]
        
        # Create the L-shaped chamber polygon
        chamber = Polygon(points, fill=False, color='black', linewidth=2)
        self.ax.add_patch(chamber)

    
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

                if speed == 0:
                    self.particles[i].set_radius(0.0)
        
        # Update time text with the requested format
        time_seconds = float(time_str[1:]) * 0.1  # Convert from frame number to seconds
        self.time_text.set_text(f'Time: {time_seconds:.1f} seconds')
        
        return self.particles + [self.time_text]
    
    def makeMP4(self, output_mp4=None, fps=10, sample_interval=5):
        """Create MP4 video only - much faster and smaller than GIF"""
        if output_mp4 is None:
            output_mp4 = f"{os.path.splitext(self.filename)[0]}_animation.mp4"
        
        print(f"Creating MP4 animation with {len(self.times)} frames...")
        
        # Sample frames to reduce the total number (much faster)
        if sample_interval > 1:
            sampled_frames = list(range(0, len(self.times), sample_interval))
            total_frames = len(sampled_frames)
            print(f"Sampling every {sample_interval} frames: {total_frames} total frames")
        else:
            sampled_frames = list(range(len(self.times)))
            total_frames = len(self.times)
        
        # Create animation with sampled frames
        anim = animation.FuncAnimation(
            self.fig, 
            self.update, 
            frames=total_frames,
            interval=1000/fps,
            blit=True
        )
        
        try:
            print("Saving MP4 video...")
            anim.save(output_mp4, writer='ffmpeg', fps=fps, 
                     bitrate=1800, dpi=100,
                     extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])
            print(f"✓ MP4 video saved as {output_mp4}")
            print(f"  Frames: {total_frames}, FPS: {fps}, Estimated size: 1-5MB")
        except Exception as e:
            print(f"✗ Error saving MP4: {e}")
            if "ffmpeg" in str(e).lower():
                print("Please install ffmpeg:")
                print("Windows: Download from https://ffmpeg.org/ and add to PATH")
                print("Mac: brew install ffmpeg")
                print("Linux: sudo apt-get install ffmpeg")
            return False
        
        plt.close()
        return True

def create_MP4_from_file(filename, right_height=0.09, fps=10, sample_interval=1):
    """Create only MP4 video from simulation file"""
    if not os.path.exists(filename):
        print(f"✗ File {filename} not found!")
        return False
    
    try:
        animator = ParticleAnimation(filename, right_height=right_height)
        output_mp4 = f"{os.path.splitext(filename)[0]}_animation.mp4"
        success = animator.makeMP4(output_mp4, fps=fps, sample_interval=sample_interval)
        return success
    except Exception as e:
        print(f"✗ Error creating animation: {e}")
        return False

def create_MP4s_for_all_simulations(fps=10):
    """Create MP4 videos for all simulation files"""
    simulation_files = [
        ("simulation_L_0.03.txt", 0.03),
        ("simulation_L_0.05.txt", 0.05),
        ("simulation_L_0.07.txt", 0.07),
        ("simulation_L_0.09.txt", 0.09)
    ]
    
    for file, height in simulation_files:
        if os.path.exists(file):
            print(f"\nProcessing {file} with height {height}...")
            create_MP4_from_file(file, right_height=height, fps=fps)
        else:
            print(f"✗ File {file} not found, skipping...")

# Function to check if ffmpeg is available
def check_ffmpeg():
    try:
        result = subprocess.run(['ffmpeg', '-version'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ ffmpeg is available for MP4 creation")
            return True
        else:
            print("✗ ffmpeg not found. MP4 videos cannot be created.")
            print("Please install ffmpeg:")
            print("Windows: Download from https://ffmpeg.org/ and add to PATH")
            print("Mac: brew install ffmpeg")
            print("Linux: sudo apt-get install ffmpeg")
            return False
    except:
        print("✗ ffmpeg not found. MP4 videos cannot be created.")
        return False

if __name__ == "__main__":
    # Check if ffmpeg is available
    if not check_ffmpeg():
        exit(1)
    
    # Process all simulation files (recommended)
    create_MP4s_for_all_simulations(fps=30)
    
    print("\n✓ All MP4 videos created successfully!")
    print("MP4 files are much faster to create and smaller than GIF files.")