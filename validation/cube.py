import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


def create_small_cube(x, y, z, size):
  """Create vertices and faces for a small cube at position (x,y,z) with given size"""
  vertices = [
    [x, y, z],
    [x + size, y, z],
    [x + size, y + size, z],
    [x, y + size, z],
    [x, y, z + size],
    [x + size, y, z + size],
    [x + size, y + size, z + size],
    [x, y + size, z + size]
  ]

  faces = [
    [vertices[0], vertices[1], vertices[2], vertices[3]],  # bottom
    [vertices[4], vertices[5], vertices[6], vertices[7]],  # top
    [vertices[0], vertices[1], vertices[5], vertices[4]],  # front
    [vertices[2], vertices[3], vertices[7], vertices[6]],  # back
    [vertices[0], vertices[3], vertices[7], vertices[4]],  # left
    [vertices[1], vertices[2], vertices[6], vertices[5]]  # right
  ]

  return vertices, faces


if __name__ == "__main__":
  L = 6
  small_size = 1  # Size of each small cube

  fig = plt.figure(figsize=(12, 10))
  ax = fig.add_subplot(111, projection="3d")

  # Create outer cube vertices
  outer_vertices = [
    [0, 0, 0],
    [L, 0, 0],
    [L, L, 0],
    [0, L, 0],
    [0, 0, L],
    [L, 0, L],
    [L, L, L],
    [0, L, L]
  ]

  outer_faces = [
    [outer_vertices[0], outer_vertices[1], outer_vertices[2], outer_vertices[3]],
    [outer_vertices[4], outer_vertices[5], outer_vertices[6], outer_vertices[7]],
    [outer_vertices[0], outer_vertices[1], outer_vertices[5], outer_vertices[4]],
    [outer_vertices[2], outer_vertices[3], outer_vertices[7], outer_vertices[6]],
    [outer_vertices[0], outer_vertices[3], outer_vertices[7], outer_vertices[4]],
    [outer_vertices[1], outer_vertices[2], outer_vertices[6], outer_vertices[5]]
  ]

  # Draw the outer cube with very transparent fill
  outer_cube = Poly3DCollection(
    outer_faces,
    alpha=0.02,  # Even more transparent
    edgecolor="black",
    linewidth=3  # Thicker outer edges
  )
  ax.add_collection3d(outer_cube)

  # Draw surface cubes with more defined edges
  for x in range(0, L, small_size):
    for y in range(0, L, small_size):
      for z in range(0, L, small_size):
        # Check if the cube is on the surface
        is_on_surface = (
            x == 0 or x == L - small_size or
            y == 0 or y == L - small_size or
            z == 0 or z == L - small_size
        )

        if is_on_surface:
          _, faces = create_small_cube(x, y, z, small_size)
          # Make cubes more opaque with dark edges
          small_cube = Poly3DCollection(
            faces,
            alpha=0.6,  # More opaque for better visibility
            edgecolor="black",  # Black edges for contrast
            linewidth=2.5,  # Thicker edges
            facecolor="lightblue"
          )
          ax.add_collection3d(small_cube)

  # Draw grid lines on the surface for extra definition
  # Front face (z=0) - solid grid lines
  for x in range(0, L + 1, small_size):
    ax.plot([x, x], [0, L], [0, 0], 'k-', linewidth=1.2, alpha=0.8)
  for y in range(0, L + 1, small_size):
    ax.plot([0, L], [y, y], [0, 0], 'k-', linewidth=1.2, alpha=0.8)

  # Back face (z=L) - dashed grid lines
  for x in range(0, L + 1, small_size):
    ax.plot([x, x], [0, L], [L, L], 'k--', linewidth=1.0, alpha=0.5)
  for y in range(0, L + 1, small_size):
    ax.plot([0, L], [y, y], [L, L], 'k--', linewidth=1.0, alpha=0.5)

  # Top face (y=L) - grid lines
  for x in range(0, L + 1, small_size):
    ax.plot([x, x], [L, L], [0, L], 'k-', linewidth=0.8, alpha=0.5)
  for z in range(0, L + 1, small_size):
    ax.plot([0, L], [L, L], [z, z], 'k-', linewidth=0.8, alpha=0.5)

  # Right face (x=L) - grid lines
  for y in range(0, L + 1, small_size):
    ax.plot([L, L], [y, y], [0, L], 'k-', linewidth=0.8, alpha=0.5)
  for z in range(0, L + 1, small_size):
    ax.plot([L, L], [0, L], [z, z], 'k-', linewidth=0.8, alpha=0.5)

  # Draw diagonal lines on the front face for better visual structure
  ax.plot([0, L], [0, L], [0, 0], 'r-', linewidth=2.0, alpha=0.6)
  ax.plot([0, L], [L, 0], [0, 0], 'r-', linewidth=2.0, alpha=0.6)

  # Highlight the outer edges with even thicker lines
  outer_edges = [
    # Front face edges
    [[0, 0, 0], [L, 0, 0]],
    [[L, 0, 0], [L, L, 0]],
    [[L, L, 0], [0, L, 0]],
    [[0, L, 0], [0, 0, 0]],
    # Back face edges
    [[0, 0, L], [L, 0, L]],
    [[L, 0, L], [L, L, L]],
    [[L, L, L], [0, L, L]],
    [[0, L, L], [0, 0, L]],
    # Connecting edges
    [[0, 0, 0], [0, 0, L]],
    [[L, 0, 0], [L, 0, L]],
    [[L, L, 0], [L, L, L]],
    [[0, L, 0], [0, L, L]]
  ]

  for edge in outer_edges:
    x_vals = [edge[0][0], edge[1][0]]
    y_vals = [edge[0][1], edge[1][1]]
    z_vals = [edge[0][2], edge[1][2]]
    ax.plot(x_vals, y_vals, z_vals, 'black', linewidth=3.5, alpha=0.9)

  ax.set_xlim(0, L)
  ax.set_ylim(0, L)
  ax.set_zlim(0, L)

  ax.set_xlabel("x")
  ax.set_ylabel("y")
  ax.set_zlabel("z")

  ax.set_box_aspect([1, 1, 1])

  # Adjust viewing angle for better visibility
  ax.view_init(elev=25, azim=45)

  plt.tight_layout()
  plt.savefig("/home/washindeiru/Downloads/cells_3d_v2.png")
  plt.show()