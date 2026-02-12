#!/usr/bin/env python3
"""Create a simple wall (box) STL file for testing particle-wall contact."""

import struct

def create_box_stl(filename, size_x, size_y, size_z):
    """Create a box (cuboid) STL file centered at origin."""

    # Box vertices (centered at origin)
    hx, hy, hz = size_x / 2, size_y / 2, size_z / 2
    vertices = [
        (-hx, -hy, -hz),  # 0
        ( hx, -hy, -hz),  # 1
        ( hx,  hy, -hz),  # 2
        (-hx,  hy, -hz),  # 3
        (-hx, -hy,  hz),  # 4
        ( hx, -hy,  hz),  # 5
        ( hx,  hy,  hz),  # 6
        (-hx,  hy,  hz),  # 7
    ]

    # Define 12 triangles (2 per face)
    # Format: (v0, v1, v2) - vertices in counter-clockwise order when looking from outside
    triangles = [
        # Bottom face (z = -hz)
        (0, 2, 1), (0, 3, 2),
        # Top face (z = hz)
        (4, 5, 6), (4, 6, 7),
        # Front face (y = -hy)
        (0, 1, 5), (0, 5, 4),
        # Back face (y = hy)
        (2, 3, 7), (2, 7, 6),
        # Left face (x = -hx)
        (0, 4, 7), (0, 7, 3),
        # Right face (x = hx)
        (1, 2, 6), (1, 6, 5),
    ]

    def compute_normal(v0, v1, v2):
        """Compute the normal vector for a triangle."""
        # Edge vectors
        e1 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
        e2 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])

        # Cross product
        nx = e1[1] * e2[2] - e1[2] * e2[1]
        ny = e1[2] * e2[0] - e1[0] * e2[2]
        nz = e1[0] * e2[1] - e1[1] * e2[0]

        # Normalize
        length = (nx*nx + ny*ny + nz*nz) ** 0.5
        if length > 0:
            nx, ny, nz = nx/length, ny/length, nz/length

        return (nx, ny, nz)

    # Write binary STL
    with open(filename, 'wb') as f:
        # 80-byte header
        header = b'Binary STL - Wall Box for ZDEM-X' + b'\0' * (80 - 32)
        f.write(header)

        # Number of triangles
        f.write(struct.pack('<I', len(triangles)))

        # Write each triangle
        for tri in triangles:
            v0, v1, v2 = [vertices[i] for i in tri]
            normal = compute_normal(v0, v1, v2)

            # Normal vector
            f.write(struct.pack('<fff', *normal))
            # Vertex 1
            f.write(struct.pack('<fff', *v0))
            # Vertex 2
            f.write(struct.pack('<fff', *v1))
            # Vertex 3
            f.write(struct.pack('<fff', *v2))
            # Attribute byte count (unused)
            f.write(struct.pack('<H', 0))

    print(f"Created {filename} with {len(triangles)} triangles")
    print(f"Box size: {size_x} x {size_y} x {size_z}")

def create_plane_stl(filename, size_x, size_y):
    """Create a flat plane (2 triangles) STL file centered at origin in XY plane."""

    hx, hy = size_x / 2, size_y / 2

    # Plane vertices (z = 0)
    vertices = [
        (-hx, -hy, 0.0),  # 0
        ( hx, -hy, 0.0),  # 1
        ( hx,  hy, 0.0),  # 2
        (-hx,  hy, 0.0),  # 3
    ]

    # 2 triangles, normal pointing up (positive z)
    triangles = [
        (0, 1, 2),
        (0, 2, 3),
    ]

    normal = (0.0, 0.0, 1.0)

    # Write binary STL
    with open(filename, 'wb') as f:
        header = b'Binary STL - Plane for ZDEM-X' + b'\0' * (80 - 29)
        f.write(header)
        f.write(struct.pack('<I', len(triangles)))

        for tri in triangles:
            v0, v1, v2 = [vertices[i] for i in tri]
            f.write(struct.pack('<fff', *normal))
            f.write(struct.pack('<fff', *v0))
            f.write(struct.pack('<fff', *v1))
            f.write(struct.pack('<fff', *v2))
            f.write(struct.pack('<H', 0))

    print(f"Created {filename} with {len(triangles)} triangles")
    print(f"Plane size: {size_x} x {size_y}")

if __name__ == "__main__":
    # Create a thin box as wall (2m x 2m x 0.1m)
    create_box_stl("wall_box.stl", 2.0, 2.0, 0.1)

    # Also create a simple plane
    create_plane_stl("wall_plane.stl", 4.0, 4.0)
