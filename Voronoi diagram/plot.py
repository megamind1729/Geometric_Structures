import matplotlib.pyplot as plt
import numpy as np

# Function to return endpoints of the truncated edge
def truncate(line, min_x, max_x, min_y, max_y):
    parts = line.split()
    assert len(parts) == 6  # The line should have 6 parts
    u = float(parts[4])
    v = float(parts[5])
    assert (u != 0 or v != 0)
    assert parts[0] != 'null'    # The x_coordinate of first point should not be null
    x1 = float(parts[0])
    if (parts[1] == 'null'):    # If the y_coordinate of first point is null, then the edge is vertical from positive infinity
        assert(u == 0)
        assert(x1 >= min_x and x1 <= max_x)
        if (parts[3] == 'null'):    # If the line does not end
            assert(parts[2] == 'null')
            return x1, max_y, x1, min_y
        else:
            y2 = float(parts[3])
            return x1, max_y, x1, max(float(parts[3]), min_y)
    else:
        y1 = float(parts[1])
        if (parts[2] == 'null'):    # If half-edge does not end
            assert(parts[3] == 'null')
            return truncate_half_edge(x1, y1, u, v, min_x, min_y, max_x, max_y)
        else:   # If half-edge ends, that is, finite length edge
            x2, y2 = float(parts[2]), float(parts[3])
            return truncate_segment(x1, y1, x2, y2, min_x, min_y, max_x, max_y)

def truncate_half_edge(x1, y1, u, v, min_x, min_y, max_x, max_y):
    if (u == 0):    # If the half-edge is vertical
        if(x1 < min_x or x1 > max_x or y1 < min_y):    # If the half-edge is outside the bounding box
            return None
        else:
            return x1, min(max_y, y1), x1, min_y
    elif(v == 0):   # If the half-edge is horizontal (can be towards both sides)
        if(y1 < min_y or y1 > max_y):
            return None
        elif(y1 > max_y and u > 0):
            return None
        elif(y1 < min_y and u < 0):
            return None
        else:
            if(u > 0):
                return max(min_x, x1), y1, max_x, y1
            else:
                return min(max_x, x1), y1, min_x, y1
    else:   # If the half-edge is neither vertical nor horizontal
        if((u > 0 and x1 > max_x) or (u < 0 and x1 < min_x) or (v > 0 and y1 > max_y) or (v < 0 and y1 < min_y)):
            return None
        elif(u > 0):
            x1, y1, x2, y2 = truncate_segment(x1, y1, max_x, y1 + (max_x - x1) * v / u, min_x, min_y, max_x, max_y)
        else:
            x1, y1, x2, y2 = truncate_segment(x1, y1, min_x, y1 + (min_x - x1) * v / u, min_x, min_y, max_x, max_y)
        return x1, y1, x2, y2

def truncate_segment(x1, y1, x2, y2, min_x, min_y, max_x, max_y):
    # Check if the line segment is completely outside the bounding box
    if (x1 < min_x and x2 < min_x) or (x1 > max_x and x2 > max_x) or (y1 < min_y and y2 < min_y) or (y1 > max_y and y2 > max_y):
        return None
    
    # Truncate the line segment to the bounding box
    if x1 != x2:
        m = (y2 - y1) / (x2 - x1)
        if x1 < min_x:
            y1 = y1 + m * (min_x - x1)
            x1 = min_x
        if x2 < min_x:
            y2 = y2 + m * (min_x - x2)
            x2 = min_x
        if x1 > max_x:
            y1 = y1 + m * (max_x - x1)
            x1 = max_x
        if x2 > max_x:
            y2 = y2 + m * (max_x - x2)
            x2 = max_x
    else:
        y1 = max(y1, min_y)
        y2 = max(y2, min_y)
        y1 = min(y1, max_y)
        y2 = min(y2, max_y)

    if y1 != y2:
        m = (x2 - x1) / (y2 - y1)
        if y1 < min_y:
            x1 = x1 + m * (min_y - y1)
            y1 = min_y
        if y2 < min_y:
            x2 = x2 + m * (min_y - y2)
            y2 = min_y
        if y1 > max_y:
            x1 = x1 + m * (max_y - y1)
            y1 = max_y
        if y2 > max_y:
            x2 = x2 + m * (max_y - y2)
            y2 = max_y
    else:
        x1 = max(x1, min_x)
        x2 = max(x2, min_x)
        x1 = min(x1, max_x)
        x2 = min(x2, max_x)

    return x1, y1, x2, y2


# Read points from file
points = []
with open('input.txt', 'r') as file:
    n = int(file.readline())
    for _ in range(n):
        x, y = map(float, file.readline().split())
        points.append((x, y))

# Separate the points into X and Y coordinates for plotting
x_coords, y_coords = zip(*points)

# Find the dimensions of the bounding box
padding = 20
min_x = min(x_coords) - padding
max_x = max(x_coords) + padding
min_y = min(y_coords) - padding
max_y = max(y_coords) + padding

# Create the plot
plt.scatter(x_coords, y_coords, color='red')

# Plot the bounding box
bbox_edges = [((min_x, min_y), (max_x, min_y)),
              ((max_x, min_y), (max_x, max_y)),
              ((max_x, max_y), (min_x, max_y)),
              ((min_x, max_y), (min_x, min_y))]
for edge in bbox_edges:
    plt.plot(*zip(*edge), color='green', linestyle='dashed')

# Set the limits of the plot to the bounding box dimensions
margin = 10
plt.axis('square')
plt.xlim(min_x - margin, max_x + margin)
plt.ylim(min_y - margin, max_y + margin)

# Read half edges from file and store them in a list
edges = []
with open('edge.txt', 'r') as file:
    for line in file:
        tuple = truncate(line, min_x, max_x, min_y, max_y) 
        if(tuple):
            x1, y1, x2, y2 = tuple
            edges.append(((x1, y1), (x2, y2)))

# Plot the edges
for edge in edges:
    plt.plot(*zip(*edge), color='blue')

# Show the plot
plt.show()