from PIL import Image
import numpy as np
from random import randint
import os
import sys

dimensions = int(sys.argv[1])
R_color = int(sys.argv[2])
G_color = int(sys.argv[3])
B_color = int(sys.argv[4])

R = np.full((dimensions, dimensions), R_color)
G = np.full((dimensions, dimensions), G_color)
B = np.full((dimensions, dimensions), B_color)

# Convert the pixels into an array using numpy
array_R = np.array(R, dtype=np.uint8)
array_G = np.array(G, dtype=np.uint8)
array_B = np.array(B, dtype=np.uint8)

RGB = np.dstack([array_R, array_G, array_B])

# Use PIL to create an image from the new array of pixels
new_image = Image.fromarray(RGB)
new_image.save('base_color_map.jpg')

dir_path = os.path.dirname(os.path.realpath(__file__))

os.rename(dir_path + "\\base_color_map.jpg", dir_path + "\\generated_maps\\base_color_map.jpg")
 