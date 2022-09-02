from PIL import Image
import numpy as np
from random import randint
import os
import sys

dimensions = int(sys.argv[1])

pixels = np.zeros((dimensions, dimensions))

for i in range(0, dimensions):
    for j in range(0, dimensions):
        pixels[i][j] = randint(0, 255)

# Convert the pixels into an array using numpy
array = np.array(pixels, dtype=np.uint8)

# Use PIL to create an image from the new array of pixels
new_image = Image.fromarray(array)
new_image.save('displacement_map.jpg')
 
dir_path = os.path.dirname(os.path.realpath(__file__))

os.rename(dir_path + "\\displacement_map.jpg", dir_path + "\\generated_maps\\displacement_map.jpg")
 