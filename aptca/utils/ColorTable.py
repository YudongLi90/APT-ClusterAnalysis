"""
A color table and related methods to distinguish neighboring items. 
"""
import numpy as np
from itertools import cycle

COLOR_TABLE = [[0,0,0], [1,0,103], [213, 255, 0], [255, 0, 86], [158, 0, 142],
                       [14, 76, 161], [255, 229, 2], [0, 95, 57], [0, 255, 0], [149, 0, 58],
                       [255, 147, 126], [164, 36, 0], [0, 21, 68], [145, 208, 203], [98, 14, 0]
                       , [107, 104, 130], [0, 0, 255], [0, 125, 181], [106, 130, 108]
                       , [0, 174, 126], [194, 140, 159], [190, 153, 112], [0, 143, 156]
                       , [95, 173, 78], [255, 0, 0], [255, 0, 246], [255, 2, 157], [104, 61, 59]
                       , [255, 116, 163], [150, 138, 232], [152, 255, 82], [167, 87, 64]
                       , [1, 255, 254], [255, 238, 232], [254, 137, 0], [189, 198, 255]
                       , [1, 208, 255], [187, 136, 0], [117, 68, 177], [165, 255, 210]
                       , [255, 166, 254], [119, 77, 0], [122, 71, 130], [38, 52, 0]
                       , [0, 71, 84], [67, 0, 44], [181, 0, 255], [255, 177, 103]
                       , [255, 219, 102], [144, 251, 146], [126, 45, 210], [189, 211, 147]
                       , [229, 111, 254], [222, 255, 116], [0, 255, 120], [0, 155, 255]
                       , [0, 100, 1], [0, 118, 255], [133, 169, 0], [0, 185, 23]
                       , [120, 130, 49], [0, 255, 198], [255, 110, 65], [232, 94, 190]]


def gen_color_table(style='int'):
    """
    Generate a color table. Three color code styles can be chosed: int, float, hex.

    style:
        int:
            [r, g, b], where all color range from 0 to 255;
        float:
            [r, g, b], where all color range from 0.0 to 1.0;
        hex:
            "BBBBBB", where color is represented by a hex number string. 

    return a color table (list) in specified style.
    """
    def rgb2hex(color_rgb):
        """
        A little helper function to convert RGB color to hex color.
        """
        r = color_rgb[0]
        g = color_rgb[1]
        b = color_rgb[2]
        hex_color = "{:02x}{:02x}{:02x}".format(r,g,b)
        return hex_color

    def rgb2float(color_rgb):
        r = color_rgb[0]
        g = color_rgb[1]
        b = color_rgb[2]
        return [r/255.0, g/255.0, b/255.0]

    if style == 'int':
        color_table = COLOR_TABLE
    elif style == 'float':
        color_table = list(map(rgb2float, COLOR_TABLE))
    elif style == 'hex':
        color_table = list(map(rgb2hex, COLOR_TABLE))

    return color_table

# Test
if __name__ == '__main__':
    print(len(COLOR_TABLE))
    print(gen_color_table('hex'))