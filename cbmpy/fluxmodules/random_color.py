'''
This code is adapted from https://gist.github.com/adewes/5884820

Thanks to adewes for putting it online

TODO: clarify copyright - if I would not be so lazy I could have come up with it
myself
'''

import random
 
myrandom = random.Random(135) 

def get_random_color():
    x = myrandom.uniform(0,1.0)
    y = myrandom.uniform(0.7,1.0)
    z = myrandom.uniform(0.7,1.0)
    return (x, y, z)
 
def color_distance(c1,c2):
    return min(abs(c1[0] - c2[0]), 1-abs(c1[0] - c2[0])) + abs(c1[1] - c2[1]) + abs(c1[2] - c2[2])
 
def generate_new_color(existing_colors):
    """ creates a new random color (in hsv scheme) that is most different
    to the existing colors.
    """
    max_distance = None
    best_color = None
    for i in range(0,100):
        color = get_random_color()
        if not existing_colors:
            return color
        best_distance = min([color_distance(color,c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    return best_color
 