from os import error
from random import randrange
import numpy as np
import sympy as sp

"""
Briad Word Style:
just have them as integers...
    s1, s2, s3, = 1, 2, 3

Can write a method to communicate with sage if needed

"""

class Braid():
    """
    """
    def __init__(self, n, *ops):
        """
        """
        # Establish Group/No. of Strands
        self.braid_group = n

        # BraidWord and tracked strands developed at the same time
        self.braid_word = []
        for i in ops:
            # Check for invalid braid operation in selected group
            if abs(i) >= n:
                raise ValueError("Braid operations cannot exceed the containing braid group")
            # Build braid word
            self.braid_word.append(i)

        self.strand_labels = self.track_strands()

    def track_strands(self):
        """
        This calculates the labels of the strands when labelled from the 'bottom' of the braid up.
        """
        # Initial Positions (AT BOTTOM)
        label_positions = [i + 1 for i in range(self.braid_group)]
        strand_labels = []
        # Going through the sigmas backwards
        for op in reversed(self.braid_word):
            index = abs(op) - 1
            # Grab undercrossing string
            if np.sign(op) == -1:
                strand_labels.append(label_positions[index + 1])
            else:
                strand_labels.append(label_positions[index])
            # Swap
            label_positions[index], label_positions[index + 1] = label_positions[index + 1], label_positions[index]
        
        self.end_labels = label_positions

        strand_labels.reverse()
        return strand_labels

    def __str__(self):
        """
        Have this return with s's as the sigmas....
        """
        return f"Braid: {self.braid_word}\nLabels: {self.strand_labels}"

def alex_poly(braid, k):
    """
    Needs the number of caps and cups to calculate this stage....
        n-braid == (r + 2*k)-braid
    
    This cannot neccesarily take directly from the braid class object as the labels will be incorrect
        --> No y label for the first strand.


    What we need to do is establish the labelling truly for the strands in the braid
    Can we do this from the info given
    """
    # Grab label list
    label_list = braid.strand_labels

    n = braid.braid_group

    end_pos = braid.end_labels
    eq = []
    visited = []
    for i in range(n):
        same = []
        looped = False
        bot_list = True
        while i not in visited or bot_list == False:
            if bot_list:
                same.append(i + 1)
                visited.append(i)

            # If we have done a cap or round the back
            if not looped:
                if i < n - 2 * k:            
                    bot_list = not bot_list
                elif n % 2 == i % 2:
                    i += 1
                elif n % 2 != i % 2:
                    i -= 1
                looped = True
                
            # Straight Across
            else:
                if bot_list:
                    i = end_pos.index(i + 1)
                else:
                    i = end_pos[i] - 1
                bot_list = not bot_list
                looped = False
        
        if len(same) > 0:
            eq.append(same)

    # Using establlished above equivalences to relabel label list
    for l in eq:
        label_list = [l[0] if x in l else x for x in label_list]
    
    print(label_list)

    # Get reduced burau matrix using the tracked strands
    mat = red_burau(braid, label_list)

    r = n - 2*k

    # Perform nesseccery matrix manipulations
    x = sp.symbols("x")
    # Sub x from first r-1 diag entries
    for i in range(r - 1):
        mat[i, i] -= x
    
    for j in range(r - 1, braid.braid_group - 2, 2):
        mat.col_del(j)
    
    for l in range(r, braid.braid_group - 1, 2):
        mat.row_del(l)

    # sp.pprint(mat)

    # Determinant
    print("\nDETERMINANT:  ", mat.det())

def red_burau(braid, label_list):
    """
    implement each

    returns reduced burau matrix for current braid
    """
    # Initial matrix that will be multiplied by to create final braid burau matrix
    mat = sp.eye(braid.braid_group)

    # Going through each operation
    for i in range(len(braid.braid_word)):
        # Establish sigma and its corresponding label
        op = braid.braid_word[i]
        if label_list[i] == 1:
            label = sp.symbols("y")
        else:
            label = sp.symbols("t" + str(label_list[i]))

        # idnetity matrix to be turned into burau
        burau = sp.eye(braid.braid_group)

        # Inverse or not check...
        row = abs(op) - 1
        if np.sign(op) == 1:
            if row != 0:
                burau[row, row - 1] = label
            burau[row, row] = -label
            burau[row, row + 1] = 1
        else:
            if row != 0:
                burau[row, row - 1] = -1
            burau[row, row] = -label**-1
            burau[row, row + 1] = label**-1

        # Mulitply each time
        mat = mat*burau
    
    # Delete last row and column
    mat.row_del(braid.braid_group - 1)
    mat.col_del(braid.braid_group - 1)

    sp.pprint(mat)
    return mat

test = Braid(5, 3, 2, 2, 4, -1, -1, -2, -3, -4)
alex_poly(test, 1)

