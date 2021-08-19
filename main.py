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
        braid_group is the number of strands in braid.
        briad_word is the list of artin opertaors as integers.
        undercrossing_labels is the list of UNDERCROSSING strands at each operation
        top_labels is n-length list of final strand positions (AT TOP)
        bot_labels is n-length list of initial strand positions (AT BOTTOM)
        """
        # Establish Group/No. of Strands
        self.braid_group = n

        # BraidWord and tracked strands developed at the same time
        self.braid_word = []
        for i in ops:
            # Check for invalid braid operation in selected group
            if abs(i) >= n:
                raise ValueError("Braid operations cannot exceed the containing braid group.")
            # Build braid word
            self.braid_word.append(i)

        self.undercrossing_labels = self.track_strands()

    def track_strands(self):
        """
        Tracks starnd positions through the braid itself

        This calculates the labels of the strands when labelled from the 'bottom' of the braid up.
        """
        # Initial Positions (AT BOTTOM)
        label_positions = [i + 1 for i in range(self.braid_group)]
        self.bot_labels = label_positions

        undercrossing_labels = []
        # Going through artin opertaors backwards (i.e. from bottom of braid)
        for op in reversed(self.braid_word):
            index = abs(op) - 1
            # Grab undercrossing string
            if np.sign(op) == -1:
                undercrossing_labels.append(label_positions[index + 1])
            else:
                undercrossing_labels.append(label_positions[index])
            # Swap
            label_positions[index], label_positions[index + 1] = label_positions[index + 1], label_positions[index]

        # Final Positions (AT TOP)
        self.top_labels = label_positions

        undercrossing_labels.reverse()
        return undercrossing_labels

    def draw(self):
        """
        Use MatPLotLib? or Bokeh? to draw a braid.
        """
        return

    def __str__(self):
        """
        Have this return with s's as the sigmas....
        """
        return f"Braid: {self.braid_word}\nLabels: {self.undercrossing_labels}"

class Braid_Kernal(Braid):
    """
    This is a sub class of braid and provides the framework to
    """
    def __init__(self, n, k, *ops):
        """
        """
        # Check Valid k input:
        if k > n / 2:
            raise ValueError("Number of caps cannot exceed half of total number of strands.")
        self.caps = k
        super().__init__(n, *ops)

    def track_strands(self):
        """
        first calls super() method to track braid labels
        then works out which strands are the same around the kernal.
        """
        label_list = super().track_strands()

        n = self.braid_group
        k = self.caps

        end_pos = self.top_labels
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

        return label_list

    def reduced_burau(self):
        """
        """
        mat = sp.eye(self.braid_group)
        label_list = self.undercrossing_labels

        # Going through each operation
        for i in range(len(self.braid_word)):
            # Establish sigma and its corresponding label
            op = self.braid_word[i]
            if label_list[i] == 1:
                label = sp.symbols("y")
            else:
                label = sp.symbols("t" + str(label_list[i]))

            # idnetity matrix to be turned into burau
            burau = sp.eye(self.braid_group)

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
        mat.row_del(self.braid_group - 1)
        mat.col_del(self.braid_group - 1)

        sp.pprint(mat)

        return mat
    
    def alexander_polynomail(self):
        
        M = self.reduced_burau()

        n = self.braid_group
        k = self.caps

        r = n - 2*k

        S = M[:, r-1:n-1:2]


    def draw(self):
        """
        Use MatPLotLib? or Bokeh? to draw a braid.
        """
        print("draw kernal")
        super().draw()

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
    label_list = braid.undercrossing_labels

    n = braid.braid_group

    end_pos = braid.top_labels
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

    # This doesnt work becuase with deletion the index xhanges!!!!
    # for j in range(r - 1, int((n - 2) / 2)):
    #     mat.col_del(j)
    # for l in range(r, int((n - 1) / 2)):
    #     mat.row_del(l)

    mat.col_del(r-1)
    mat.col_del(r)


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

# test = Braid(5, 3, 2, 2, 4, -1, -1, -2, -3, -4)
# alex_poly(test, 1)

new = Braid_Kernal(5, 1, 3, 2, 2, 4, -1, -1, -2, -3, -4)
new.reduced_burau()
print(new)
