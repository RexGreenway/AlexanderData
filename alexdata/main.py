
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy.simplify.simplify import simplify

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

    def draw(self, style = "comp", linewidth = 3, gap_size = 3, color = "rainbow", save = False):
        """
        Uses MatPLotLib to draw a braid.

        style : "comp" or "ext"
        linewidth : int
        gap_size : int
        color : "rainbow" or all fixed colour from 
            'b' as blue
            'g' as green
            'r' as red
            'c' as cyan
            'm' as magenta
            'y' as yellow
            'k' as black
            'w' as white
        """
        n = self.braid_group
        braid = self.braid_word

        fig, ax = plt.subplots(figsize = (4, 8))

        t = np.arange(0, 1.05, 0.05)

        # line break
        lb = np.array([[np.NaN, np.NaN]])

        x = 0
        layer = []

        strands = [np.empty((0, 2)) for _ in range(n)]

        # Start positons
        for s in range(1, n + 1):
            S = np.array([[s], [0]]).T
            strands[s - 1] = np.concatenate((strands[s - 1], S))

        # Calc crossing points
        for count, op in enumerate(braid):
            i = abs(op)

            if style == "comp":
                if i in layer or i + 1 in layer:
                    layer = []
                    x += 2
                layer.append(i)
                layer.append(i + 1)

            # Crossing Matrices
            R1 = np.array([0.5*t**2 + i, -t - x]).T
            R2 = np.array([-0.5*t**2 + 1 + i, t - 2 - x]).T
            L1 = np.array([-0.5*t**2 + 1 + i, -t - x]).T
            L2 = np.array([0.5*t**2 + i, t - 2 - x]).T

            if np.sign(op) == 1:    # overcrossing
                strands[i - 1] = np.concatenate((strands[i - 1], R1, np.flipud(R2)))
                strands[i] = np.concatenate((strands[i], L1[:-gap_size, :], lb, np.flipud(L2[:-gap_size, :])))
            elif np.sign(op) == -1:     # undercrossing
                strands[i - 1] = np.concatenate((strands[i - 1], R1[:-gap_size, :], lb, np.flipud(R2[:-gap_size, :])))
                strands[i] = np.concatenate((strands[i], L1, np.flipud(L2)))
            
            # swap lines for each operation
            strands[i - 1], strands[i] = strands[i], strands[i - 1]

            if style == "ext" or count + 1 == len(braid):
                x += 2

        # End positons
        for s in range(1, n + 1):
            S = np.array([[s], [-(x)]]).T
            strands[s - 1] = np.concatenate((strands[s - 1], S))
            if color == "rainbow":
                ax.plot(strands[s - 1][:, 0], strands[s - 1][:, 1], lw = linewidth)
            else:
                ax.plot(strands[s - 1][:, 0], strands[s - 1][:, 1], lw = linewidth, c = color)

        # Figure Details
        fig.suptitle("Braid:   " + str(braid))
        ax.axis("off")
        ax.set_aspect(2*n / x)
        plt.tight_layout()

        # SAVE CHECK
        if save:
            plt.savefig("test.svg")
        else:
            plt.show()

    def __str__(self):
        """
        Have this return with s's as the sigmas....
        """
        return f"Braid: {self.braid_word}\nLabels: {self.undercrossing_labels}"


## Braid Kernel ##
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

        self.eq = eq
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
            if np.sign(op) == -1:
                if row != 0:
                    burau[row, row - 1] = label
                burau[row, row] = -label
                burau[row, row + 1] = 1
            else:
                if row != 0:
                    burau[row, row - 1] = 1
                burau[row, row] = -label**-1
                burau[row, row + 1] = label**-1
            
            sp.pprint(burau)
            print("\n")

            # Mulitply each time
            mat = mat*burau
        
        sp.pprint(mat.applyfunc(simplify))
        print("\n")

        # Delete last row and column
        mat.row_del(self.braid_group - 1)
        mat.col_del(self.braid_group - 1)

        sp.pprint(mat.applyfunc(simplify))

        return mat
    
    def alexander_polynomial(self):
        """
        """

        M = self.reduced_burau()

        n = self.braid_group
        k = self.caps
        r = n - 2*k

        # remove x
        x = sp.symbols("x")
        for i in range(r - 1):
            M[i, i] -= x

        # remove columns
        M = sp.Matrix([[M[:, :(r-1)], M[:, r:(n-2):2], M[:, (n-2):]]])
        # remove rows
        M = sp.Matrix([M[:r, :], M[r+1:n-1:2, :], M[n-1:, :]])

        print("\n")
        sp.pprint(M.applyfunc(simplify))
        print("\n")
        # print(M.det())
        return M.det()
    
    def alexander_data(self):
        """
        """
        det = self.alexander_polynomial()
        det_simp = sp.simplify(det)
        sp.pprint(det_simp)

        x = sp.symbols("x")
        y = sp.symbols("y")
        for i in det_simp.free_symbols:
            if i != x and i != y:
                t = sp.symbols(str(i))
                s = sp.symbols("s" + str(i)[-1])
                det_simp = det_simp.subs(t, s**2)
        
        sp.pprint(det_simp)
        # This works at leat for this specific det but the method is GREEDY.
        sp.pprint(sp.collect(det_simp, [x**2*y**2, x**2*y, x*y**2, x*y, x, y]))


    def draw(self, style = "ext", linewidth = 3, gap_size = 5, color = "rainbow", save = False):
        """
        Use MatPLotLib? or Bokeh? to draw a braid.
        style : "comp" or "ext"
        linewidth : int
        gap_size : int
        color : 
            'b' as blue
            'g' as green
            'r' as red
            'c' as cyan
            'm' as magenta
            'y' as yellow
            'k' as black
            'w' as white
        """
        n = self.braid_group
        braid = self.braid_word
        caps = self.caps

        # Est. figure
        fig, ax = plt.subplots(figsize = (6, 8))

        # Initial Values
        t = np.arange(0, 1.05, 0.05)
        x = 0
        layer = []

        # line break
        lb = np.array([[np.NaN, np.NaN]])

        # Strands
        strands = [np.empty((0, 2)) for _ in range(n)]

        # Start positons
        for s in range(1, n + 1):
            S = np.array([[s], [0]]).T
            strands[s - 1] = np.concatenate((strands[s - 1], S))

        # Calc crossing points
        for count, op in enumerate(braid):
            i = abs(op)

            if style == "comp":
                if i in layer or i + 1 in layer:
                    layer = []
                    x += 2
                layer.append(i)
                layer.append(i + 1)

            # Crossing Matrices
            R1 = np.array([0.5*t**2 + i, -t - x]).T
            R2 = np.array([-0.5*t**2 + 1 + i, t - 2 - x]).T
            L1 = np.array([-0.5*t**2 + 1 + i, -t - x]).T
            L2 = np.array([0.5*t**2 + i, t - 2 - x]).T

            if np.sign(op) == 1:    # overcrossing
                strands[i - 1] = np.concatenate((strands[i - 1], R1, np.flipud(R2)))
                strands[i] = np.concatenate((strands[i], L1[:-gap_size, :], lb, np.flipud(L2[:-gap_size, :])))
            elif np.sign(op) == -1:     # undercrossing
                strands[i - 1] = np.concatenate((strands[i - 1], R1[:-gap_size, :], lb, np.flipud(R2[:-gap_size, :])))
                strands[i] = np.concatenate((strands[i], L1, np.flipud(L2)))

            # swap lines for each operation
            strands[i - 1], strands[i] = strands[i], strands[i - 1]

            if style == "ext" or count + 1 == len(braid):
                x += 2

        # Enpoints, Loopping + caps
        for k in range(1, n + 1):
            # Loops
            if k <= n - 2*caps:
                # Initial values..
                loop_x = np.arange(0, k + 0.05, 0.05)
                loop_y = np.sqrt(k**2 - loop_x**2)
                # Form underloop
                u_x = np.concatenate([loop_x[::-1], -loop_x])
                u_y = np.concatenate([-loop_y[::-1], -loop_y])
                B = np.array([u_x, u_y - x]).T
                # Form overloop
                o_x = np.concatenate([-loop_x[::-1], loop_x])
                o_y = np.concatenate([loop_y[::-1], loop_y])
                T = np.array([o_x, o_y]).T
                # Add to strand matrixes
                strands[k - 1] = np.concatenate([strands[k - 1], B, T])
            # Caps
            elif (k % 2) == (n % 2):
                # Initial values..
                loop_x = np.arange(0, 0.55, 0.05)
                loop_y = np.sqrt(0.5**2 - loop_x**2)
                # Form underloop
                u_x = np.concatenate([loop_x[::-1], -loop_x])
                u_y = np.concatenate([-loop_y[::-1], -loop_y])
                B = np.array([u_x + k - 0.5, u_y - x]).T
                # Form overloop
                o_x = np.concatenate([-loop_x[::-1], loop_x])
                o_y = np.concatenate([loop_y[::-1], loop_y])
                T = np.array([o_x + k - 0.5, o_y]).T
                # Add to strand matrixes
                strands[k - 1] = np.concatenate([strands[k - 1], B, lb, T])
            # Endpoints
            else:
                S = np.array([[k], [-x]]).T
                strands[k - 1] = np.concatenate((strands[k - 1], S))

        # PLOTTING (plot sam,e colours in final version using available equivalence infomation!!!!!)
        if color == "rainbow":
            for group in range(len(self.eq)):
                for s in self.eq[group]:
                    ax.plot(strands[s - 1][:, 0], strands[s - 1][:, 1], lw = linewidth, c = "C" + str(group))
        else:
            for s in range(n):
                ax.plot(strands[s][:, 0], strands[s][:, 1], lw = linewidth, c = color)

        # Figure Details
        fig.suptitle("Braid:   " + str(braid) + "\nNumber of Caps:   " + str(caps))
        ax.axis("off")
        ax.set_aspect(2*n / x)
        plt.tight_layout()
    
        # SAVE CHECK
        if save:
            plt.savefig("test.svg")
        else:
            plt.show()

# test = Braid(5, 3, 2, 2, -4, -1, -1, -2, -3, -4)
# test.draw()
# alex_poly(test, 1)

new = Braid_Kernal(5, 1, 3, 2, 2, -4, -1, -1, -2, -3, -4)
new.draw()
# new.alexander_polynomial()
new.alexander_data()
# print(new)
