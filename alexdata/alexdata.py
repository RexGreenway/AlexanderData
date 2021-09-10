"""
AlexanderData Library
=====================
Created in progress toward completion of COMP702 - Dissertation
Module, MA Data Science and Artificial Intelligence @ University
of Liverpool.

This library was written for the purpose of calculating the 
Alexander Data, an invariant for textiles. To achieve this goal
this package allows for the creation of Braid objects - 
specifically the Braid Kernel, a representation of repeating
textile structures. From this object the user can calculate the
entire Alexander Data, or each of the internel stages piece-meal.

"""

# Dependent Libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


__all__ = [
    "Braid",
    "Braid_Kernel"
]


class Braid():
    """
    Braid class object with internel starnd tracking and drawing
    functionality.

    Attributes
    ----------
    braid_group : int
        The number of strands in the braid.
    braid_word : list
        Sequence of Artin operators that define the braid 
        (Negative values indicate undercrossing strands).
    
    undercrossing_labels : list
        List of underscrossing strands' labels for each 
        Artin operation repectively.
    
    Notes
    -----
    abcd
    """
    def __init__(self, n, *ops):
        """
        Initialises the Braid class object, including calling
        the internal tracking function to generate strand positions.
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

        self.undercrossing_labels = self._track_strands()

    def _track_strands(self):
        """
        Protected Method
        Tracks starnd positions through the given braid.

        top_labels : list
            List of 'briad_group' length, indicating strand
            names / labels present at the start, or 'top', 
            of the braid. (Tracked internally from sequence
            of operations)
        bot_labels : list
            List of 'briad_group' length, indicating strand
            names / labels present at the end, or 'bottom',
            of the Briad. (Default [1, 2, ..., n])

        Returns
        -------
        undercrossing_labels : list
            List of underscrossing strands' labels for each 
            Artin operation repectively.

        Notes
        -----
        This calculates the labels of the strands when labelled from the
        'bottom' of the braid up as is the set convention for labelling
        used by Professor Morton in his work.
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
        Draws the Braid as a MatpLotlib figure. 

        Parameters
        ----------
        style : "comp" or "ext"
            "comp" renders the image of the braid in a compact style
            with crossings parallel to one another if possible. "ext",
            for extended, shows the crossings in series.
        linewidth : int (Default = 3)
            Thickness of the strands in the figure.
        gap_size : int (Default = 3)
            Amount of space shown at crossings for undercrossing strands.
        color : str
            Multicolor strands defined by "rainbow". Single fixed colour
            for all strands can be chosen from:
                {'b': blue,
                'g': green,
                'r': red,
                'c': cyan,
                'm': magenta,
                'y': yellow,
                'k': black,
                'w': white}

        Notes
        -----
        abcd

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
        Returns descriptive information about the Braid object.
        """
        return f"Braid: {self.braid_word}\nLabels: {self.undercrossing_labels}"


## Braid Kernel ##
class Braid_Kernel(Braid):
    """
    Braid Kernel class object for the creation of braid-based kernels
    and the subsequent calculation of the Alexander Data.

    Attributes
    ----------
    (See parent class of Braid for braid specific attributes.)

    eq : list[list]
        List of equivalence classes for strand labels, detailing the
        strands that connect to one another through loops of the Kernel.

    Notes
    -----
    abcd

    """
    def __init__(self, n, k, *ops):
        """
        Initialises the Braid_Kernel class object with internal strand
        tracking through the braid and around the loops of teh Kernel. 
        """
        # Check Valid k input:
        if k > n / 2:
            raise ValueError("Number of caps cannot exceed half of total number of strands.")
        self.caps = k
        super().__init__(n, *ops)

    def _track_strands(self):
        """
        Protected Method
        Tracks strand positions through the braid and the loops of the
        braid-based Kernel.

        Returns
        -------
        label_list : list
            Adjusted list of undercrossing strands accounting for
            equivalences aross Kernel loops and caps.

        Notes
        -----
        This method first calls the parent method of track strands to
        establish crossing specific strand movement before accounting
        for Kernel loops and caps.
        
        """
        # Call parent label tracking method for undercrossing labels in braid.
        label_list = super()._track_strands()

        n = self.braid_group
        k = self.caps

        end_pos = self.top_labels
        eq = []
        visited = []
        for i in range(n):
            same = []
            looped = False
            bot_list = True
            direction = +1
            while i not in visited or bot_list == False:
                if bot_list:
                    # Append to same: strand label and direction tuple.
                    same.append((i + 1, direction))
                    visited.append(i)
                
                # If we have done a cap or round the back
                if not looped:
                    # Loops
                    if i < n - 2 * k:
                        bot_list = not bot_list
                    # Caps
                    elif n % 2 == i % 2:
                        i += 1
                        direction = -direction
                    elif n % 2 != i % 2:
                        i -= 1
                        direction = -direction
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
            label_list = [l[0][0] if x in [k[0] for k in l] else x for x in label_list]

        self.eq = eq
        return label_list

    def reduced_burau(self, print_result = True):
        """
        Produces the Reduced Burau Matrix for the Braid_Kernel object.

        Parameters
        ----------
        print_result : boolean (Default = True)
            Prints the final result of this method. Defaulted as True
            so as to be shown if mthod called individually.

        Returns
        -------
        mat : SymPy Matrix
            The reduced Burau Matrix of the current Briad_Kernel object.

        Notes
        -----
        This method performs the substitution of the first strand "t1" with
        "y" in line with the Kernel. This means that all fabric specific
        strands within this code are shown as "t2" onwards. 

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

            # Mulitply each time
            mat = mat*burau

        # Delete last row and column
        mat.row_del(self.braid_group - 1)
        mat.col_del(self.braid_group - 1)

        # Prints Reduced Burau
        if print_result:
            sp.pprint(mat.applyfunc(sp.simplify))

        return mat
    
    def alexander_polynomial(self, print_result = True):
        """
        Produces the Alexander polynomial for the Braid_Kernel

        Parameters
        ----------
        print_result : boolean (Default = True)
            Prints the final result of this method. Defaulted as True
            so as to be shown if mthod called individually.

        Returns
        -------
        det : SymPy Equation
            Determinant of the Briad Kernel in terms of "x", "y", and
            strands "t?".

        Notes
        -----
        abcd

        """

        M = self.reduced_burau(print_result = False)

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

        # Prints modified red-burau Determinant
        if print_result:
            print(M.det())

        return M.det()
    
    def alexander_data(self, print_result = True):
        """
        Produces the Alexander Data of the Braid_Kernel.

        Parameters
        ----------
        print_result : boolean (Default = True)
            Prints the final result of this method. Defaulted as True
            so as to be shown if mthod called individually.
        
        Returns
        -------

        Notes
        -----

        """
        # Gets Alexander Poly. / determinant
        det = self.alexander_polynomial(print_result = False)
        det = sp.simplify(det)

        x = sp.symbols("x")
        y = sp.symbols("y")

        n = self.braid_group
        k = self.caps
        r = n - 2*k

        U = y
        V = x

        # Finding linking no.s and constructing U and V
        # finds strands that cross the y strand
        strands = []
        prev = 0
        for index, op in enumerate(self.braid_word):
            if op == -1 and prev == -1:
                # grab strand from undercrossing labels
                strands.append(self.undercrossing_labels[index])
            prev = op

        for i, g in enumerate(self.eq):
            # skip y strand
            if g[0][0] == 1:
                continue
            a = 0
            b = 0
            for s in g:
                # If NOT capped part then add up directions. A
                if s[0] <= r:
                    a += s[1]
                # If crossing y strand.
                if s[0] in strands:
                    b += s[1]

            s = sp.symbols("s" + str(i + 1))
            U *= s**a
            V *= s**b

        print("U = ", U)
        print("V = ", V)

        # Subbing t for s^2.
        for i in det.free_symbols:
            if i != x and i != y:
                t = sp.symbols(str(i))
                s = sp.symbols("s" + str(i)[-1])
                det = det.subs(t, s**2)
        
        # PRINTS ALEX POLY with subbed si
        sp.pprint(det)
        # This works at leat for this specific det but the method is GREEDY.
        # sp.pprint(sp.collect(det, [x**2*y**2, x**2*y, x*y**2, x*y, x, y]))

        # max power is (r + k - 1)
        data = sp.zeros(r + k - 1)
        for i in range(r + k - 1):
            for j in range(r + k - 1):
                # print("\nRow:", (j + 1), "Col: ", (i + 1))
                UV = U**(j + 1)*V**(i + 1)
                # sp.pprint(UV)
                term = det.coeff(x**(i + 1)*y**(j + 1))
                # sp.pprint(term)
                if term != 0:
                    inpt = 0
                    for t in sp.Add.make_args(term):
                        inpt += t / UV
                        # sp.pprint(inpt.simplify())
                data[j, i] = inpt

        # sp.pprint(data.applyfunc(sp.simplify))
                
        return U, V

    def draw(self, style = "ext", linewidth = 3, gap_size = 5, color = "rainbow", save = False):
        """
        Draws the Braid_Kernel as a MatpLotlib figure. 

        Parameters
        ----------
        style : "comp" or "ext"
            "comp" renders the image of the braid in a compact style
            with crossings parallel to one another if possible. "ext",
            for extended, shows the crossings in series.
        linewidth : int (Default = 3)
            Thickness of the strands in the figure.
        gap_size : int (Default = 3)
            Amount of space shown at crossings for undercrossing strands.
        color : str
            Multicolor strands defined by "rainbow". Single fixed colour
            for all strands can be chosen from:
                {'b': blue,
                'g': green,
                'r': red,
                'c': cyan,
                'm': magenta,
                'y': yellow,
                'k': black,
                'w': white}

        Notes
        -----
        abcd
        
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
                    ax.plot(strands[s[0] - 1][:, 0], strands[s[0] - 1][:, 1], lw = linewidth, c = "C" + str(group))
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


if __name__ == "__main__":

    # new = Braid(5, 3, 2, 2, -4, -1, -1, -2, -3, -4)   
    new = Braid_Kernel(5, 1, 3, 2, 2, -4, -1, -1, -2, -3, -4)
    # print(new.undercrossing_labels)
    # new.alexander_polynomial()
    new.alexander_data()

    new.draw()
    # print(new)
