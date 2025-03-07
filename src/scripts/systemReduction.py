import numpy as np
import re
import sys
from fractions import Fraction
from math import lcm

def read_mixed_delimiter_file(file_path):
    """
    Reads a matrix from a text file with mixed delimiters.
    
    Parameters:
        file_path (str): Path to the file.
    
    Returns:
        tuple: (2D NumPy array, Column labels, Row labels)
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Regular expression to split on spaces, commas, semicolons, or tabs
    split_pattern = r'[,\s;]+'

    # Read column labels from the first line
    column_labels = re.split(split_pattern, lines[0].strip())

    # Read numerical data
    data = [list(map(float, re.split(split_pattern, line.strip()))) for line in lines[1:]]

    # Generate row labels as R1, R2, ...
    row_labels = [f"R{i+1}" for i in range(len(data))]

    return np.array(data, dtype=float), column_labels, row_labels

def compute_label_scaling_factors(row_combinations):
    """
    Computes scaling factors to clear fractions from row labels.

    Parameters:
        row_combinations (numpy.ndarray): The row tracking matrix.

    Returns:
        list: Scaling factors for each row.
    """
    m = row_combinations.shape[0]
    scale_factors = np.ones(m, dtype=int)

    for i in range(m):
        # Convert row to fractions and get LCM of denominators
        fractions = [Fraction.from_float(val).limit_denominator(1000) for val in row_combinations[i, :]]
        denominators = [frac.denominator for frac in fractions if frac.numerator != 0]

        if denominators:
            scale_factors[i] = lcm(*denominators)

    return scale_factors

def rref_without_column_pivoting(file_path):
    """
    Reads a matrix from a file, performs RREF without column pivoting,
    keeping the first row as column labels and updating row labels with operations.
    Ensures diagonal terms are negative (-1).

    Parameters:
        file_path (str): Path to the text file containing the matrix.

    Returns:
        tuple: (RREF matrix as NumPy array, Updated column labels, Updated row labels, Rank)
    """
    # Read file with flexible delimiter handling
    A, column_labels, row_labels = read_mixed_delimiter_file(file_path)
    
    m, n = A.shape

    # Initialize row labels tracking
    row_combinations = np.eye(m, dtype=float)  # Identity matrix to track row transformations

    row = 0  # Current row to process
    for col in range(n):
        if row >= m:
            break

        # Find best pivot in this column (no column swaps)
        best_row = np.argmax(np.abs(A[row:, col])) + row  # Find row with largest value in this column

        # Skip if the pivot is zero (column is already reduced)
        if A[best_row, col] == 0:
            continue

        # Swap rows to move best pivot into position
        if best_row != row:
            A[[row, best_row], :] = A[[best_row, row], :]
            row_combinations[[row, best_row], :] = row_combinations[[best_row, row], :]  # Swap row tracking matrix

        # Normalize the pivot row
        pivot = A[row, col]
        A[row, :] /= pivot
        row_combinations[row, :] /= pivot

        # Eliminate other entries in this column
        for i in range(m):
            if i != row:
                factor = A[i, col]
                A[i, :] -= factor * A[row, :]
                row_combinations[i, :] -= factor * row_combinations[row, :]

        row += 1  # Move to the next row

    # Detect rank by counting nonzero rows
    rank = np.count_nonzero(np.any(A != 0, axis=1))

    # Move zero rows to the bottom
    zero_rows = np.where(~np.any(A, axis=1))[0]  # Indices of zero rows
    nonzero_rows = np.where(np.any(A, axis=1))[0]  # Indices of nonzero rows

    if len(zero_rows) > 0:
        A = np.vstack([A[nonzero_rows, :], A[zero_rows, :]])
        row_combinations = np.vstack([row_combinations[nonzero_rows, :], row_combinations[zero_rows, :]])

    # Ensure all diagonal elements are -1
    for i in range(min(m, n)):
        if A[i, i] > 0:
            A[i, :] *= -1
            row_combinations[i, :] *= -1

    # Compute scaling factors for row labels
    scale_factors = compute_label_scaling_factors(row_combinations)

    # Apply scaling factors to both the row labels and the matrix itself
    A *= scale_factors[:, np.newaxis]
    row_combinations *= scale_factors[:, np.newaxis]

    # Convert row tracking matrix into readable labels
    new_row_labels = []
    for i in range(m):
        label_terms = [
            f"{int(coeff)}*{original_label}" if coeff != 1 else f"{original_label}"
            for coeff, original_label in zip(row_combinations[i], row_labels) if abs(coeff) > 1e-10
        ]
        new_label = " + ".join(label_terms)
        new_row_labels.append(new_label if new_label else "0")

    return A, column_labels, new_row_labels, rank


def print_matrix(A, row_labels, column_labels):
    """
    Prints the matrix in a nicely formatted way.
    
    Parameters:
        A (numpy.ndarray): The matrix to print.
        row_labels (list): Labels for each row.
        column_labels (list): Labels for each column.
    """
    col_widths = [max(len(label), 10) for label in column_labels]
    row_label_width = max(len(label) for label in row_labels)

    # Print column headers
    print(" " * (row_label_width + 2) + " ".join(label.ljust(w) for label, w in zip(column_labels, col_widths)))

    # Print rows
    for row_label, row in zip(row_labels, A):
        row_values = " ".join(f"{val:10.5f}" for val in row)
        print(f"{row_label.ljust(row_label_width)}  {row_values}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python rref_pivot.py <matrix_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    rref_matrix, updated_labels, updated_row_labels, matrix_rank = rref_without_column_pivoting(file_path)

    print_matrix( rref_matrix, updated_row_labels, updated_labels )

    print("\nMatrix Rank:", matrix_rank)

    secondarySpecies = updated_labels[:matrix_rank]
    primarySpecies = updated_labels[matrix_rank:]

    print( "primary species: ", primarySpecies )
    print( "secondary species: ", secondarySpecies )