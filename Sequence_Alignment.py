# Information #############################################################
# Creator: Marc Olivier Bertoz
# Date of creation: 27 October 2021
# Goal of the Application: Main goal is to try to create a Dynamic Programming application that allow to perform sequence alignment, using the Streamlit package
###########################################################################

# Loading packages ########################################################
#from matplotlib.pyplot import step
import streamlit as st
import numpy as np

###########################################################################

# Functions definition ####################################################

# Create null matrix

def create_matrix(rows,columns):
    matrix_null = np.array(np.zeros((rows,columns), dtype=int))
    return matrix_null

# Fill first row and column with gap value

def fill_gap(gap_value, matrix,row_length,column_length):
    matrix_gap = matrix
    n=1
    # Filling the first column
    for i in range(column_length):
        matrix_gap[n+i,0] = matrix_gap[n+(i-1),0] + (gap_value)
    
    # Filling the first row
    for i in range(row_length):
        matrix_gap[0,n+i] = matrix_gap[0,n+(i-1)] + (gap_value)

    return matrix_gap

# Computing the scores

def compute_score(match_value,mismatch_value,gap_value,matrix,row_length,column_length, sequence1,sequence2):
    matrix_score = matrix
    n=1
    for i in range(column_length):
        for j in range(row_length):
            scores=[]
            
            # Computing the gap in horizontal
            gap_horizontal = matrix_score[j,n+i] + (gap_value)

            scores.append(gap_horizontal)

            # Computing the gap in vertical
            gap_vertical = matrix_score[n+j,i] +(gap_value)

            scores.append(gap_vertical)

            # Computing the match or mismatch

            if sequence1[i] == sequence2[j]:
                value = matrix_score[j,i] + (match_value)
                scores.append(value)
            else:
                value = matrix_score[j,i] + (mismatch_value)
                scores.append(value)
            

            matrix_score[n+j,n+i]= np.amax(scores)

    return matrix_score

# Get the Optimal value from matrix

def get_optimal_value(matrix):
    n_rows = len(matrix)
    n_cols = len(matrix[0])
    optimal_value = matrix[n_rows-1,n_cols-1]

    return optimal_value

# Get lenght of matrix
def length_matrix(matrix):
    n_rows = len(matrix)
    n_cols = len(matrix[0])
    return n_rows,n_cols

# Retrive the sequence
def post_sequence(matrix,sequence_solution,n_rows,n_cols,match_value,mismatch_value,gap_value):
        
    if n_rows == 0 or n_cols== 0:
        return sequence_solution
    
    # Computing the various score
    # Here will be followed almost the same process as for filling up the matrix,
    # but starting from the end of the matrixand going back recurvely
    
    # Declaring the statements

    # Gap in horizontal
    Gap_Horizontal = matrix[n_rows-1][n_cols] + gap_value
    
    # Gap in vertical
    Gap_Vertical = matrix[n_rows][n_cols-1] + gap_value
    
    # Match or mismatch (diagonal)
    Final_Diagonal = 0
    Pos_gap_d = matrix[n_rows-1][n_cols-1]

    if Pos_gap_d + match_value == matrix[n_rows][n_cols]:
        Final_Diagonal = Pos_gap_d + match_value
    elif Pos_gap_d + mismatch_value == matrix[n_rows][n_cols]:
        Final_Diagonal = Pos_gap_d + mismatch_value


    # Get the best value
    Best_Score = max(Gap_Horizontal,Gap_Vertical,Final_Diagonal)

    # Applying the recursion, meaning calling the same funciton within it, but using new data, 
    # as moving backwards, the matrix becomes smaller and smaller.

    if Best_Score == Gap_Horizontal:
        sequence_solution.append('Gap_Horizontal')
        return post_sequence(matrix[:n_rows][:n_cols+1],sequence_solution,n_rows-1,n_cols,match_value,mismatch_value,gap_value)

    elif Best_Score == Gap_Vertical:
        sequence_solution.append('Gap_Vertical')
        return post_sequence(matrix[:n_rows+1][:n_cols],sequence_solution,n_rows,n_cols-1,match_value,mismatch_value,gap_value)

    elif Best_Score == Final_Diagonal:
        sequence_solution.append('Align')
        return post_sequence(matrix[:n_rows][:n_cols],sequence_solution,n_rows-1,n_cols-1,match_value,mismatch_value,gap_value)

# Create new vector with correct sequence
def final_sequences(sequence1,sequence2,sequence_solution):

    for i in range(len(sequence_solution)):
        if sequence_solution[i] == "Align":
            continue
        elif sequence_solution[i] == "Gap_Vertical":
            sequence2.insert(i,'-')
        elif sequence_solution[i] == "Gap_Horizontal":
            sequence1.insert(i,'-')

    return sequence1,sequence2

###########################################################################
# Vector for computing alignement
sequence = []

# Front-End ###############################################################
st.title("Sequence Alignment - Dynamic Programming")
st.write("This is an experiment, with the main goal to create an application that allow to compute the best **global** sequence alignment of two DNA sequences using dynamic programming. " )

# Columns for the Scoring
st.header("Scoring point:")
left_score_column, middle_score_column, right_score_column = st.columns(3)


# Columns for the sequences
st.header("Sequences:")
left_sequence_column, right_sequence_column = st.columns(2)



# Display columns ###########################################################################

# Columns for scores
with left_score_column:
    match=st.number_input("Match:",step=1)
with middle_score_column:
    mismatch=st.number_input("Mismatch:", step=1)
with right_score_column:
    gap=st.number_input("Gap:", step=1)

# Columns for input
with left_sequence_column:
    first_sequence=st.text_area("Write the first sequence as in the eexample: ACTG")
with right_sequence_column:
    second_sequence=st.text_area("Write the second sequence as in the eexample: ACTG")

##############################################################################################

# Creating button
pulsante_calcolo=st.button("Compute")

# Button commands
if pulsante_calcolo:
    array1 = list(first_sequence)
    array2 = list(second_sequence)
    array1_length = len(array1)
    array2_length = len(array2)
   

    # Creating the matrix starting from a null matriz (zeros) to the score matrix
    zeros_matrix = create_matrix(array2_length+1,array1_length+1)

    gap_matrix = fill_gap(gap,zeros_matrix,array1_length,array2_length)

    score_matrix = compute_score(match,mismatch,gap,gap_matrix,array2_length,array1_length,array1,array2)

        
    # Getting the optimal value
    optimal_value = get_optimal_value(score_matrix)

       
    # Computing the postion of the optimal value, that will be used for traceback
    n_rows,n_cols= length_matrix(score_matrix)

    new_n_rows,new_n_cols = n_rows-1,n_cols-1

    # Getting the sequence and computing the final two sequenced aligned
    sequenza_finale= post_sequence(score_matrix,sequence,new_n_rows,new_n_cols,match,mismatch,gap)[::-1]

    sequenza1,sequenza2 = final_sequences(array1,array2,sequenza_finale)

    # Printing out the Two aligned sequence with also the sequence vector
    # Columns for results
    matrix_column, sequences_column = st.columns(2)

    with matrix_column:
        st.header("Score Matrix:")
        score_matrix
        st.header("Optimal Value: {}".format(optimal_value))
    with sequences_column:
        st.header("Traceback sequence:")
        st.markdown("Sequence: {}".format(sequenza_finale))
        st.markdown("First sequence: {}".format(sequenza1))
        st.markdown("Second sequence: {}".format(sequenza2))

###########################################################################