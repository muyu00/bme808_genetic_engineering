import numpy as np

# task 1
def longestoverlap(s1,s2):
    l1 = len(s1)
    l2 = len(s2)
    n = min(l1,l2)
    while n>0:
        if s1[l1-n:l1] == s2[0:n]:
            return n
        n=n-1
    return 0

# task 2
# file of strings are opened
file1 = open("exercise2.txt", "r")
alllines = file1.read()
# lines are printed out
print(alllines)
# lines are spilt -- all added into a list
lines = alllines.split('\n')

# remove an item from the collection (last item)
lines.pop()
print(lines)

# determine how many lines there are to compare
n = len(lines)
print('number of sequences:', n)

# complete the following code
# determine how much each fragment overlaps with the others
# store information in matrix overlap
overlap = [[0 for i in range(n)] for j in range(n)]
for i in range(n):
      for j in range(n):
           if i != j: # no self overlaps, therefore will always be zero
               overlap[i][j] = longestoverlap(lines[i], lines[j])

print("overlap matrix:")
print(np.matrix(overlap)) # print the overlap matrix

# implementation of the greedy algorithm
# determine the largest value in the overlap
# initialize path array
path = np.zeros(n)

i, j = 0,0
max_num = 0
count = 0

while count < len(overlap):
    # initialize for the loops
    i, j = 0,0
    max_num = 0
    
    # loop through to find the largest number

    for i in range(n):
        for j in range(n):
            # check for the largest number in the overlap matrix
            if overlap[i][j] > max_num:
                max_num = overlap[i][j]
                idx_i = i
                idx_j = j
    
    # update the path to be stored in the path array
    # if we are at the last node, set the row to -1 to indicate the start
    if count == n-1:
        path[idx_i] = -1
        break

    path[idx_i] = idx_j
    print(f"the largest overlap found between sequence {idx_i} and sequence {idx_j} is: {max_num}")

    # remove all positive values in row idx_i
    # remove all positive values in column idx_j
    for k in range(len(overlap)):
        overlap[idx_i][k] = 0
        overlap[k][idx_j] = 0

    count += 1

print("\npath of the fragments:")
print(list(range(n)))


# task 4
# determine the start position (which fragment is starting the matching)
# note: path[i] == j means fragment i will overlap with j
# looping the sequences
# list comprehension
start = [i for i in range(len(path)) if i not in path][0]

# initialize variables
current = start  # current fragment we are processing
indent = 0       # indentation for printing overlaps visually

# print the first fragment of the overlap
print(f"{current} {lines[current]}")

# printing out the path array and fragments with overlaps
# while destination != -1
# format and print in order of the path
while path[int(current)] != -1:
    next_node = path[int(current)]                                                     # next fragment in path
    ov = longestoverlap(lines[int(current)], lines[int(next_node)])                    # overlap length 
    indent = indent + len(lines[int(current)]) - ov                                    # increase indentation by overlap
    print(f"{int(next_node)} {" " * indent} {lines[int(next_node)]} ({ov})")
    
    current = next_node                

