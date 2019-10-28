import sys
from math import log, exp
from compsci260lib import get_fasta_dict

def posterior_decoding(input_file, f_hmm_file):
    # read the state names and number
    states = f_hmm_file.readline().split()
    K = len(states)
    
    # read the initial probabilities
    probs = f_hmm_file.readline().split()
    initial_probs = [float(prob) for prob in probs]
    
    # read the transition matrix
    transitions = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(trans_prob) for trans_prob in matrix_row_arry]
        transitions[i] = matrix_row
    
    # read the emission symbols
    emission_symbols = f_hmm_file.readline().split()
    
    # read the emission probability matrix
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(emit_prob) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row
    
    f_hmm_file.close()
    
    seq_dict = get_fasta_dict(input_file)
    emit_str = seq_dict.values()[0]  # there's only 1
    
    print "Done reading sequence of length ", len(emit_str)
    for a in range(len(emit_probs)):
        for c in range(len(transitions[0])):
            transitions[a][c] = log(transitions[a][c])
        for b in range(len(emit_probs[0])):
            emit_probs[a][b] = log(emit_probs[a][b])
    
    # Run the forward algorithm
    forward = run_forward(states, initial_probs, transitions, emission_symbols, emit_probs, emit_str)
    
    # Run the backward algorithm
    backward = run_backward(states, initial_probs, transitions, emission_symbols, emit_probs, emit_str)
    
    # Calculate the posterior probabilities
    # Initializing the posterior 2D matrices
    posterior = [[float(0) for _ in range(K)] for _ in range(len(emit_str))]
    for i in range(len(emit_str)):
        # Did not normalize the probabilities (i.e., did not divide by P(X)),
        # because we will only use these probabilities to compare
        # posterior[i][0] versus posterior[i][1]
        for k in range(K):
            posterior[i][k] = forward[i][k] + backward[i][k]
    
    # Print out the decoded results
    print "start\tstop\tstate"
    total_count = 0
    state2_count = 0
    curr_state = 0
    start = 0
    for i in range(len(posterior)):
        max = float('-inf')
        max_state = 0
        for k in range(K):
            if(posterior[i][k] > max):
                max = posterior[i][k]
                max_state = k
        if(max_state != curr_state or i == len(posterior) - 1):
            total_count += 1
            if curr_state == 1:
                state2_count += 1
            if i == len(posterior) - 1:
                print start + 1, "\t", i + 1, "\tstate", curr_state + 1
            else:
                print start + 1, "\t", i, "\tstate", curr_state + 1
            start = i
        curr_state = max_state
    
    print "Total number of regions reported:", total_count
    print "Total number of state 2 regions:", state2_count
    
def run_forward(states, initial_probs, transitions, 
    emission_symbols, emit_probs, emit_str):
    """Calculates the forward probability matrix"""
    K = len(states)
    L = len(emit_str)
    forward = [[float(0) for _ in range(K)] for _ in range(L)]

    # Initialize
    emit_index = get_emit_index(emit_str[0], emission_symbols)
    for k in range(K):
        forward[0][k] = log(initial_probs[k]) + emit_probs[k][emit_index]
    
    # Iterate
    for i in range(1, L): #  the index
        emit_index = get_emit_index(emit_str[i].upper(), emission_symbols)
        
        # Compute the forward probabilities for the states
        for l in range(K): # the state
            logp = transitions[0][l] + forward[i - 1][0]
            esum = 1
            for k in range(1, K):
                esum += exp((transitions[k][l] + forward[i - 1][k]) - logp)
            forward[i][l] = emit_probs[l][emit_index] + log(esum) + logp

    return forward        
        
def run_backward(states, initial_probs, transitions, 
    emission_symbols, emit_probs, emit_str):
    """Calculates the backward probability matrix"""
    K = len(states)
    L = len(emit_str)
    backward = [[float(0) for _ in range(K)] for _ in range(L)]
    
    # Initialize
    for k in range(K):
        backward[L-1][k] = log(1)  # which is zero, but just to be explicit...
        
    # Iterate
    for i in range(L-2, -1, -1):
        emit_index = get_emit_index(emit_str[i+1].upper(), emission_symbols)
        # Compute the backward probabilities for the states
        for l in range(K): # the state
            logp = transitions[0][l] + backward[i + 1][0]
            esum = 1
            for k in range(1, K):
                esum += exp((transitions[k][l] + backward[i + 1][k]) - logp)
            backward[i][l] = emit_probs[l][emit_index] + log(esum) + logp
    
    return backward        
        
def get_emit_index(input_val, alphabet):
    """does a linear search to find the index of a character"""
    for i in range(len(alphabet)):
        if alphabet[i] == input_val:
            return i
    sys.exit("Could not find character " + input_val)


if __name__ == '__main__':
    hmm_file = "../HMM.methanococcus.txt"
#     input_file = "../artificial.genome.fasta"
    input_file = "../bacterial.genome.fasta"
    
    f_in_file = open(input_file)
    f_hmm_file = open(hmm_file)
    
    if f_in_file is None:
        sys.exit("Can't open HMM file: " + hmm_file)
    if f_hmm_file is None:
        sys.exit("Can't open file: " + input_file)
    
    posterior_decoding(input_file, f_hmm_file)
