import sys
from math import log, exp
from compsci260lib import get_fasta_dict

def posterior_decoding(input_file, f_hmm_file):
    # read the state names and number
    states = f_hmm_file.readline().split()
    K = len(states)
    
    # read the initial probabilities
    probs = f_hmm_file.readline().split()
    initial_probs = [log(float(prob)) for prob in probs]
        
    # read the transition matrix
    transitions = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [log(float(trans_prob)) for trans_prob in matrix_row_arry]
        transitions[i] = matrix_row
        
    # read the emitted symbols
    emitted_symbols = f_hmm_file.readline().split()
    
    # read the emission probability matrix
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [log(float(emit_prob)) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row
    
    f_hmm_file.close()
    
    seq_dict = get_fasta_dict(input_file)
    emit_str = seq_dict.values()[0]  # there's only 1
    
    print "Done reading sequence of length ", len(emit_str)
    
    # Run the forward algorithm
    forward = run_forward(states, initial_probs, transitions, emitted_symbols, emit_probs, emit_str)
    
    # Run the backward algorithm
    backward = run_backward(states, initial_probs, transitions, emitted_symbols, emit_probs, emit_str)
    
    # Calculate the posterior probabilities
    # Initializing the posterior 2D matrices
    posterior = [[float(0) for _ in range(K)] for _ in range(len(emit_str))]
    for i in range(len(emit_str)):
        # Did not normalize the probabilities (i.e., did not divide by P(X)),
        # because we will only use these probabilities to compare
        # posterior[i][0] versus posterior[i][1]
        for k in range(K):
            posterior[i][k] = forward[i][k] + backward[i][k]
    
    # Print the decoded results
    printout = []
    
    for i in range(len(posterior)):
        max_at_index = float('-inf')
        max_k = float('-inf')
        for j in range(len(posterior[i])):
            if posterior[i][j] > max_at_index:
                max_at_index = posterior[i][j]
                max_k = j
        
        printout.append(max_k)
    
    # You could print it out, but it's very long:
    # print "".join([str(x) for x in printout])
    
    prev_state = -1
    prev_state_index = - 1
    
    print
    print "start\tstop\tstate"
    
    for i in range(len(printout)):
        if i == 0:
            prev_state = printout[i]
            prev_state_index = 0
        else:
            if prev_state != printout[i]:
                print str(prev_state_index+1) + "\t" + str(i) + "\t" + "state " + str(prev_state+1) 
                prev_state = printout[i]
                prev_state_index = i 
    
    # print the last sequence here!
    if prev_state_index < len(printout) - 1:
        print str(prev_state_index+1) + "\t" + str(len(printout)) + "\t" + "state " + str(prev_state+1)


def run_forward(states, initial_probs, transitions, 
    emitted_symbols, emit_probs, emit_str):
    """Calculates the forward probability matrix"""
    
    K = len(states)
    L = len(emit_str)
    
    forward = [[float(0) for _ in range(K)] for _ in range(L)]
    
    # Initialize
    emit_index = get_emit_index(emit_str[0], emitted_symbols)
    for k in range(K):
        forward[0][k] = initial_probs[k] + emit_probs[k][emit_index]   # we already log-transformed all our probabilities
    
    # Iterate
    for i in range(1, L):
        emit_index = get_emit_index(emit_str[i].upper(), emitted_symbols)
        # Compute the forward probabilities for the states
        for k in range(K):
            forward[i][k] = emit_probs[k][emit_index]
            # transitions for previous position
            forward[i][k] += logsumexp(
                [forward[i-1][l] + transitions[l][k] for l in range(K)]
                )
    return forward


def run_backward(states, initial_probs, transitions, 
    emitted_symbols, emit_probs, emit_str):
    """Calculates the backward probability matrix"""
    
    K = len(states)
    L = len(emit_str)
    
    backward = [[float(0) for _ in range(K)] for _ in range(L)]
    
    # Initialize
    for k in range(K):
        backward[L-1][k] = log(1)  # which is zero, but just to be explicit...
    
    # Iterate
    for i in range(L-2, -1, -1):
        emit_index = get_emit_index(emit_str[i+1].upper(), emitted_symbols)
        # Compute the backward probabilities for the states
        for k in range(K):
            # transitions to the next position
            backward[i][k] = logsumexp(
                [transitions[k][l] + emit_probs[l][emit_index] + backward[i+1][l] for l in range(K)]
                )
    return backward


def get_emit_index(input_val, alphabet):
    """does a linear search to find the index of a character"""
    for i in range(len(alphabet)):
        if alphabet[i] == input_val:
            return i
    sys.exit("Could not find character " + input_val)


def logsumexp(logs):
    """logs is a list of float values. This function calculates log(\sum_i exp(logs_i))"""
    
    max_logs = max(logs)
    
    # for numerical stability, subtract max_logs from every number in logs 
    sum_exp_logs = sum([exp(logsi - max_logs) for logsi in logs])
    
    return log(sum_exp_logs) + max_logs


if __name__ == '__main__':
    hmm_file = "../HMM.methanococcus.txt"
    input_file = "../artificial.genome.fasta"
    input_file = "../bacterial.genome.fasta"

    f_hmm_file = open(hmm_file)
    
    if f_hmm_file is None:
        sys.exit("Can't open file: " + input_file)

    posterior_decoding(input_file, f_hmm_file)
    
