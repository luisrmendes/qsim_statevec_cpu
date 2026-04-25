// Grover's algorithm on 3 qubits

include "stdgates.inc";

qubit[3] q;
bit[2] c;

h q[0];
h q[1];
h q[1];
cx q[0],q[1];
h q[1];
h q[0];
x q[0];
h q[1];
x q[1];
h q[1];
cx q[0],q[1];
h q[1];
x q[0];
h q[0];
x q[1];
h q[1];
