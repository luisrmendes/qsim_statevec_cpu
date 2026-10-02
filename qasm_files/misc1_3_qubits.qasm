// Multi-operation test on 3 qubits
include "stdgates.inc";

qubit[3] q;

h q[0];
x q[1];
y q[2];
s q[0];
t q[1];
z q[2];
cx q[1], q[0];
cz q[2], q[1];
h q[2];
