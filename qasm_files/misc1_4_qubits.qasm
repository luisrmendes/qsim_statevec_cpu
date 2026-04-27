// Multi-operation test on 4 qubits
include "stdgates.inc";

qubit[4] q;

h q[0];
h q[1];
x q[2];
z q[3];
s q[2];
t q[3];
cx q[2], q[0];
cz q[3], q[1];
y q[0];
h q[3];
