// Multi-operation test on 3 qubits with mixed controls
include "stdgates.inc";

qubit[3] q;

x q[0];
y q[1];
z q[2];
cx q[0], q[1];
cx q[2], q[0];
