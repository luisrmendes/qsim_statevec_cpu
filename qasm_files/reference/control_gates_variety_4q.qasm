OPENQASM 3.0;
include "stdgates.inc";
qubit[4] q;
x q[0];
cx q[0],q[2];
cz q[2],q[3];
ccx q[0],q[1],q[3];
