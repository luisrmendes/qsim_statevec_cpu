OPENQASM 3.0;
include "stdgates.inc";
gate sy a { ry(pi/2) a; }
qubit[3] q;
h q[0];
s q[1];
t q[2];
sx q[0];
sy q[1];
z q[2];
cx q[0],q[1];
cz q[1],q[2];
ccx q[0],q[1],q[2];
