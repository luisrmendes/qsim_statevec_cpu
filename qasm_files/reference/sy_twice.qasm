OPENQASM 3.0;
include "stdgates.inc";
gate sy a { ry(pi/2) a; }
qubit[3] q;
sy q[0];
sy q[0];
