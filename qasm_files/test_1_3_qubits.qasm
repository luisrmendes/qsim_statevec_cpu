OPENQASM 3.0;
// Results: 
//   Qubit 1 -> 100%
//   Qubit 2 -> 50%
//   Qubit 3 -> 100%

include "stdgates.inc";

qubit[3] q;

x q[0];
h q[1];
y q[2];