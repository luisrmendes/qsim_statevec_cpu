// Results: 
//   Qubit 1 -> 100%
//   Qubit 2 -> 0%
//   Qubit 3 -> 0%
//   Qubit 4 -> 0%
//   Qubit 5 -> 0%

qubit[5] q;

z q[0];
h q[1];
x q[2];
y q[3];
z q[4];
x q[3];
x q[2];
h q[1];
x q[0];