#!/bin/bash
for j in {1..10}
do
for size in 12 40 80 200 400
do
mpiexec -np 4 --host master,wn1,wn2,wn3 ./jompmpi $size $size >> resultjompmpi.doc

echo acabado intento $j con vector de $size
done
sleep 2
done
