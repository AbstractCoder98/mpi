#!/bin/bash
for i in {1..10}
do
for n in 12 40 80 200 400
do
./test $n >> resuljacobiS.doc

echo acabado
done
done
