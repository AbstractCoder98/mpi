#!/bin/bash
for j in {1..10}
do
for n in 12 40 80 200 400
do
./ompj $n >> rompjaco.doc

echo acabado
done
done
