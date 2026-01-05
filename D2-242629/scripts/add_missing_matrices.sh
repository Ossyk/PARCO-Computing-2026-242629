#!/bin/bash



./src/random_matrix_generator 32000 10 1 matrices/randomN640k.mtx
./src/random_matrix_generator 64000 10 1 matrices/random1M280k.mtx
./src/random_matrix_generator 128000 10 1 matrices/random2M560k.mtx
./src/random_matrix_generator 16000 10 1 matrices/randomN320k.mtx

cd matrices

wget http://sparse-files.engr.tamu.edu/MM/Janna/Flan_1565.tar.gz
tar -xzf Flan_1565.tar.gz
mv Flan_1565/Flan_1565.mtx .
rm -rf Flan_1565/
rm Flan_1565.tar.gz

wget http://sparse-files.engr.tamu.edu/MM/Williams/webbase-1M.tar.gz
tar -xzf webbase-1M.tar.gz
mv webbase-1M/webbase-1M.mtx .
rm -rf webbase-1M/
rm webbase-1M.tar.gz

cd ..

echo "matrices are ready to be used"