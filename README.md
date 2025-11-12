   1) 

git clone https://github.com/Ossyk/PARCO-Computing-2026-242629.git

   2)  

qsub -I -q short_cpuQ -l select=1:ncpus=64:mem=64gb,walltime=01:00:00

   3)  access the project folder
   
   4) **OPTION 1** 

gcc src/matrix_generator.c src/csr_utils.c -o matrix_generator

gcc -std=c99 -lm src/spMV_seq.c src/csr_utils.c -o src/spMV_seq

gcc -fopenmp -std=c99 -lm src/spMV_parall.c src/csr_utils.c -o src/spMV_parall

gcc -fopenmp -std=c99 -lm src/simd_vs_parallfor.c src/csr_utils.c -o src/simd_vs_parallfor
   
      (do you get permission error or a  bash ^M error?)
   
   5) 

src/spMV_seq matrices/m1.csr

src/spMV_parall 4 static 1 matrices/m4.csr

src/simd_vs_parallfor  matrices/m3.csr 4 dynamic
   
   6) 

scripts/run_linux.sh matrices/m5.csr

**OPTION 2 : run directly the interactive script (provate entrambi please, cancellando la directory e facendo git clone e run directly il prossimo cmd)**

scripts/interactive_script_l.sh
   
      (do you get permission error or a  bash ^M error?)




