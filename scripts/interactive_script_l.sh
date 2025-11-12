#!/bin/bash

# Get the script's directory and go to project root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR/.."

echo "Welcome to 242629's script!"
echo



while true; do
    echo "=== Main Menu ==="
    echo "1) Create your matrix"
    echo "2) Run a sequential code"
    echo "3) Run a parallel code"
    echo "4) Run default benchmark code"
    echo "5) Run personalized benchmark code"
    echo "6) See available matrices"    
    echo "7) Exit"

    echo -n "Enter your choice (1-7): "
    
    read choice
    echo
    
    case $choice in
        1)
            echo "Creating you square CSR matrix..."
            echo "Enter the rows and cols number:"
	    read size
	    echo "Enter sparsity:"
 	    read sparsity 
	    echo "Enter matrix's name (other than m1-m5 and already chosen names):"
	    read name

	    gcc src/matrix_generator.c src/csr_utils.c -o src/matrix_generator

	    src/matrix_generator $size $size $sparsity matrices/$name.csr
            echo
            ;;
        2)
            echo "Running sequential code..."
            echo "Enter matrix name (from m1 to m5 or created)"
	    read name
	    gcc -std=c99 -lm src/spMV_seq.c src/csr_utils.c -o src/spMV_seq
	    src/spMV_seq matrices/$name.csr
            echo
            ;;
        3)
            echo "Running parallel code..."
            echo "enter nb of threads:"
		read threads
		echo "enter scheduler:"
		read scheduler
		echo "enter chunks:"
		read chunks
		echo "enter matrix name:"
		read name
	        gcc -fopenmp -std=c99 -lm src/spMV_parall.c src/csr_utils.c -o src/spMV_parall 
		src/spMV_parall $threads $scheduler $chunks matrices/$name.csr
		echo
            ;;
        
	4)    
	    echo "running default benchmark test.."
	    echo "Enter matrix name (from m1 to m5 or created)"
	    read name

	    gcc -fopenmp -std=c99 -lm src/spMV_parall.c src/csr_utils.c -o src/spMV_parall 

	    #chmod +x scripts/run_linux.sh
	    #sed -i 's/\r$//' scripts/run_linux.sh
	    scripts/run_linux.sh matrices/$name.csr
	    echo
            ;;
       
	 5)echo "Running personalised benchmarking test.."
		echo "enter nb of threads:"
		read threads
		echo "enter scheduler:"
		read scheduler
		echo "enter chunks:"
		read chunks
		echo "enter matrix name:"
		read matrix

	        gcc -fopenmp -std=c99 -lm src/spMV_parall.c src/csr_utils.c -o src/spMV_parall 
		scripts/run_linux_personalized.sh $threads $scheduler $chunks matrices/$name.csr
		echo
	    ;;	
	6) echo "Available matrices :"
		cd matrices
		ls
		echo
		cd ..
	;;


	7)  echo "GoodBYE!"
	    exit 0
	    ;;
	
	*)
            echo "Invalid choice. Please try again."
            echo
            ;;
    esac
done
