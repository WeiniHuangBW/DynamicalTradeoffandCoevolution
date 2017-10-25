# DynamicalTradeoffandCoevolution
This project contains the code and library files for Manuscript "Dynamical trade-offs arise from antagonistic coevolution and decrease intraspecific diversity".

# How to compile and run the code with g++
the first command will generate a file mtwist.o, the second command will generate mtwist.a, the third command will give a final executable file and the last command will run the code in g++ (please see the detailed explaination of parameters inside CoevolutionPreyPredator.c as well as in our manuscript)

1. gcc -c mtwist.c
2. ar rs mtwist.a mtwist.o
3. g++ -o runFile.out CoevolutionPreyPredator.c mtwist.a
4. ./runFile.out parameters
