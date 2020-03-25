NTL =  -lntl -lgmp
ImQF = McCurley.o ImQF.o ZZlib.o SmithNF.o ImQFClassGroup.o
ReQF = Buchmann.o ReQF.o ReQFD.o ZZlib.o SmithNF.o ReQFClassGroup.o RGCD.o

table1: table1.o $(ImQF)
	g++ table1.o $(ImQF) $(NTL)

table2: table2.o $(ReQF)
	g++ table2.o $(ReQF) $(NTL)

fig1: fig1.o $(ImQF)
	g++ fig1.o $(ImQF) $(NTL)

fig2: fig2.o $(ReQF)
	g++ fig2.o $(ReQF) $(NTL)

