all: BEHR

functions.o: functions.c
	gcc -arch i386 -g -c functions.c  

pquad.o: pquad.f
	gfortran-mp-4.3 -c pquad.f

ptswts.o: ptswts.f
	gfortran-mp-4.3 -c ptswts.f

BEHR: BEHR.c functions.o pquad.o ptswts.o  	
	gcc -arch i386 -o BEHR BEHR.c functions.o pquad.o ptswts.o -lm 

clean:
	rm -f *~ *# *.bak *.o core* outputs.txt BEHR
	rm -f *_R.txt *_HR.txt *_C.txt *_prob.txt *_draws.txt

tar: 
	make clean 
	-tar -cvf BEHR.tar .
