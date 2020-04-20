#
#  Makefile for querytopo program
#
CFLAGS = -lm
SRC = querytopo.cpp 
OBJ = querytopo.o
ULIBS = /home/sonntag/Libcpp/libjohn2.a


/home/sonntag/bin/querytopo : querytopo
	mv querytopo /home/sonntag/bin/querytopo

querytopo: $(OBJ) $(ULIBS)
	g++ $(CFLAGS) -L/home/sonntag/Libcpp -o querytopo $(OBJ) -ljohn2
	
querytopo.o: querytopo.cpp
	g++ -c querytopo.cpp

$(ULIBS): FORCE
	cd /home/sonntag/Libcpp; $(MAKE)

FORCE:
