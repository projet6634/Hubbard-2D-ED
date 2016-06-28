.PHONY: all clean

default: all

all:
	(cd src; $(MAKE); cd ..; cp src/main.x bin)	
clean:
	(cd src; $(MAKE) clean)
