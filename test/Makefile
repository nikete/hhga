.PHONY: all clean

hhga:=../bin/hhga

all: test #clean

test: $(hhga)
	prove -v t

$(hhga):
	cd .. && $(MAKE) bin/hhga

#clean:
#	echo
