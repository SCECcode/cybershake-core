include ../../Compilers.mk

BIN_LIST = CreateFaultList.class

all: $(BIN_LIST)

CreateFaultList.class: CreateFaultList.java
	javac -classpath .:.. CreateFaultList.java

clean :
	-rm -f $(BIN_LIST) *.class  *~
