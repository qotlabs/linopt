# location of the Python header files
 
PYTHON_VERSION = 3.5
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
 
# location of the Boost Python include files and library
 
# BOOST_INC = /usr/include
# BOOST_LIB = /usr/lib

STATES_LIB = ~/c++/linopt/lib/states.o
 
# compile mesh classes
TARGET = pystates
 
$(TARGET).so: $(TARGET).o
	g++ -shared -Wl,--export-dynamic $(TARGET).o $(STATES_LIB) -o $(TARGET).so -std=gnu++11 -lboost_python-py35 -L/usr/lib/python$(PYTHON_VERSION)/config -lpython3.5m  
 
$(TARGET).o: $(TARGET).cpp
	g++ -I$(PYTHON_INCLUDE) -std=gnu++11 -fPIC -c $(TARGET).cpp