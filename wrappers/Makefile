AR := ar
ARFLAGS := crs
DOXYGEN := doxygen
CC := g++
LIBLINOPT_PATH := ../lib
LDLIBS := -lm -lboost_python-py35 -lpython3.5m
CFLAGS := -std=c++0x -fopenmp -O3 -march=native -ffast-math -fPIC \
          -Wall -Wextra -Wno-unused-result -Wno-int-in-bool-context \
		  -I /usr/include -I /usr/include/eigen3 -I/usr/include/python3.5 \
		  -D_USE_MATH_DEFINES \
		  -fpermissive
TARGET := pylinopt
VPATH = . $(LIBLINOPT_PATH)

SRCS := $(wildcard $(addsuffix /*.cpp, $(VPATH)))
HDRS := $(wildcard $(addsuffix /*.h, $(VPATH)))
OBJS := $(notdir $(SRCS:.cpp=.o))

.PHONY: all clean help static shared

all:	$(TARGET).so
static: clean $(TARGET).a
shared: clean $(TARGET).so

$(TARGET).a:	$(OBJS)
	$(AR) $(ARFLAGS) $@ $^

$(TARGET).so:	$(OBJS)
	$(CC) $(CFLAGS) -shared -Wl,-soname,$@ -o $@ $^ $(LDLIBS)
	
%.o: %.cpp
	$(CC) $(CFLAGS) $(addprefix -I , $(VPATH)) -c -MMD $<
include $(wildcard *.d)

doc:	doxygen.conf $(SRCS) $(HDRS)
	$(DOXYGEN) $<

clean:
	$(RM) -r *.o *.d *.s *.a *.so* doc/

help:
	@echo "Available commands:"
	@echo "    help - display this help message."
	@echo "    $(TARGET).a - build static version of the library."
	@echo "    $(TARGET).so - build shared version of the library."
	@echo "    static - full rebuild of the static version."
	@echo "    shared - full rebuild of the shared version."
	@echo "    doc - build documentation using doxygen."
	@echo "    clean - clean build files."

