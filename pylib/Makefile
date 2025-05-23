AR := ar
ARFLAGS := crs
DOXYGEN := doxygen
CC := g++
PYTHON := python3
LIBLINOPT_PATH := ../lib
LDLIBS := -lm \
		  -L ${LIBLINOPT_PATH} -llinopt
CFLAGS := -std=c++17 -fopenmp -O3 -march=native -ffast-math -fPIC \
          -Wall -Wextra \
          -I /usr/include -I ${LIBLINOPT_PATH} \
          -isystem /usr/include/eigen3 \
		  $(shell $(PYTHON) -m pybind11 --includes) \
          -D_USE_MATH_DEFINES \
          -fpermissive
TARGET := linopt
VPATH = .

SRCS := $(wildcard $(addsuffix /*.cpp, $(VPATH)))
HDRS := $(wildcard $(addsuffix /*.h, $(VPATH)))
OBJS := $(notdir $(SRCS:.cpp=.o))

.PHONY: all clean doc help

all: $(TARGET).so

$(LIBLINOPT_PATH)/liblinopt.a:
	$(MAKE) -C $(LIBLINOPT_PATH)

$(TARGET).so: $(LIBLINOPT_PATH)/liblinopt.a $(OBJS)
	$(CC) $(CFLAGS) -shared -Wl,-soname,$@ -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(addprefix -I , $(VPATH)) -c -MMD $<
include $(wildcard *.d)

clean:
	$(RM) -r *.o *.d *.s *.a *.so*

help:
	@echo "Available commands:"
	@echo "    help - display this help message."
	@echo "    $(TARGET).so - build shared version of the library."
	@echo "    clean - clean build files."

