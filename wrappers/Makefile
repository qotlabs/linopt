AR := ar
ARFLAGS := crs
DOXYGEN := doxygen
CC := g++
PYTHON := python3
LIBLINOPT_PATH := ../lib
MINIEIGEN_DCONV_PATH := ./minieigen/src/double-conversion
LDLIBS := -lm \
		  $(shell pkg-config $(PYTHON) --libs) \
		  -L ${LIBLINOPT_PATH} -llinopt 
CFLAGS := -std=c++0x -fopenmp -O3 -march=native -ffast-math -fPIC \
          -Wall -Wextra -Wno-unused-parameter -Wno-int-in-bool-context\
          -I /usr/include \
          $(shell pkg-config eigen3 --cflags) \
          $(shell pkg-config $(PYTHON) --cflags) \
          -D_USE_MATH_DEFINES \
          -fpermissive
TARGET := pylinopt
VPATH = .

SRCS := $(wildcard $(addsuffix /*.cpp, $(VPATH)))
HDRS := $(wildcard $(addsuffix /*.h, $(VPATH)))
OBJS := $(notdir $(SRCS:.cpp=.o))

DCONV_SRCS := $(wildcard $(addsuffix /*.cc, $(MINIEIGEN_DCONV_PATH)))
DCONV_SRCS := $(filter-out %test.cc, $(DCONV_SRCS))
OBJS += $(notdir $(DCONV_SRCS:.cc=.o))

.PHONY: all clean doc help

all:	$(TARGET).so

$(LIBLINOPT_PATH)/liblinopt.a:	..
	$(MAKE) -C $(LIBLINOPT_PATH) static

$(TARGET).so:	$(LIBLINOPT_PATH)/liblinopt.a $(OBJS)
	$(CC) $(CFLAGS) -shared -Wl,-soname,$@ -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(addprefix -I , $(VPATH)) -c -MMD $<
include $(wildcard *.d)

%.o: $(MINIEIGEN_DCONV_PATH)/%.cc
	$(CC) $(CFLAGS) -c $<

doc:	doxygen.conf $(SRCS) $(HDRS)
	$(DOXYGEN) $<

clean:
	$(RM) -r *.o *.d *.s *.a *.so* doc/

help:
	@echo "Available commands:"
	@echo "    help - display this help message."
	@echo "    $(TARGET).so - build shared version of the library."
	@echo "    doc - build documentation using doxygen."
	@echo "    clean - clean build files."

