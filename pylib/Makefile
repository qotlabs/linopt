AR := ar
ARFLAGS := crs
DOXYGEN := doxygen
CC := g++
PYTHON := python3
LIBLINOPT_PATH := ../lib
LDLIBS := -lm \
		  $(shell pkg-config $(PYTHON) --libs) \
		  -L ${LIBLINOPT_PATH} -llinopt
CFLAGS := -std=c++11 -fopenmp -O3 -march=native -ffast-math -fPIC \
          -Wall -Wextra -Wno-unused-parameter -Wno-int-in-bool-context\
          -I /usr/include -I ${LIBLINOPT_PATH} \
          $(shell pkg-config eigen3 --cflags) \
          $(shell pkg-config $(PYTHON) --cflags) \
          -D_USE_MATH_DEFINES \
          -fpermissive
TARGET := linopt
VPATH = .

SRCS := $(wildcard $(addsuffix /*.cpp, $(VPATH)))
HDRS := $(wildcard $(addsuffix /*.h, $(VPATH)))
OBJS := $(notdir $(SRCS:.cpp=.o))

.PHONY: all clean doc help

all:	$(TARGET).so

$(LIBLINOPT_PATH)/liblinopt.a:
	$(MAKE) -C $(LIBLINOPT_PATH)

$(TARGET).so:	$(LIBLINOPT_PATH)/liblinopt.a $(OBJS)
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
	@echo "    $(TARGET).so - build shared version of the library."
	@echo "    doc - build documentation using doxygen."
	@echo "    clean - clean build files."

