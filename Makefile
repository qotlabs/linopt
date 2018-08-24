CC := g++
LIBLINOPT_PATH := ./lib
LDLIBS := -L ${LIBLINOPT_PATH} -lm -llinopt
CFLAGS := -std=c++0x -fopenmp -O3 -march=native -ffast-math \
          -Wall -Wextra -Wno-unused-result -Wno-int-in-bool-context \
		  -I /usr/include -I /usr/include/eigen3 -I . -I ${LIBLINOPT_PATH} \
		  -D_USE_MATH_DEFINES
VPATH = .
TARGET := example

SRCS := $(wildcard $(addsuffix /*.cpp, $(VPATH)))
HDRS := $(wildcard $(addsuffix /*.h, $(VPATH)))
OBJS := $(notdir $(SRCS:.cpp=.o))

.PHONY: all clean

all:	$(TARGET)

$(LIBLINOPT_PATH)/liblinopt.a:
	$(MAKE) -C $(LIBLINOPT_PATH)

$(TARGET):	$(LIBLINOPT_PATH)/liblinopt.a $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS) 

%.o: %.cpp
	$(CC) $(CFLAGS) -c -MD $<
include $(wildcard *.d)

clean:
	$(RM) -r *.o *.d
