CC := g++
LIBLINOPT_PATH := "lib"
LDLIBS := -L ${LIBLINOPT_PATH} -lm -llinopt
CFLAGS := -std=c++0x -fopenmp -O3 -march=native -ffast-math \
          -Wall -Wextra -Wno-unused-result -Wno-int-in-bool-context \
		  -I /usr/include -I /usr/include/eigen3 -I . -I ${LIBLINOPT_PATH} \
		  -D_USE_MATH_DEFINES
VPATH = .
TARGET := linopt

SRCS := $(wildcard $(addsuffix /*.cpp, $(VPATH)))
HDRS := $(wildcard $(addsuffix /*.h, $(VPATH)))
OBJS := $(notdir $(SRCS:.cpp=.o))

.PHONY: all clean liblinopt

all:	liblinopt $(TARGET)

$(TARGET):	$(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS) 

liblinopt:
	$(MAKE) -C $(LIBLINOPT_PATH)

%.o: %.cpp
	$(CC) $(CFLAGS) -c -MD $<
include $(wildcard *.d)

clean:
	$(RM) -r *.o *.d
	$(MAKE) -C $(LIBLINOPT_PATH) clean
