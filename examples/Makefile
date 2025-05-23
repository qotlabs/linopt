CC := g++
LIBLINOPT_PATH := ../lib
LDLIBS := -L ${LIBLINOPT_PATH} -lm -llinopt
CFLAGS := -std=c++17 -fopenmp -O3 -march=native -ffast-math \
          -Wall -Wextra -Wno-unused-result -Wno-int-in-bool-context \
		  -I /usr/include -isystem /usr/include/eigen3 -I . -I ${LIBLINOPT_PATH} \
		  -D_USE_MATH_DEFINES
VPATH = .

SRCS := $(wildcard $(addsuffix /*.cpp, $(VPATH)))
HDRS := $(wildcard $(addsuffix /*.h, $(VPATH)))
BINS := $(notdir $(SRCS:.cpp=.example))

.PHONY: all clean

all: $(LIBLINOPT_PATH)/liblinopt.a $(BINS)

$(LIBLINOPT_PATH)/liblinopt.a:
	$(MAKE) -C $(LIBLINOPT_PATH)

%.example: %.cpp
	$(CC) $(CFLAGS) -MD $< -o $@ $(LDLIBS)
include $(wildcard *.d)

clean:
	$(RM) -r *.o *.d *.example
