AR := ar
ARFLAGS := crs
DOXYGEN := doxygen
CC := g++
PYTHON := python
LIBLINOPT_PATH := ../lib
LDLIBS := -lm \
		  $(shell $(PYTHON) -c "import sys; print '-lboost_python-py%d%d'%(sys.version_info[0],sys.version_info[1])") \
		  $(shell pkg-config $(PYTHON) --libs)
CFLAGS := -std=c++0x -O3 -march=native -ffast-math -fPIC \
          -Wall -Wextra -Wno-unused-result -Wno-unused-parameter \
		  -I /usr/include \
		  $(shell pkg-config eigen3 --cflags) \
		  $(shell pkg-config $(PYTHON) --cflags) \
		  -D_USE_MATH_DEFINES \
		  -fpermissive
TARGET := pylinopt
VPATH = . $(LIBLINOPT_PATH)

SRCS := $(wildcard $(addsuffix /*.cpp, $(VPATH)))
HDRS := $(wildcard $(addsuffix /*.h, $(VPATH)))
OBJS := $(notdir $(SRCS:.cpp=.o))

.PHONY: all clean help static shared

all:	$(TARGET).so

$(TARGET).so:	$(OBJS) minieigen.a
	$(CC) $(CFLAGS) -shared -Wl,-soname,$@ -o $@ $^ $(LDLIBS)

MINIEIGEN_BUILD_DIR := minieigen/build/temp.linux-x86_64-2.7/src
MINIEIGEN_OBJS := $(shell find $(MINIEIGEN_BUILD_DIR) -type f -name '*.o')

minieigen.a: $(MINIEIGEN_OBJS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) $(addprefix -I , $(VPATH) minieigen/src minieigen/src/double-conversion minieigen) -c -MMD $<
include $(wildcard *.d)

doc:	doxygen.conf $(SRCS) $(HDRS)
	$(DOXYGEN) $<

clean:
	$(RM) -r *.o *.d *.s *.a *.so* $(LIBLINOPT_PATH)/*.o $(LIBLINOPT_PATH)/*.d doc/

help:
	@echo "Available commands:"
	@echo "    help - display this help message."
	@echo "    $(TARGET).so - build shared version of the library."
	@echo "    doc - build documentation using doxygen."
	@echo "    clean - clean build files."

