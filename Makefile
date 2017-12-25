AR := ar
ARFLAGS := crs
DOXYGEN := doxygen
CC := g++
LDLIBS := -lm
CFLAGS := -std=c++0x -fopenmp -O3 -march=native -ffast-math \
          -Wall -Wextra -Wno-unused-result \
		  -I /usr/include -D _USE_MATH_DEFINES
VPATH = .
OBJS_MAIN := main.o
TARGET := linopt

SRCS := $(wildcard $(addsuffix /*.cpp, $(VPATH)))
HDRS := $(wildcard $(addsuffix /*.h, $(VPATH)))
OBJS := $(filter-out $(OBJS_MAIN), $(notdir $(SRCS:.cpp=.o)))

.PHONY: all clean help static shared

all:	$(TARGET)
static: clean $(TARGET).a
shared: clean $(TARGET).so

$(TARGET):	$(OBJS_MAIN) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS) 

$(TARGET).a:	$(OBJS)
	$(AR) $(ARFLAGS) $@ $^

$(TARGET).so:	CFLAGS += -fPIC 
$(TARGET).so:	$(OBJS)
	$(CC) $(CFLAGS) -shared -Wl,-soname,$@ -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(addprefix -I , $(VPATH)) -c -MMD $<
include $(wildcard *.d)

doc:	doxygen.conf $(SRCS) $(HDRS)
	$(DOXYGEN) $<

clean:
	$(RM) -r *.o *.d *.s *.a *.so* doc/ $(TARGET)

help:
	@echo "Available commands:"
	@echo "    help - display this help message."
	@echo "    $(TARGET) - build main program."
	@echo "    $(TARGET).a - build static version of the library."
	@echo "    $(TARGET).so - build shared version of the library."
	@echo "    static - full rebuild of the static version."
	@echo "    shared - full rebuild of the shared version."
	@echo "    doc - build documentation using doxygen."
	@echo "    clean - clean build files."

