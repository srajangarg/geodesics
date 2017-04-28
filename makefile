CC := c++
CFLAGS := -Wall -g -std=c++11 -fno-strict-aliasing
LIBS := -lX11 -lXi -lXmu -lglut -lGLU -lGL -lm
ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

SRCS := $(shell ls -1 $(ROOT_DIR)/DGP/*.cpp | sed 's/ /\\ /g') \
        $(shell ls -1 $(ROOT_DIR)/DGP/Graphics/*.cpp | sed 's/ /\\ /g') \
        $(shell ls -1 $(ROOT_DIR)/*.cpp | sed 's/ /\\ /g')
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)

.PHONY: format clean

geodesics: $(OBJS)
	@echo Linking geodesics
	@$(CC) $(CFLAGS) -o geodesics $(OBJS) $(LIBS)

%.o : %.cpp
	@echo Compiling $<
	@$(CC) $(CFLAGS) -c $< -o $@
	@$(CC) -std=c++11 -MM -MQ '$@' $< > $*.d

-include $(DEPS)

clean:
	@$(RM) $(OBJS) $(DEPS) *~ $(MAIN)
	@echo Cleaning finished

format:
	$(shell for FILE in `git ls-files | grep -E "\.(cpp|h|hpp|c)$$" | grep -Ev "DGP"`; do clang-format-3.7 -i $$FILE; done)
