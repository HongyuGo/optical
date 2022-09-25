TARGET = app
MODE=DEBUG

CC = gcc
SRCS  = $(shell find ./src     -type f -name *.c)
HEADS = $(shell find ./include -type f -name *.h)
OBJS = $(SRCS:.c=.o)
DEPS = Makefile.depend

INCLUDES = -I./include

ifeq ($(MODE), DEBUG)
CFLAGS = -O0 -ggdb -Wall -Wextra -Werror -Wno-unused-parameter -Wno-attributes -fno-omit-frame-pointer -fno-optimize-sibling-calls $(INCLUDES)
else
CFLAGS = -O2 -Wextra -Werror -Wno-unused-parameter -Wno-attributes $(INCLUDES)
endif



LDFLAGS = -lm


all: $(TARGET)

$(TARGET): $(OBJS) $(HEADS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

run: all
	@./$(TARGET)

.PHONY: depend clean checkmm format install
depend:
	$(CC) $(INCLUDES) -MM $(SRCS) > $(DEPS)
	@sed -i -E "s/^(.+?).o: ([^ ]+?)\1/\2\1.o: \2\1/g" $(DEPS)

clean:
	$(RM) $(OBJS) $(TARGET)

install:
	sudo apt-get update && sudo apt-get upgrade && sudo apt-get -y install gcc g++ valgrind make clang-format clangd-11 bear

format:
	clang-format -i -style=file include/*.h src/*.c

checkmm:
	valgrind --trace-children=yes \
   --leak-check=full \
   --track-origins=yes \
   --soname-synonyms=somalloc=*jemalloc* \
   --error-exitcode=1 \
   ./$(TARGET)

-include $(DEPS)
