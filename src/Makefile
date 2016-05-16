# Makefile template for shared library

CC = gcc
CFLAGS = -fPIC -Wall -Wextra -O2 -g -std=gnu99
LDFLAGS = -shared -Bsymbolic -lc
RM = rm -f
TARGET_LIB = histk.so
SRCS = histk.c
OBJS = $(SRCS:.c=.o)

.PHONY: all
all: ${TARGET_LIB}

$(TARGET_LIB): $(OBJS)
	$(CC) ${LDFLAGS} -o $@ $^

.PHONY: clean
clean:
	-${RM} ${TARGET_LIB} ${OBJS}
