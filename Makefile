# wildcard
SRC = $(wildcard *.c)
LDFLAGS = $(shell pkg-config --libs raylib)
CFLAGS = $(shell pkg-config --cflags raylib)

main: $(SRC)
	gcc -g -Wall -Werror -o main $(SRC) -lm $(LDFLAGS) $(CFLAGS)