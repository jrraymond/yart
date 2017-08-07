CC=nuclang++
CFLAGS=-Wall -Wextra -O0 -g -std=c++1z

yart: yart.m.cpp
	$(CC) $(CFLAGS) yart.m.cpp -o yart

.PHONY: clean
clean:
	rm *o yart

