CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm

symnmf: symnmf.o tester_for_c.o symnmf.h
	$(CC) -o symnmf symnmf.o tester_for_c.o $(CFLAGS)

symnmf.o: symnmf.c
	$(CC) -c symnmf.c $(CFLAGS)

tester_for_c.o: tester_for_c.c
	$(CC) -c tester_for_c.c $(CFLAGS)

clean:
	rm -f *o
	rm symnmf