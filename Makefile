ifdef SystemRoot
   RM = del /Q
   FixPath = $(subst /,\,$1)
   EXT = .exe
else
   ifeq ($(shell uname), Linux)
      RM = rm -f
      FixPath = $1
      EXT =
   endif
endif

all: lths$(EXT)

lths$(EXT): lths.o
	gcc -o lths$(EXT) lths.c -I. -std=c99 -O2

lths.o: lths.c
	gcc -c lths.c

clean:
	$(RM) lths.o lths$(EXT)
