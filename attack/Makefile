# Makefile for Orthros Distinguisher Discovery Tools

CC = gcc
CFLAGS = -g -Wall -O3 -std=c99
LDFLAGS = -lm

# Targets
TARGETS = difflin

.PHONY: all clean help

all: $(TARGETS)
	@echo "Build complete"

# difflin - Differential-linear correlation verification tool
difflin: difflin.c orthros.c utils.c
	$(CC) $(CFLAGS) -o difflin difflin.c $(LDFLAGS)

clean:
	rm -f $(TARGETS)
	rm -f *.o *.dSYM core core.*

help:
	@echo "Orthros Distinguisher Discovery - Makefile"
	@echo "=========================================="
	@echo ""
	@echo "Targets:"
	@echo "  all     - Build difflin tool"
	@echo "  difflin - Build difflin tool"
	@echo "  clean   - Remove built files"
	@echo "  help    - Show this help"
	@echo ""
	@echo "Example:"
	@echo "  ./difflin 4 4 00000000000000000000110111100000 \\"
	@echo "            00000000000000008008000000000000 \\"
	@echo "            00000000ff0f00000000000000000000 20 8 0"
