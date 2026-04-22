# ================================================================
#  Makefile — polynomials project
#
#  Targets:
#    make          →  build ./poly_demo  (default)
#    make run      →  build & run
#    make clean    →  remove build artefacts
#    make rebuild  →  clean, then build
# ================================================================

CC      = gcc
TARGET  = poly_demo

SRCS    = polynomials.c main.c
OBJS    = $(SRCS:.c=.o)

# ── compiler flags ───────────────────────────────────────────────
CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic \
          -Wno-unused-parameter \
          -g

# ── linker flags ─────────────────────────────────────────────────
LDFLAGS = -lm

# ─────────────────────────────────────────────────────────────────
# Default target
# ─────────────────────────────────────────────────────────────────
.PHONY: all run clean rebuild

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build successful → ./$(TARGET)"

# ─────────────────────────────────────────────────────────────────
# Compile each .c to .o
# ─────────────────────────────────────────────────────────────────
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

# ─────────────────────────────────────────────────────────────────
# Convenience targets
# ─────────────────────────────────────────────────────────────────
run: all
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

rebuild: clean all
