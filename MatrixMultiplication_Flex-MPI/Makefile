CC := mpicc
PKG_CONFIG ?= pkg-config

includedir := /usr/local/include
librarydir := /usr/local/lib
UNAME := $(shell uname -m)
GLPK=-I/usr/include/ -L/usr/lib/$(UNAME)-linux-gnu


sources := mat_x_mpi.c mat_x_adm.c

binaries := mat_x_mpi mat_x_adm


objects := $(sources:.c=.o)
depends := $(sources:.c=.d)

CPPFLAGS := -I$(includedir) -MMD
CFLAGS := -I. -I$(includedir) -std=gnu99 -Wall -Wextra -Werror=uninitialized -O2 -g
.PHONY: all clean

all: $(binaries)

clean:
	$(RM) $(binaries)
	$(RM) $(objects)
	$(RM) $(depends)


# forces the creation of object files
# necessary for automatic dependency handling
$(objects): %.o: %.c

mat_x: CFLAGS +=
mat_x: LDLIBS += -lempi -lpapi $(GLPK) -lglpk -L$(librarydir) -licc -Wl,--no-undefined,-rpath-link=$(librarydir),-rpath=$(librarydir)

mat_x_mpi: CFLAGS +=
mat_x_mpi: LDLIBS += -lempi -lpapi $(GLPK) -lglpk -L$(librarydir) -licc -Wl,--no-undefined,-rpath-link=$(librarydir),-rpath=$(librarydir)

mat_x_adm: CFLAGS +=
mat_x_adm: LDLIBS += -lempi -lpapi $(GLPK) -lglpk -L$(librarydir) -licc -Wl,--no-undefined,-rpath-link=$(librarydir),-rpath=$(librarydir)

-include $(depends)
