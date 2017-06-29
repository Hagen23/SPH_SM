CC = g++
CFLAGS = -std=c++11 -Wall
MATH3D = Math3D/
INCLUDES = -I Math3D/ -I SPH_SM/
LDFLAGS = -lGL -lglut -lGLU
DEBUGF = $(CFLAGS) -ggdb
SOURCES = *.cpp Math3D/*.cpp SPH_SM/*.cpp
OUTF = build/
MKDIR_P = mkdir -p

all: build

directory: $(OUTF)

build: directory $(SOURCES) $(OUTF)sph_sm

$(OUTF):
	$(MKDIR_P) $(OUTF)

$(OUTF)sph_sm: $(OUTF)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(OUTF)sph_sm $(SOURCES) $(LDFLAGS)

debug: directory $(SOURCES) $(OUTF)sph_sm_d

$(OUTF)sph_sm_d: $(OUTF)
	$(CC) $(DEBUGF) $(INCLUDES) -o $(OUTF)sph_sm_d $(SOURCES) $(LDFLAGS) 

clean:
	@[ -f $(OUTF)sph_sm ] && rm $(OUTF)sph_sm || true
	@[ -f $(OUTF)sph_sm_d ] && rm $(OUTF)sph_sm_d || true

rebuild: clean build

redebug: clean debug