FC = gfortran
FCFLAGS = -O3 
LIBS = 

SOURCES = Spectra.f
COMBINESRC = combine.f90
COMBINE = combine.x
OBJECTS = $(SOURCES:.f=.o)
TARGET = traject.x

all: $(TARGET) $(COMBINE)

$(TARGET) : $(OBJECTS)
	$(FC) $(FCFLAGS) $(OBJECTS) $(LIBS) -o $@ 

$(COMBINE) :
	$(FC) $(FCFLAGS) $(COMBINESRC) -o $(COMBINE)

clean : 
	rm -f $(OBJECTS) $(TARGET) $(COMBINE)
 
%.o : %.f
	$(FC) $(FCFLAGS) -c $<

# DO NOT DELETE
