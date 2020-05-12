
MPICXX := mpicxx
CXXFLAGS := -O3

srcs := bici.cc header/tinyxml2.cc
hdrs := $(wildcard header/*.h)

exe_serial := bici
exe_mpi := bici_mpi

$(exe_serial): $(srcs) $(hdrs)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(srcs) -o $(exe_serial)

$(exe_mpi): $(srcs) $(hdrs)
	$(MPICXX) $(CXXFLAGS) $(CPPFLAGS) -DUSE_MPI $(TARGET_ARCH) $(srcs) -o $(exe_mpi)

.PHONY : all
all: $(exe_serial) $(exe_mpi)

.PHONY : clean
clean:
	rm -f $(exe_serial) $(exe_mpi)
