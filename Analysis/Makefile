
CXX := mpicxx
CXXFLAGS := -O3
CPPFLAGS := -DUSE_MPI

srcs := bici.cc header/tinyxml2.cc
hdrs := $(wildcard header/*.h)
exe := bici

$(exe): $(srcs) $(hdrs)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(srcs) -o bici

.PHONY : clean

clean:
	rm -f $(exe)
