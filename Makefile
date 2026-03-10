CXX = g++
CXXFLAGS = -fopenmp

TARGET = main

SRC = \
main.cpp \
fft_seq/fft_seq.cpp \
fft_omp/fft_omp.cpp \
fft_async/fft_async.cpp \
utils/comparer.cpp

OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

run: all
	./$(TARGET)
