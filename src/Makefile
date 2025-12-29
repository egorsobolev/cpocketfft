
OBJ = pocketfft_wrapper.o

all: libcpocketfft.so libcpocketfft.a

libcpocketfft.so: $(OBJ)
	gcc -O0 -shared -o $@ $^

libcpocketfft.a: $(OBJ)
	ar rcs $@ $^

%.o:	%.cpp
	gcc -O3 -fPIC -I.. -c -o $@ $<

.PHONY: clean
clean:
	rm -f libcpocketfft.so
	rm -f libcpocketfft.a
	rm -f *.o
