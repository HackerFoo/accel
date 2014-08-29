accel: accel.c
	gcc -o accel accel.c

.PHONY: clean
clean:
	rm -f accel
