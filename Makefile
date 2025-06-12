
all: clean utils_backup.dll

# linkage :
#gcc -shared -o utils_backup.dll utils_backup.o
utils_backup.dll: utils_backup.o
	gcc -shared -fPIC -o libhubbard.dll utils_backup.c -lm

# compilation :
utils_backup.o: utils_backup.c
	gcc -Wall -g -O2 -c utils_backup.c


clean:
ifeq ($(OS),Windows_NT)
	del /Q *.o *.so 2>nul || exit 0
else
	rm -f *.o *.so
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete
endif
