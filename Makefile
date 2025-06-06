
all: clean utils.dll

# linkage :
utils.dll: utils_backup.o
	gcc -shared -o utils.dll utils_backup.o

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
