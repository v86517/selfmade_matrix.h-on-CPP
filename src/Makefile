CC = g++
OS = $(shell uname)
ifeq ($(OS), Linux)
	TEST_LIB = -lgtest -pthread
	LEAK = valgrind --trace-children=yes --track-fds=yes --track-origins=yes --leak-check=full --show-leak-kinds=all
else
	TEST_LIB = -lgtest
	LEAK = leaks -atExit --
endif

CFLAGS = -c -Wall -Werror -Wextra -std=c++17

TEST_FLAG = --coverage

all: s21_matrix_oop.a

s21_matrix_oop.a:
	$(CC) $(CFLAGS) -lstdc++ -std=c++17 s21_matrix_oop.cpp -o s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

gcov_s21_matrix_oop.a:
	rm -f s21_matrix_oop.a
	$(CC) $(TEST_FLAG) $(CFLAGS) -lstdc++ -std=c++17 s21_matrix_oop.cpp -o s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

clean:
	rm -rf *.o *.out *.a s21_matrix_test
	rm -rf *.gcda *.gcno report gcov_report.* gcov_report *.info

test: gcov_s21_matrix_oop.a
		$(CC) $(CFLAGS) test.cpp
	$(CC) $(TEST_FLAG) *.o s21_matrix_oop.a -o s21_matrix_test $(TEST_LIB) 
	./s21_matrix_test

gcov_report: test
	lcov -t "test" -o s21_test.info --no-external -c -d .
	genhtml -o report s21_test.info
	open ./report/index.html

style-fix:
	cp ../materials/linters/.clang-format .
	clang-format -i *.cpp
	rm -f .clang-format

style-check:
	cp ../materials/linters/.clang-format .
	clang-format -n *.cpp
	rm -f .clang-format

leaks_check: test
	$(LEAK) ./s21_matrix_test

.PHONY: all s21_matrix_oop.a gcov_s21_matrix_oop.a test gcov_report style-check style-fix clean leaks_check