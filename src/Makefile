CC = g++
CPPFLAGS = -lgtest -std=c++17
CHECKFLAGS = -Wall -Werror -Wextra

all: test gcov_report

test: 
	@${CC} test.cpp s21_matrix_oop.cpp ${CHECKFLAGS} ${CPPFLAGS} -o test
	@./test

%.o: %.cpp
	@${CC} -std=c++17 ${CHECKFLAGS} -c $^

%.a: %.o
	@ar rcs $@ $^
	@ranlib $@

clean:
	@/bin/rm -rf *.o *.a test *.gcno *gcda report *.info test *.out

gcov_report: test.cpp s21_matrix_oop.cpp
	@$(CC) --coverage $^ ${CPPFLAGS} ${CHECKFLAGS} -o test
	@./test
	@lcov -t "test" -o test.info --no-external -c -d ./
	@genhtml -o report test.info
	@open report/index.html

style:
	python3 ../materials/linters/cpplint.py --extensions=cpp *.cpp *.h

cppcheck:
	cppcheck --enable=all --suppress=missingIncludeSystem . *.cpp

leaks: test
	CK_FORK=no leaks --atExit -- ./test

valgrind:
	valgrind --tool=memcheck --leak-check=full ./test