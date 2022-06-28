build: pathFinder

pathFinder: pathFinder.cpp
	@echo "Compiling \"pathFinder\"." && g++ -std=c++14 -O3 -pthread pathFinder.cpp -o pathFinder

clean:
	rm -f pathFinder
