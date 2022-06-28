build: pathFinder pathFinder_indvid

pathFinder: pathFinder.cpp
	@echo "Compiling \"pathFinder\"." && g++ -std=c++14 -O3 -pthread pathFinder.cpp -o pathFinder

pathFinder_individ: pathFinder.cpp
	@echo "Compiling \"pathFinder_indvid\"." && g++ -std=c++14 -O3 -pthread pathFinder_indvid.cpp -o pathFinder_indvid

clean:
	rm -f pathFinder pathFinder_indvid
