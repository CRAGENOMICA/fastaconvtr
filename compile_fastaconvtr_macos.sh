rm -rf build
export CC=/usr/bin/clang
export CXX=/usr/bin/clang++

#To compile fastaconvtr
sh ./build.sh
cp ./build/fastaconvtr ./bin
#gcc ./sources/*.c -lm -o ./bin/fastaconvtr -Wall -O3 -g -O0 -lz
