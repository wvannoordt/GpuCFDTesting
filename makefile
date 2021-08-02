CC_HOST   := g++
CC_DEVICE := nvcc -ccbin=${CC_HOST}

SRC_FILES_HOST := $(shell find src -name *.cpp)
SRC_FILES_CUDA := $(shell find src -name *.cu)


ifndef PTL_DIR
PTL_DIR := ${PTL}
ifeq (${PTL_DIR},)
$(error "Cannot find PTL on the current system! Please install PTL at https://github.com/wvannoordt/PTL and set the environment variable PTL to the install location")
endif
endif


COMPILE_TIME_OPT :=

INCLUDE := -I. -I${PTL}/include -I./src -I/usr/local/cuda/include
LINKS   := -L${PTL}/lib -lPTL -L/usr/local/cuda/lib64 -lcudadevrt -lcudart

TARGET := gpu-tgv

HOST_OBJ   := $(addprefix obj/,$(addsuffix .o,$(basename $(notdir ${SRC_FILES_HOST}))))
DEVICE_OBJ := $(addprefix obj/,$(addsuffix .o,$(basename $(notdir ${SRC_FILES_CUDA}))))

DLINK :=
ifneq (${DEVICE_OBJ},)
DLINK := dlink
endif

HOST_FLAGS   := -x c++ -g -Wno-unknown-pragmas -std=c++11 -Werror -c ${COMPILE_TIME_OPT}
DEVICE_FLAGS := -x cu -rdc=true -Xcompiler ${COMPILE_TIME_OPT} -dc
DLINK_FLAGS  := -Xcompiler -dlink

run: main
	cd exampleRun; ./${TARGET} input.ptl
	
main: setup ${HOST_OBJ} ${DEVICE_OBJ} ${DLINK}
	${CC_HOST} ./obj/*.o -o bin/${TARGET} ${LINKS}
	cp bin/${TARGET} exampleRun

${HOST_OBJ}: obj/%.o : src/%.cpp
	${CC_HOST} ${HOST_FLAGS} ${INCLUDE} $< -o $@

${DEVICE_OBJ}: obj/%.o : src/%.cu
	${CC_DEVICE} ${DEVICE_FLAGS} ${INCLUDE} $< -o $@

dlink:
	${CC_DEVICE} ${DLINK_FLAGS} ${DEVICE_OBJ} -o obj/cu.dlink.o
	
setup:
	mkdir -p obj bin

clean:
	rm -rf obj bin
	rm -f ./exampleRun/${TARGET}