
build:
	make -C src

install:
	@mkdir -p lib
	@mkdir -p include
	mv src/libhyp.a lib/.
	cp src/*.h include/.

clean:
	make -C src clean
	@rm -vf libhyp.a
	@rm -rvf lib include
