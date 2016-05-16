all:
	make -C src
clean:
	make clean -C src
image:
	docker build -t histk .
test: image
	docker run -it histk ruby test.rb
