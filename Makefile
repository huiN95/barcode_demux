build:
	cargo build --release

install:
	cp target/release/barcode_demux /usr/bin/

build-and-install:
	cargo build --release
	cp target/release/barcode_demux /usr/bin/
