emtrey :
	go build src/emtrey.go

clean :
	rm -f emtrey

install :
	cp emtrey /usr/bin/
