

scrape: wtfgenes
	bin/scrape-gaf.js

wtfgenes:
	git clone git@github.com:evoldoers/wtfgenes.git

clean:
	rm -rf web

build: clean wtfgenes scrape
