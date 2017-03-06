

scrape: wtfgenes lib
	bin/scrape-gaf.js

wtfgenes:
	git clone git@github.com:evoldoers/wtfgenes.git

lib: wtfgenes wtfgenes/lib
	cp -rf wtfgenes/lib lib
