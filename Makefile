

scrape: wtfgenes web lib
	bin/scrape-gaf.js

wtfgenes:
	git clone git@github.com:evoldoers/wtfgenes.git

web: wtfgenes wtfgenes/web
	cp -rf wtfgenes/web web

lib: wtfgenes wtfgenes/lib
	cp -rf wtfgenes/lib lib
