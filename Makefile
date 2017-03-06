

scrape: wtfgenes
	bin/scrape-gaf.js

wtfgenes:
	git clone git@github.com:evoldoers/wtfgenes.git
