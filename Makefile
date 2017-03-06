

wtfgenes:
	git clone https://github.com/evoldoers/wtfgenes.git

web: wtfgenes wtfgenes/web
	cp -rf wtfgenes/web web

lib: wtfgenes wtfgenes/lib
	cp -rf wtfgenes/lib lib
