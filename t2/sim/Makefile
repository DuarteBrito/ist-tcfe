all: ngspice1.log ngspice2.log ngspice3.log ngspice4.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	sed -n '/^circ2_TAB/,/^circ2_END/{p;/^circ2_END/q}' ngspice2.log | grep -v circ2_TAB | grep -v circ2_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > circ2_tab.tex;
	ps2pdf trans3.ps
	ps2pdf trans4.ps

ngspice1.log: ../sim/circ1.net
	ngspice -b $< -o $@

ngspice2.log: ../sim/circ2.net
	ngspice -b $< -o $@

ngspice3.log: ../sim/circ3.net
	ngspice -b $< -o $@

ngspice4.log: ../sim/circ4.net
	ngspice -b $< -o $@

clean:
	@rm -f *.log *.net *.tex *.pdf *.ps *~

.PHONY: all clean
